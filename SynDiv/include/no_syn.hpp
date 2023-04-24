#ifndef NO_SYN_HPP
#define NO_SYN_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include <tuple>
#include <regex>
#include <string.h>
#include <malloc.h>
#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "save.hpp"
#include "GzChunkReader.hpp"
#include "ThreadPool.hpp"
#include "cal.hpp"

using namespace std;

namespace NOSYN
{
    // 调试代码
    bool debug = false;
    // bool debug = true;
    
    // kseq.h 打开文件
    KSEQ_INIT(gzFile, gzread)



    /**
     * @brief 更改全局参数
     * 
     * @param debugTmp     是否调试代码
     * 
     * @return 0
    */
    int change_global(
        const bool & debugTmp
    )
    {
        debug = debugTmp;
        
        return 0;
    }


    // 所有索引基因组的染色体长度
    class LENGTH
    {
    private:
        // 染色体长度文件
        vector<string> lengthsVec_;
        vector<string> lengthsTitles_;

        // 所有样品的长度信息
        unordered_map<string, unordered_map<string, uint32_t> > sampleChrLenMap;  // map<sample, map<chr, chrLen> >
    public:
        LENGTH(
            vector<string> lengthsVec, 
            vector<string> lengthsTitles
        ) {
            lengthsVec_ = lengthsVec;
            lengthsTitles_ = lengthsTitles;
        }
        ~LENGTH() {}

        /**
         * @brief 打开文件
         * 
         * @return void
        **/
        void index_lengths()
        {
            for (size_t i = 0; i < lengthsVec_.size(); i++)
            {
                // 染色体长度文件
                string lengthFileName = lengthsVec_[i];

                // title
                string sampleName = "";


                // 如果含有title，用列表里的，如果不含有，根据文件名提取
                if (lengthsTitles_.size() > 0)
                {
                    sampleName = lengthsTitles_[i];
                }
                else  // 否则提取
                {
                    vector<string> lengthFileNameVec = split(lengthFileName, "/");  // 路径拆分
                    string baseName = lengthFileNameVec[lengthFileNameVec.size() - 1];  // 文件名
                    vector<string> baseNameVec = split(baseName, "."); // 按 '.' 拆分
                    sampleName = baseNameVec[0];  // 样品名为列表中第一个元素
                }

                // 初始化字典
                sampleChrLenMap[sampleName];

                // open file
                GzChunkReader GzChunkReaderClass(lengthFileName);
                // read line
                string line;
                while (GzChunkReaderClass.read_line(line))
                {
                    // 跳过空行
                    if (line.empty())
                    {
                        continue;
                    }

                    // 拆分
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    // 检查列数是否符合规定
                    if (lineVec.size() != 2)
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Error: the number of columns is not 2 -> " << line << endl;
                        exit(1);
                    }

                    // 染色体号和长度信息
                    string chromosome = lineVec[0];
                    uint32_t chrLen = stol(lineVec[1]);
                    
                    // 存储
                    sampleChrLenMap[sampleName][chromosome] = chrLen;
                }
            }
        }

        /**
         * @brief 获取染色体长度
         * 
         * @param sample 样品名
         * @param chr    染色体号
         * 
         * @return uint32_t
        **/
        uint32_t get_length(
            const string & sample, 
            const string & chr
        )
        {
            // 查找染色体长度
            // 找样品
            auto fintIter1 = sampleChrLenMap.find(sample);
            // 如果没有，报错
            if (fintIter1 == sampleChrLenMap.end())
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: '" << sample << "' not present in chromosome length file." << endl;
                exit(1);
            }
            // 找染色体
            auto findIter2 = fintIter1->second.find(chr);
            if (findIter2 == fintIter1->second.end())
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: the length of '" << chr << "' was not found -> " << sample << endl;
                exit(1);
            }

            // 染色体长度
            return findIter2->second;
        }
    };


    // 将 CAL::SYNCOOR 作为基类
    class NOSYNCOOR:public CAL::SYNCOOR
    {
    private:

    public:
        // 记录非共线性的坐标
        map<string, map<uint32_t, map<string, vector<tuple<uint32_t, uint32_t> > > > > chrStartSampleLociVecMap; // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >

        // 输出字符串
        string outTxt;

        NOSYNCOOR() {}
        NOSYNCOOR(string fileName) {
            init(fileName);
        }
        ~NOSYNCOOR() {}

        void find_no_syn_coor(
            LENGTH & lengthClass
        );
    };


    /**
     * @brief 寻找非共线性的坐标
     * 
     * @param lengthClass  样品染色体长度class
     * 
     * @return void
    **/
    void NOSYNCOOR::find_no_syn_coor(
        LENGTH & lengthClass
    )
    {   
        // 遍历共线性坐标
        for (auto iter1 = coorChrLociSampleLociMap.begin(); iter1 != coorChrLociSampleLociMap.end(); iter1++)  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
        {
            string chrTmp = iter1->first;  // 染色体号

            // 非共线性的数量，判断是否需要更改参考基因组的坐标，防止重复添加坐标
            uint32_t preNoSynNum = 0;

            uint32_t refStart;  // 参考基因组的起始
            uint32_t refEnd;  // 参考基因组的终止


            // 用于判断样品在小区域中的非共线性坐标是否被找到
            unordered_map<string, bool> sampleNoSynFindBoolMap;  // map<sample, bool>   true->找到了   false->没找到
            // 初始化每个样品的非共线性坐标是否被找到
            for (auto iter1 : sampleNameVec)  // map<sample, bool>
            {
                if (iter1 == "reference")
                {
                    continue;
                }

                sampleNoSynFindBoolMap[iter1] = false;
            }


            for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)  // map<refStart, map<sample, tuple(start, end)> >
            {
                uint32_t refStartTmp;  // 参考基因组的起始
                uint32_t refEndTmp;  // 参考基因组的终止

                // 节点的起始和终止位置
                auto findIter1 = iter2->second.find("reference");
                if (findIter1 != iter2->second.end())
                {
                    refStartTmp = get<0>(findIter1->second);
                    refEndTmp = get<1>(findIter1->second);
                }
                else
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Error: 'reference' not present in 'coorChrLociSampleLociMap' -> " << iter2->first << endl;
                    exit(1);
                }
                
                // 非共线性的数量，判断是否需要更改参考基因组的坐标，防止重复添加坐标
                uint32_t noSynNum = 0;

                // 如果上一个节点全是共线性的，更新参考基因组的坐标
                if (preNoSynNum  == 0)
                {
                    refStart = refStartTmp;
                    refEnd = refEndTmp;
                }

                for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++)  // map<sample, tuple(start, end)>
                {
                    string sampleTmp = iter3->first;  // 样品名
                    uint32_t StartTmp = get<0>(iter3->second);  // 起始
                    uint32_t EndTmp = get<1>(iter3->second);  // 终止

                    // 参考基因组的直接跳过 
                    if (sampleTmp == "reference")
                    {
                        continue;
                    }

                    // 记录为非共线性的
                    if (StartTmp == 0 || EndTmp == 0)
                    {
                        // 非共线性数量加1
                        noSynNum++;

                        // 找非共线性的位置
                        // 先初始化字典
                        if (chrStartSampleLociVecMap.find(chrTmp) == chrStartSampleLociVecMap.end())  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
                        {
                            chrStartSampleLociVecMap[chrTmp];
                        }
                        if (chrStartSampleLociVecMap[chrTmp].find(refStart) == chrStartSampleLociVecMap[chrTmp].end())  // map<refStart, map<sample, vector<tuple<start, end> > > >
                        {
                            chrStartSampleLociVecMap[chrTmp][refStart];
                        }

                        // 如果该sample坐标被找到了，跳过
                        if (sampleNoSynFindBoolMap[sampleTmp])  // map<sample, findBool>
                        {
                            continue;
                        }

                        // 临时坐标
                        uint32_t noSynStartTmp = 1;
                        uint32_t noSynEndTmp = 1;

                        auto iter2Tmp = iter2;  // 临时节点
                        
                        if (iter2 == iter1->second.begin())  // 如果是第一个节点
                        {
                            noSynStartTmp = 1;
                        }
                        else
                        {
                            iter2Tmp = iter2;
                            iter2Tmp--;  // 向上移动一位

                            // 如果上一个节点还是0，且不为第一个节点，继续前移
                            while (iter2Tmp != iter1->second.begin() && get<1>(iter2Tmp->second[sampleTmp]) == 0)
                            {
                                iter2Tmp--;  // 向上移动一位
                            }

                            // 赋值前一个坐标
                            if (get<1>(iter2Tmp->second[sampleTmp]) != 0)
                            {
                                noSynStartTmp = get<1>(iter2Tmp->second[sampleTmp]) + 1;
                            }
                        }

                        
                        iter2Tmp = iter2;  // 初始化迭代器
                        iter2Tmp++;  // 向下移动一位

                        // 如果下一个节点还是0，且不为最后一个节点，继续后移
                        while (iter2Tmp != iter1->second.end() && get<0>(iter2Tmp->second[sampleTmp]) == 0)
                        {
                            iter2Tmp++;  // 向下移动一位
                        }

                        // 如果是最后一个节点，则为染色体长度
                        if (iter2Tmp == iter1->second.end())
                        {
                            // 染色体长度
                            noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);
                        }
                        else
                        {
                            // 下一个节点的起始坐标 -1
                            noSynEndTmp = get<0>(iter2Tmp->second[sampleTmp]) - 1;
                        }

                        // 起始必须小于终止
                        if (noSynStartTmp < noSynEndTmp)
                        {
                            // 初始化字典
                            if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                            {
                                chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                            }
                            
                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(noSynStartTmp, noSynEndTmp));  // 记录非共线性坐标

                            sampleNoSynFindBoolMap[sampleTmp] = true;  // 记录该样品被找到了
                        }
                    }
                    else
                    {
                        /* 判断是否为最后一个节点 */
                        auto iter2Tmp = iter2;  // 临时节点
                        iter2Tmp++;  // 向下移动一位

                        // 如果是第一个节点，判断是否从 1 开始的比对。如果和参考基因组是共线性的话，不用找染色体前半部分，
                        if (iter2 == iter1->second.begin() && StartTmp == 0)
                        {
                            // 初始化字典
                            if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                            {
                                chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                            }

                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(1, StartTmp - 1));  // 赋值
                        }
                        else if (iter2Tmp == iter1->second.end())  // 如果是最后一个节点，判断是否到达染色体结尾
                        {
                            // 染色体长度
                            uint32_t noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);

                            if (EndTmp < noSynEndTmp)
                            {
                                // 初始化字典
                                
                                if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                                {
                                    chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                                }

                                chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(EndTmp + 1, noSynEndTmp));  // 赋值
                            }
                        }
                        
                        // 更新该样品被没找到
                        sampleNoSynFindBoolMap[sampleTmp] = false;
                    }
                }
                // 更新非共线性的数量
                preNoSynNum = noSynNum;
            }
        }


        // 结果转为字符串
        for (auto iter1 : chrStartSampleLociVecMap)  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
        {
            for (auto iter2 : iter1.second)  // map<refStart, map<sample, vector<tuple<start, end> > > >
            {
                // 空的节点跳过
                if (iter2.second.size() == 0)
                {
                    continue;
                }
                
                outTxt += iter1.first + "\t" + to_string(iter2.first);

                for (auto iter3 : iter2.second)  // map<sample, vector<tuple<start, end> > >
                {
                    outTxt += "\t" + iter3.first + ":";

                    for (size_t i = 0; i < iter3.second.size(); i++)  // vector<tuple<start, end> >
                    {
                        if (i > 0)
                        {
                            outTxt += ";";
                        }
                        
                        outTxt += to_string(get<0>(iter3.second[i])) + "-" + to_string(get<1>(iter3.second[i]));
                    }

                    outTxt += "\t";
                }

                outTxt += "\n";
            }
        }
    }
}

#endif