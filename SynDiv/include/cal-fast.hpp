#ifndef CAL_HPP
#define CAL_HPP
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
#include "GzChunkReader.hpp"
#include "save.hpp"
#include "ThreadPool.hpp"

using namespace std;

namespace CAL
{
    // 调试代码
    bool debug = false;
    // bool debug = true;

    // 线程数
    int threads = 30;

    // kseq.h 打开文件
    KSEQ_INIT(gzFile, gzread)



    /**
     * @brief 更改全局参数
     * 
     * @param debugTmp     是否调试代码
     * @param threadsTmp   线程数
     * 
     * @return 0
    */
    int change_global(
        const bool & debugTmp, 
        const int & threadsTmp
    )
    {
        debug = debugTmp;
        threads = threadsTmp;

        // 如果调试代码，线程数设置为1
        if (debug)
        {
            threads = 1;
        }
        
        return 0;
    }


    /**
     * @brief 解析参数
     * 
     * @param alignsTitles     aligns 文件 title
     * @param alignsVec        aligns 路径Vec
     * @param syriConfigFileName   syri 输出的配置文件
     * 
     * @return make_tuple(alignsMap, sampleSampleSyriMap)       map<sampleName, alignsPath>, map<sample1, map<sample2, syriOutPath> >
    */
    tuple<map<string, string>, map<string, map<string, string> > > aligns_parameter(
        const vector<string> & alignsTitles, 
        const vector<string> & alignsVec, 
        const string & syriConfigFileName
    )
    {
        /* ************************************* aligns ************************************* */
        map<string, string> alignsMap;  // map<sampleName, alignsPath>

        uint32_t indexTmp = 0;  // 记录 title 的索引
        for (auto iter1 : alignsVec)
        {
            string alignsPathTmp = iter1;
            string ailgnsTitleTmp;

            // 如果含有 title
            if (alignsTitles.size() > 0)
            {
                ailgnsTitleTmp = alignsTitles[indexTmp];
            }
            else
            {
                vector<string> alignsPathVecTmp = split(alignsPathTmp, "/");  // 路径拆分
                ailgnsTitleTmp = alignsPathVecTmp[alignsPathVecTmp.size() - 1];  // 最后一个
                std::regex reg1(".aligns");  // 替换
                ailgnsTitleTmp = regex_replace(ailgnsTitleTmp, reg1, "");  // 'An-1.aligns' 删除 '.aligns'
            }

            // 赋值
            alignsMap[ailgnsTitleTmp] = alignsPathTmp;

            indexTmp++;  // 索引迭代
        }

        /* ************************************* syri.out ************************************* */
        map<string, map<string, string>> sampleSampleSyriMap;  // map<sample1, map<sample2, syriOutPath> >

        // open file
        GzChunkReader GzChunkReaderClass(syriConfigFileName);

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

            // 检查列出是否为3
            if (lineVec.size() != 3)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: The '" << syriConfigFileName << "' dataframe does not contain three columns." << endl;
                exit(1);
            }

            // 样品和路径信息
            const string& sample1 = lineVec[0];
            const string& sample2 = lineVec[1];
            const string& syriPath = lineVec[2];

            // 初始化
            auto it = sampleSampleSyriMap.find(sample1);
            if (it == sampleSampleSyriMap.end())
            {
                it = sampleSampleSyriMap.emplace(sample1, map<string, string>{}).first;
            }

            it->second[sample2] = syriPath;
        }
        
        return make_tuple(alignsMap, sampleSampleSyriMap);
    }



    struct CALSTRUCTURE
    {
        // sampleName
        string sampleName;

        // 记录位点共线性的数量
        map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >

        // 记录位点共线性对应的样品名和位置
        map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
    };




    // 坐标转换
    class SYNCOOR
    {
    private:
        // 共线性坐标文件  'coor' 输出结果
        string fileName_;

        // 判断 'sampleNameVec' 中是否含有所有的sample
        bool sampleNameVecBool = false;
    public:
        // 共线性坐标
        map<string, map<int, map<string, tuple<uint32_t, uint32_t> > > > coorChrLociSampleLociMap;  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >

        // 样品名
        vector<string> sampleNameVec;

        SYNCOOR() {}
        SYNCOOR(string fileName) {
            fileName_ = fileName;
        }
        ~SYNCOOR() {}

        // 初始化接口，用于继承类使用
        void init(string fileName) {
            fileName_ = fileName;
        }

        // 构建索引
        void open_coor()
        {
            // open file
            GzChunkReader GzChunkReaderClass(fileName_);

            // read line
            string line;
            while (GzChunkReaderClass.read_line(line))
            {
                // 跳过空行
                if (line.empty())
                {
                    continue;
                }

                // 制表符拆分
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // 先检查列数是不是三的倍数
                if (lineVec.size() % 3 != 0)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "'Error: The number of columns is not a multiple of three -> " << fileName_ << endl;
                    exit(1);
                }

                // 每三列记录一次
                int idxTmp = lineVec.size() / 3;
                // 临时键 
                string chrTmp = "";
                uint32_t refStartTmp = 0;
                for (size_t i = 0; i < idxTmp; i++)
                {
                    string sampleNameTmp = lineVec[3*i + 0];
                    uint32_t qryStartTmp = 0;
                    uint32_t qryEndTmp = 0;

                    // 判断是否是数字，不是一律为0
                    if (isdigit(lineVec[3*i + 1][0]))
                    {
                        qryStartTmp = stol(lineVec[3*i + 1]);
                    }
                    if (isdigit(lineVec[3*i + 2][0]))
                    {
                        qryEndTmp = stol(lineVec[3*i + 2]);
                    }
                    
                    // 染色体号
                    if (i == 0)
                    {
                        chrTmp = sampleNameTmp;
                        refStartTmp = stoi(lineVec[3*i + 1]);

                        // 初始化 map 
                        if (coorChrLociSampleLociMap.find(chrTmp) == coorChrLociSampleLociMap.end())
                        {
                            coorChrLociSampleLociMap[chrTmp];
                        }
                        if (coorChrLociSampleLociMap[chrTmp].find(refStartTmp) == coorChrLociSampleLociMap[chrTmp].end())
                        {
                            coorChrLociSampleLociMap[chrTmp][refStartTmp];
                        }
                        
                        sampleNameTmp = "reference";
                    }

                    // 记录 'sampleName'，只记录第一行的
                    if (!sampleNameVecBool)
                    {
                        sampleNameVec.push_back(sampleNameTmp);
                    }
                    
                    // 添加到总的哈希表中
                    coorChrLociSampleLociMap[chrTmp][refStartTmp][sampleNameTmp] = make_tuple(qryStartTmp, qryEndTmp);
                }

                // 只记录第一行的 'sampleName'
                sampleNameVecBool = true;
            }
        }
    };


    // 坐标转换
    class COORTRANS
    {
    private:
        // show-align 输出结果 'xx.aligns'
        string aliFileName_;

        // sampleName
        string sampleName;

        // open file
        GzChunkReader* GzChunkReaderClass_;
        bool endBool = false;  // 记录文件是否遍历完

        // 记录是否要循环接下来的行
        bool chrBool = true;  // 记录该染色体对是否符合要求，目前是染色体号是否一致，不一致跳过
        bool aliBool = true;  // 记录该alignment的起始和终止是否是想要的
        bool findChrBool = true;  // 记录要找的染色体号和当前染色体号是否一致

        // 遍历文件时更改的值
        // 染色体号
        string chr = "";
        // alignment的方向、起始和终止
        string aliRefStrand = "";
        uint32_t aliRefStart = 0;
        uint32_t aliRefEnd = 0;
        string aliQryStrand = "";
        uint32_t aliQryStart = 0;
        uint32_t aliQryEnd = 0;
        // 文件流指针所在行的比对信息
        uint32_t refStart = 0;
        string refSeq = "";
        uint32_t refEnd = 0;
        uint32_t qryStart = 0;
        string qrySeq = "";
        uint32_t qryEnd = 0;
        // 记录遍历到了 '_ref_seq' 以及对应的坐标
        uint32_t refLoci = 0;
        int idxTmp = -1;
        uint32_t qryLoci = 0;

        tuple<string, uint32_t, uint32_t> get_alignment_loc(
            string infoTmp
        );
        int next_loci();  // 下一个比对段的坐标
    public:
        COORTRANS() {}
        COORTRANS(string aliFileName) {
            aliFileName_ = aliFileName;
            GzChunkReaderClass_ = new GzChunkReader(aliFileName_, 1024 * 1024 * 10);
        }
        ~COORTRANS() {
            // // 释放 GzChunkReader
            // if (GzChunkReaderClass_ != nullptr) {  
            //     delete GzChunkReaderClass_;  // 释放指针
            //     GzChunkReaderClass_ = nullptr;  // 避免二次free
            // }
        }
        
        // 找坐标
        uint32_t find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };


    /**
     * @brief 解析 '-- BEGIN alignment ' 字段
     * 
     * @param infoTmp   删除 '-- BEGIN alignment ' 后按 ' | ' 分割的字段，'+1 278 - 1703'
     * 
     * @return tuple<strand, start, end>
    */
    tuple<string, uint32_t, uint32_t> COORTRANS::get_alignment_loc(
        string infoTmp
    )
    {
        // 获取比对方向，起始和终止的信息  ref
        vector<string> infoTmp_list = split(infoTmp, " - ");  // '+1 278 - 1703'
        vector<string> infoTmp_list_tmp = split(infoTmp_list[0], "1 ");  // '+1 278'
        string strandTmp = infoTmp_list_tmp[0];
        uint32_t startTmp;
        uint32_t endTmp;
        if (strandTmp == "+")
        {
            startTmp = stoi(infoTmp_list_tmp[1]);
            endTmp = stoi(infoTmp_list[1]);
        }
        else
        {
            startTmp = stoi(infoTmp_list[1]);
            endTmp = stoi(infoTmp_list_tmp[1]);
        }

        return make_tuple(strandTmp, startTmp, endTmp);
    }

    /**
     * @brief 下一个比对段的坐标
     * 
     * @return int  0
    **/
    int COORTRANS::next_loci()
    {
        if (!endBool)
        {
            // read line
            string line;
            if (GzChunkReaderClass_->read_line(line))
            {
                // 跳过空行
                if (line.empty())
                {
                    return 0;
                }

                // 新的染色体对
                if (line.find("-- Alignments between") != string::npos)  // 新的比对染色体    -- Alignments between chr1 and chr1
                {
                    // 拆分字符串
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    string refChr = lineVec[3];  // ref染色体的临时值
                    string qryChr = lineVec[5];  // qry染色体的临时值

                    if (refChr != qryChr)
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Warning: Skipped alignment due to chromosomal mismatch -> " << refChr << " != " << qryChr << endl;
                        chrBool = false;  // 记录染色体符合规定
                    }
                    else
                    {
                        chrBool = true;  // 记录染色体不符合规定
                    }

                    // 该判断防止丢失第一个 "-- BEGIN alignment"
                    if (chr.length() == 0)
                    {
                        findChrBool = true;  // 记录是需要的染色体
                    }
                    
                    // 比对染色体信息更新
                    chr = refChr;
                }
                // 染色体必须符合规定
                else if (line.find("-- BEGIN alignment ") != string::npos && chrBool && findChrBool)  // 新的比对坐标   -- BEGIN alignment [ +1 1078 - 68996 | +1 1 - 67961 ]
                {
                    line = strip(line.erase(0, 21), '\n');  // 删除 '-- BEGIN alignment [ ' 共21个字符
                    line = line.erase(line.size()-2, 2);  // 删除 ' ]' 共2个字符，位置在最后两个

                    vector<string> lineVec = split(line, " | ");  // 分割 '+1 278 - 1703 -1 2148 - 751'

                    // 获取比对方向，起始和终止的信息  ref
                    tie(aliRefStrand, aliRefStart, aliRefEnd) = get_alignment_loc(
                        lineVec[0]
                    );

                    // 获取比对方向，起始和终止的信息  qry
                    tie(aliQryStrand, aliQryStart, aliQryEnd) = get_alignment_loc(
                        lineVec[1]
                    );
                }
                // 只循环是数字开头和包含syn的行
                // 染色体必须符合规定，比对起始和终止必须符合规定
                else if (isdigit(line[0]) != 0 && chrBool && findChrBool && aliBool)
                {
                    /* ***************************** ref ***************************** */
                    // 拆分字符串
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    // 重置坐标 (this)
                    refSeq = lineVec[1];
                    // 临时的ref序列
                    string refSeqTmp = refSeq;
                    // 删除  '.'  计算长度
                    refSeqTmp.erase(remove(refSeqTmp.begin(), refSeqTmp.end(), '.'), refSeqTmp.end());
                    if (aliRefStrand == "+")  // 正向比对
                    {
                        refStart = stoi(lineVec[0]);
                        refEnd = refStart + refSeqTmp.size() - 1;
                        refLoci = refStart - 1;  // 刷新坐标
                    }
                    else  // 反向比对
                    {
                        refEnd = stoi(lineVec[0]);
                        refStart = refEnd - refSeqTmp.size() + 1;
                        refLoci = refEnd + 1;  // 刷新坐标
                    }

                    /* ***************************** qry ***************************** */
                    // 读取下一行
                    GzChunkReaderClass_->read_line(line);

                    // 拆分字符串
                    std::istringstream iss1(line);
                    vector<string> infoVecTmp(std::istream_iterator<std::string>{iss1}, std::istream_iterator<std::string>());
                    lineVec.assign(infoVecTmp.begin(), infoVecTmp.end());

                    // 重置坐标 (this)
                    qrySeq = lineVec[1];
                    // 临时的ref序列
                    string qrySeqTmp = qrySeq;
                    // 删除  '.'  计算长度
                    qrySeqTmp.erase(remove(qrySeqTmp.begin(), qrySeqTmp.end(), '.'), qrySeqTmp.end());
                    if (aliQryStrand == "+")  // 正向比对
                    {
                        qryStart = stoi(lineVec[0]);
                        qryEnd = qryStart + qrySeqTmp.size() - 1;
                        qryLoci = qryStart - 1;  // 刷新坐标
                    }
                    else  // 反向比对
                    {
                        qryEnd = stoi(lineVec[0]);
                        qryStart = qryEnd - qrySeqTmp.size() + 1;
                        qryLoci = qryEnd + 1;  // 刷新坐标
                    }

                    /* ***************************** 刷新坐标 ***************************** */
                    idxTmp = -1;
                }
            }
            else  // 文件遍历完
            {
                // 全部坐标归零
                // 遍历文件时更改的值
                // 染色体号
                chr.clear();
                // alignment的方向、起始和终止
                aliRefStrand.clear();
                aliRefStart = 0;
                aliRefEnd = 0;
                aliQryStrand.clear();
                aliQryStart = 0;
                aliQryEnd = 0;
                // 文件流指针所在行的比对信息
                refStart = 0;
                refSeq.clear();
                refEnd = 0;
                qryStart = 0;
                qrySeq.clear();
                qryEnd = 0;
                // 记录遍历到了 '_ref_seq' 以及对应的坐标
                refLoci = 0;
                idxTmp = -1;
                qryLoci = 0;

                endBool = true;  // 记录文件是否遍历完
            }
        }

        return 0;
    }

    /**
     * @brief 寻找loci
     * 
     * @param chr_  要找的染色体
     * @param refLoci_  要找的坐标
     * @return uint32_t  0->遍历完/没找到, >0->坐标
    **/
    uint32_t COORTRANS::find_loci(
        string chr_, 
        uint32_t refLoci_
    )
    {
        /* ***************************** 遍历完，返回 ‘0’ ***************************** */
        if (endBool)
        {
            return 0;
        }
        
        /* ***************************** 没有遍历完 ***************************** */
        // 判断染色体是不是想要的
        while (chr != chr_ && !endBool)
        {
            findChrBool = false;  // 记录不是需要的染色体
            next_loci();  // 下一行
        }
        findChrBool = true;  // 记录是需要的染色体
        // 如果遍历完，返回 '0'
        if (endBool)
        {
            return 0;
        }

        // 判断alignment终止是不是小于  'refLoci_'
        while (aliRefEnd < refLoci_ && !endBool)  // 上边已经判断了染色体是否相符，所以此处不用判断
        {
            aliBool = false;  // 记录不是要的 align
            next_loci();  // 下一行
        }
        aliBool = true;  // 记录是要的 align
        // 如果遍历完，返回 '0'
        if (endBool)
        {
            return 0;
        }

        // 文件指针指向想要的行
        while (refEnd < refLoci_ && !endBool)  // 上边已经判断了染色体是否相符，所以此处不用判断
        {
            next_loci();  // 下一行
        }
        // 如果遍历完，返回 '0'
        if (endBool)
        {
            return 0;
        }

        // 找到想要的行
        if (refStart <= refLoci_ && refLoci_ <= refEnd and refLoci < refLoci_)
        {
            for (size_t idx = idxTmp + 1; idx < refSeq.length(); idx++)
            {
                // ref
                if (refSeq[idx] != '.')
                {
                    if (aliRefStrand == "+")  // 正向
                    {
                        refLoci++;
                    }
                    else  // 反向
                    {
                        refLoci--;
                    }
                }

                // qry
                if (qrySeq[idx] != '.')
                {
                    if (aliQryStrand == "+")  // 正向
                    {
                        qryLoci++;
                    }
                    else  // 反向
                    {
                        qryLoci--;
                    }
                }

                idxTmp = idx;  // 刷新索引

                // 判断是不是要找的坐标
                if (refLoci == refLoci_)
                {
                    break;  // 跳出循环
                }
            }
        }
        else if (refStart <= refLoci_ && refLoci_ <= refEnd and refLoci == refLoci_)
        {
            return qryLoci;  // 返回坐标
        }
        else  // 没找到
        {
            return 0;
        }

        return qryLoci;
    }


    // 获取共线性坐标
    class SYRIOUT
    {
    private:
        // syri.out
        string fileName_ = "";

        // open file
        GzChunkReader* GzChunkReaderClass_;
        bool endBool = false;  // 记录文件是否遍历完

        // 记录共线性信息
        string chr = "";
        uint32_t refStart = 0;
        uint32_t refEnd = 0;
        uint32_t qryStart = 0;
        uint32_t qryEnd = 0;

        // 私有函数
        int next_loci();  // 下一个比对段的坐标
    public:
        SYRIOUT() {}
        SYRIOUT(string fileName) {
            fileName_ = fileName;
            GzChunkReaderClass_ = new GzChunkReader(fileName_, 1024 * 1024 * 10);
        }
        ~SYRIOUT() {
            // // 释放 GzChunkReader
            // if (GzChunkReaderClass_ != nullptr) {  
            //     delete GzChunkReaderClass_;  // 释放指针
            //     GzChunkReaderClass_ = nullptr;  // 避免二次free
            // }
        }

        // void open_syri();

        int find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };

    /**
     * @brief 下一个共线性坐标
     * 
     * @return 0
    **/
    int SYRIOUT::next_loci()
    {
        // 没有遍历完
        if (!endBool)
        {
            // read line
            string line;
            if (GzChunkReaderClass_->read_line(line))
            {
                // 跳过空行或者不含 SYNAL 的行
                if (line.empty() || line.find("SYNAL") == string::npos)
                {
                    return 0;
                }

                // 拆分
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // 检查列数是否正常
                if (lineVec.size() != 12)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "'" << fileName_ << "': does not contain 12 columns. -> " << line << endl;
                    exit(1);
                }

                // 赋值
                chr = lineVec[0];
                refStart = stoi(lineVec[1]);
                refEnd = stoi(lineVec[2]);
                qryStart = stoi(lineVec[6]);
                qryEnd = stoi(lineVec[7]);
            }
            else  // 文件遍历完
            {
                // 全部变量清空
                chr.clear();
                refStart = 0;
                refEnd = 0;
                qryStart = 0;
                qryEnd = 0;

                endBool = true;  // 记录文件遍历完
            }
        }
        return 0;
    }

    /**
     * @brief 寻找loci
     * 
     * @param chr_  要找的染色体
     * @param refLoci_  要找的坐标
     * 
     * @return int  -1->遍历完，0->没找到, 1->坐标
    **/
    int SYRIOUT::find_loci(
        string chr_, 
        uint32_t refLoci_
    )
    {
        /* ***************************** 遍历完，返回 ‘-1’ ***************************** */
        if (endBool)
        {
            return -1;
        }
        
        /* ***************************** 没有遍历完 ***************************** */
        // 判断染色体是不是想要的
        while (chr != chr_ && !endBool)
        {
            next_loci();  // 下一行
        }
        // 如果遍历完，返回 '-1'
        if (endBool)
        {
            return -1;
        }

        // 判断refEnd是不是小于  'refLoci_'
        while (refEnd < refLoci_ && !endBool)  // 上边已经判断了染色体是否相符，所以此处不用判断
        {
            next_loci();  // 下一行
        }
        // 如果遍历完，返回 '-1'
        if (endBool)
        {
            return -1;
        }

        // 在共线性区间
        if (refStart <= refLoci_ && refLoci_ <= refEnd)
        {
            return 1;
        }
        else  // 不在区间内
        {
            return 0;
        }

        return 0;
    }



    /**
     * @brief 构建参考基因组索引
     * 
     * @param referenceFileName -> refgenome
     * 
     * @return refLenMap  map<chr, length>
    **/
    map<string, uint32_t> build_reference_index(
        const string & referenceFileName
    )
    {
        map<string, uint32_t> refLenMap;  // 参考基因组长度信息

        // 输入文件流
        gzFile gzfp = gzopen(referenceFileName.c_str(), "rb");

        // 打开文件
        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "'"
                 << referenceFileName 
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        }
        else
        {
            kseq_t *ks;
            ks = kseq_init(gzfp);
        
            while( kseq_read(ks) >= 0 )
            {
                string chromosome = ks->name.s;
                uint32_t chrLen = ks->seq.l;
                // string sequence = ks->seq.s;

                refLenMap[chromosome] = chrLen;
            }

            // 释放内存，关闭文件
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        return refLenMap;
    }




    /* ************************************** calculate syntenic diversity ************************************** */

    /**
     * @brief 计算
     * 
     * @param sampleName               样品名
     * @param refLenMap                参考基因组长度
     * @param SynCoorTmp               共线性坐标
     * @param sampleSampleSyriMap      syri.out 输出路径字典
     * @param alignsMap                show-aligns 输出路径字典
     * 
     * @return CALSTRUCTURE
    **/
    CALSTRUCTURE calculate_run(
        const string sampleName, 
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap
    )
    {
        // 保存结果
        CALSTRUCTURE CALSTRUCTURETMP;

        // 样品名
        CALSTRUCTURETMP.sampleName = sampleName;

        /* ************************************************ Find the index of sampleName ************************************************ */
        // 找 sampleName 后一个索引
        vector<string>::const_iterator findResult = find(SynCoorTmp.sampleNameVec.begin(), SynCoorTmp.sampleNameVec.end(), sampleName);
        int idxTmp = distance(SynCoorTmp.sampleNameVec.begin(), findResult) + 1;

        // 如果是最后一个样本，直接跳过
        if (idxTmp > SynCoorTmp.sampleNameVec.size() - 1)
        {
            return CALSTRUCTURETMP;
        }
        
        /* ************************************************ Syntenic Coordinate ************************************************ */
        // 共线性坐标  class的map
        map<string, SYRIOUT> SyriOutMap;  // map<sampleName, SYRIOUT>
        if (sampleName != "reference")  // 参考基因组与别的比较不需要，如果到结尾了，也跳过
        {
            // 找sample的迭代器
            map<string, map<string, string> >::const_iterator findIter1 = sampleSampleSyriMap.find(sampleName);
            if (findIter1 == sampleSampleSyriMap.end()) // 如果没有提交该样品对应的 syri.out 文件，报错
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: '" << sampleName << "' is not present in sampleSampleSyriMap." << endl;
                exit(1);
            }

            // 从索引后样本开始判断
            for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
            {
                string sampleName2 = SynCoorTmp.sampleNameVec[i];

                if (sampleName2 == "reference")  // 参考基因组没有 syri.out
                {
                    continue;
                }

                // 找sample2的迭代器
                map<string, string>::const_iterator findIter2 = findIter1->second.find(sampleName2);
                if (findIter2 == findIter1->second.end()) // 如果没有提交该样品对应的 syri.out 文件，报错
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Both '" << sampleName << "' and '" << sampleName2 << "' are not present in sampleSampleSyriMap." << endl;
                    exit(1);
                }
                string syriOutPath = findIter2->second;

                // 临时结构体
                SyriOutMap[sampleName2] = SYRIOUT(syriOutPath);
            }
        }

        /* ************************************************ Coordinate transformation ************************************************ */
        // 坐标转换 class的map
        map<string, COORTRANS> CoorTransMap;  // map<sampleName, COORTRANS>
        for (auto iter1 : SynCoorTmp.sampleNameVec)  // vector<sampleName>
        {
            if (iter1 == "reference")  // 参考基因组与别的比较不需要转换
            {
                continue;
            }

            // 查找 aligns 路径
            map<string, string>::const_iterator findIter1 = alignsMap.find(iter1);
            if (findIter1 == alignsMap.end()) // 如果没有提交该样品对应的 aligns 文件，报错
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: '" << iter1 << "' cannot be found in alignsMap." << endl;
                exit(1);
            }
            
            // 临时结构体
            CoorTransMap[iter1] = COORTRANS(findIter1->second);
        }

        /* ************************************************ Calculate syntenic diversity ************************************************ */
        // 记录每条染色体最后一个坐标，用于判断是否到达染色体末尾
        map<string, uint32_t> refChrMapTmp;  // map<chr, length>

        // 遍历总的坐标哈希表
        for (const auto& iter1 : SynCoorTmp.coorChrLociSampleLociMap)  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
        {
            // 染色体号
            string chrTmp = iter1.first;

            // 初始化字典
            auto& sampleSynOutMap = CALSTRUCTURETMP.sampleSynOutMap[chrTmp];
            auto& sampleSynOutMapTmp =  CALSTRUCTURETMP.sampleSynOutMapTmp[chrTmp];

            // if (CALSTRUCTURETMP.sampleSynOutMap.find(chrTmp) == CALSTRUCTURETMP.sampleSynOutMap.end())
            // {
            //     CALSTRUCTURETMP.sampleSynOutMap[chrTmp];

            //     // debug
            //     if (debug)
            //     {
            //         CALSTRUCTURETMP.sampleSynOutMapTmp[chrTmp];
            //     }
            // }

            // 记录上一个坐标
            uint32_t refEndTmp = 0;

            // 记录上一个 refLoci 的共线性得分，如果一样，则不存储。减少内存消耗
            uint32_t preSynNum = 0;
            map<string, uint32_t> preSynNumMap;
            for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
            {
                string sampleName2 = SynCoorTmp.sampleNameVec[i];
                preSynNumMap[sampleName2] = 0;
            }

            for (const auto& iter2 : iter1.second)  // map<refStart, map<sample, tuple(start, end)> >
            {
                // 该行的起始和终止
                uint32_t refStart = get<0>(iter2.second.at("reference"));
                uint32_t refEnd = get<1>(iter2.second.at("reference"));

                // 记录共线性的终止坐标
                refChrMapTmp[chrTmp] = refEnd;

                // 临时坐标map，所有sample的，减少查询的时间
                // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
                map<string, map<int, map<string, tuple<uint32_t, uint32_t> > > >::const_iterator findIter1 = SynCoorTmp.coorChrLociSampleLociMap.find(chrTmp);
                if (findIter1 == SynCoorTmp.coorChrLociSampleLociMap.end())
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The syntenic coordinate file does not include  '" << chrTmp << "'." << endl;
                    exit(1);
                }
                map<int, map<string, tuple<uint32_t, uint32_t> > >::const_iterator findIter2 = findIter1->second.find(refStart);
                if (findIter2 == findIter1->second.end())
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The syntenic coordinate file does not include '" << refStart << "'." << endl;
                    exit(1);
                }
                
                map<string, tuple<uint32_t, uint32_t> > SampleLociMapTmp = findIter2->second;  // map<sample, tuple<start, end> >

                // 如果和上一个坐标有间隔
                if (refStart - refEndTmp > 1)
                {
                    for (size_t refLoci = refEndTmp + 1; refLoci < refStart; refLoci++)
                    {
                        uint32_t synNumTmp = 0;  // 共线性数量，临时

                        for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
                        {
                            string sampleName2 = SynCoorTmp.sampleNameVec[i];

                            // 如果是reference与大家没坐标，那就是都没坐标
                            if (sampleName == "reference")
                            {
                                continue;
                            }
                            
                            // 如果是其它
                            uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName
                            uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName2

                            if (lociA == 0 && lociB == 0) // both del
                            {
                                synNumTmp++;  // 记为 syn

                                // debug
                                if (debug)
                                {
                                    sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                }
                                
                                // debug
                                if (debug)
                                {
                                    cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                }
                            }
                            else if (lociA != 0 && lociB == 0)  // 一个有坐标，另一个没有，one del
                            {
                                continue;
                            }
                            else if (lociA == 0 && lociB != 0)  // 一个有坐标，另一个没有，one del
                            {
                                continue;
                            }
                            else  // 都有坐标，再判断坐标的A和B是不是共线性
                            {
                                if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A和B是共线性
                                {
                                    synNumTmp++;  // 记为 syn

                                    // debug
                                    if (debug)
                                    {
                                        sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";

                                        cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                    }
                                }
                            }
                        }

                        // 如果和上一个的 synNum 不一样则添加
                        // if (synNumTmp != preSynNum || (synNumTmp == 0 && preSynNum == 0))
                        if (synNumTmp != preSynNum)
                        {
                            sampleSynOutMap[refLoci] += synNumTmp;
                            preSynNum = synNumTmp;
                        }
                    }
                }

                // 计算比对上区域的得分
                for (size_t refLoci = refStart; refLoci < refEnd + 1; refLoci++)
                {
                    uint32_t synNumTmp = 0;  // 共线性数量，临时

                    for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
                    {
                        string sampleName2 = SynCoorTmp.sampleNameVec[i];

                        // 如果在 'coor.txt' 中有坐标，代表是共线性，直接记录
                        // SampleLociMapTmp  ->  // map<sample, tuple<start, end> >
                        if (
                            get<0>(SampleLociMapTmp[sampleName]) > 0 && 
                            get<0>(SampleLociMapTmp[sampleName2]) > 0
                        )  // both syntenic
                        {
                            synNumTmp++;  // 记为 syn

                            // debug
                            if (debug)
                            {
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ";";
                            }
                        }
                        else if (
                            get<0>(SampleLociMapTmp[sampleName]) == 0 && 
                            get<0>(SampleLociMapTmp[sampleName2]) > 0
                        )  // one syntenic
                        {
                            continue;
                        }
                        else if (
                            get<0>(SampleLociMapTmp[sampleName]) > 0 && 
                            get<0>(SampleLociMapTmp[sampleName2]) == 0
                        )  // one syntenic
                        {
                            continue;
                        }
                        else  // no syntenic
                        {
                            // 如果是reference与大家没坐标，那就是都没坐标
                            if (sampleName == "reference")
                            {
                                continue;
                            }

                            // 如果是其它
                            uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName
                            uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName2

                            if (lociA == 0 && lociB == 0) // both del
                            {
                                synNumTmp++;  // 记为 syn
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";

                                // debug
                                if (debug)
                                {
                                    cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                }
                            }
                            else if (lociA != 0 && lociB == 0)  // 一个有坐标，另一个没有，one del
                            {
                                continue;
                            }
                            else if (lociA == 0 && lociB != 0)  // 一个有坐标，另一个没有，one del
                            {
                                continue;
                            }
                            else  // 都有坐标，再判断坐标对应的A和B是不是共线性
                            {
                                if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A和B是共线性
                                {
                                    synNumTmp++;  // 记为 syn

                                    // debug
                                    if (debug)
                                    {
                                        sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                        cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                    }
                                }
                            }
                        }
                    }

                    // 如果和上一个的 synNum 不一样则添加
                    // if (synNumTmp != preSynNum || (synNumTmp == 0 && preSynNum == 0))
                    if (synNumTmp != preSynNum)
                    {
                        sampleSynOutMap[refLoci] += synNumTmp;
                        preSynNum = synNumTmp;
                    }
                }

                // 更新坐标
                refEndTmp = refEnd;
            }
        }

        // 判断是否到达染色体末尾
        for (const auto& iter1 : refChrMapTmp)  // map<chr, length>
        {
            string chrTmp = iter1.first;  // 染色体
            uint32_t refPosTmp = iter1.second;  // 最后一个坐标

            // 初始化字典
            auto& sampleSynOutMap = CALSTRUCTURETMP.sampleSynOutMap[chrTmp];
            auto& sampleSynOutMapTmp =  CALSTRUCTURETMP.sampleSynOutMapTmp[chrTmp];
            
            // 检查参考基因组中是否有对应的染色体
            map<string, uint32_t>::const_iterator findIter1 = refLenMap.find(chrTmp);
            if (findIter1 == refLenMap.end())
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The reference genome does not include '" << chrTmp << "'." << endl;
                exit(1);
            }

            uint32_t refLen = findIter1->second;  // 染色体长度

            // 如果到达末尾，跳出 loop
            if (refPosTmp >= refLen)
            {
                break;
            }
            
            // 如果没到达染色体末尾
            uint32_t preSynNum = 0;
            for (uint32_t refLoci = refPosTmp + 1; refLoci < refLen + 1; refLoci++)
            {
                uint32_t synNumTmp = 0;  // 共线性数量，临时

                for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
                {
                    string sampleName2 = SynCoorTmp.sampleNameVec[i];

                    // 如果是reference与大家没坐标，那就是都没坐标
                    if (sampleName == "reference")
                    {
                        continue;
                    }
                    
                    // 如果是其它
                    uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName
                    uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName2

                    if (lociA == 0 && lociB == 0) // both del
                    {
                        synNumTmp++;  // 记为 syn

                        // debug
                        if (debug)
                        {
                            sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                        }
                        
                        // debug
                        if (debug)
                        {
                            cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                        }
                    }
                    else if (lociA != 0 && lociB == 0)  // 一个有坐标，另一个没有，one del
                    {
                        continue;
                    }
                    else if (lociA == 0 && lociB != 0)  // 一个有坐标，另一个没有，one del
                    {
                        continue;
                    }
                    else  // 都有坐标，再判断坐标的A和B是不是共线性
                    {
                        if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A和B是共线性
                        {
                            synNumTmp++;  // 记为 syn

                            // debug
                            if (debug)
                            {
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                            }
                        }
                    }
                }

                // 如果和上一个的 synNum 不一样则添加
                if (synNumTmp != preSynNum)
                {
                    sampleSynOutMap[refLoci] += synNumTmp;
                    preSynNum = synNumTmp;
                }
            }
        }

        return CALSTRUCTURETMP;
    }



    /**
     * @brief 合并结果
     * 
     * @param chromosome                         染色体号
     * @param chrLen                             染色体长度
     * @param CALSTRUCTUREVec                    多线程输出结果
     * 
     * @return tuple<chromosome, synOutMapTmp>   map<refLoci, synNum>
    **/
    tuple<string, map<uint32_t, uint32_t> > merge(
        string chromosome, 
        uint32_t chrLen, 
        const vector<CALSTRUCTURE> & CALSTRUCTUREVec
    )
    {
        // 保存结果
        map<uint32_t, uint32_t> synOutMapTmp;  //  map<refLoci, synNum>
        
        for (const auto& CALSTRUCTURETmp : CALSTRUCTUREVec) // vector<CALSTRUCTURE>
        {
            // 染色体对应的迭代器
            map<string, map<uint32_t, uint32_t> >::const_iterator iter0 = CALSTRUCTURETmp.sampleSynOutMap.find(chromosome);
            map<string, map<uint32_t, string> >::const_iterator iter3 = CALSTRUCTURETmp.sampleSynOutMapTmp.find(chromosome);

            // 临时迭代器，记录当前和下一个位点的得分
            map<uint32_t, uint32_t>::const_iterator iter1;
            map<uint32_t, uint32_t>::const_iterator iter2;
            
            if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end())  // 如果含有该染色体
            {
                iter1 = iter0->second.begin();
                iter2 = iter0->second.begin();
                iter2++;
            }
            
            // 记录上一个节点的 synNum
            uint32_t preSynNumTmp = 0;

            // 循环染色体长度，每个位置单独添加
            for (uint32_t refLoci = 1; refLoci < chrLen + 1; refLoci++)
            {
                // 判读当前位置是否指向第一个迭代器
                if (refLoci == iter1->first && iter0 != CALSTRUCTURETmp.sampleSynOutMap.end())  // 如果是，赋值该迭代器的数量
                {
                    preSynNumTmp = iter1->second;
                }
                else if (refLoci == iter2->first && iter0 != CALSTRUCTURETmp.sampleSynOutMap.end() && iter2 != iter0->second.end())  // 走到下一个节点，更新迭代器
                {
                    iter1++;  // 更新迭代器
                    iter2++;  // 更新迭代器

                    preSynNumTmp = iter1->second;  // 更新坐标
                }

                // 位点频率叠加
                synOutMapTmp[refLoci] += preSynNumTmp;

                // debug
                if (debug)
                {
                    if (preSynNumTmp != 0)
                    {
                        cerr << refLoci << " " << CALSTRUCTURETmp.sampleName << " " << preSynNumTmp << " " << iter3->second.at(refLoci) << endl;
                    }
                }
            }
        }

        return make_tuple(chromosome, synOutMapTmp);
    }


    /**
     * @brief 计算
     * 
     * @param refLenMap                参考基因组长度
     * @param SynCoorTmp               共线性坐标
     * @param sampleSampleSyriMap      syri.out 输出路径字典
     * @param alignsMap                show-aligns 输出路径字典
     * @param allSynNum                所有组合数量
     * @param outputFileName           输出文件名
     * 
     * @return 0
    **/
    int calculate(
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap, 
        const uint32_t & allSynNum, 
        const string & outputFileName
    )
    {
        /* ********************************************** calculate syntenic diversity ********************************************** */
        // 进程池
        ThreadPool pool(threads);

        // 初始化线程池
        pool.init();

        // 保存多线程的结果
        /*
            map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >
            map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
        */
        vector<future<CALSTRUCTURE> > calOutStrVecTmp;
        
        // 样品名称
        for (const auto& iter1 : SynCoorTmp.sampleNameVec)  // vector<sampleName>
        {
            string sampleNameTmp = iter1;

            // 最后一个样本跳过
            if (sampleNameTmp == SynCoorTmp.sampleNameVec[SynCoorTmp.sampleNameVec.size() - 1])
            {
                continue;
            }

            cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate: " << sampleNameTmp << endl;

            // 多线程提交并保存结果
            calOutStrVecTmp.push_back(
                pool.submit(
                    calculate_run, 
                    sampleNameTmp, 
                    ref(refLenMap), 
                    ref(SynCoorTmp), 
                    ref(sampleSampleSyriMap), 
                    ref(alignsMap)
                )
            );

            sleep(0.1);  // 线程间隔
        }

        // 多线程结果保存
        vector<CALSTRUCTURE> CALSTRUCTUREVecTmp;
        for (size_t i = 0; i < calOutStrVecTmp.size(); i++)
        {
            CALSTRUCTUREVecTmp.push_back(move(calOutStrVecTmp[i].get()));
        }
        calOutStrVecTmp.clear();
        vector<future<CALSTRUCTURE> >().swap(calOutStrVecTmp);
        malloc_trim(0);	// 0 is for heap memory


        /* ********************************************** merge the result ********************************************** */
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Merge the result." << endl;
        vector<future<tuple<string, map<uint32_t, uint32_t> > > > mergeOutVec;
        for (const auto& iter1 : refLenMap)  // map<chromosome, length>
        {
            string chromosome = iter1.first;
            uint32_t chrLen = iter1.second;

            // 多线程提交并保存结果
            mergeOutVec.push_back(
                pool.submit(
                    merge, 
                    chromosome, 
                    chrLen, 
                    ref(CALSTRUCTUREVecTmp)
                )
            );

            sleep(0.1);  // 线程间隔
        }
        

        /* ********************************************** save the result ********************************************** */
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Save the result." << endl;
        SAVE::SAVE SAVEClass(outputFileName);

        stringstream outStream; // 使用 stringstream 代替字符串拼接
        constexpr uint64_t CACHE_SIZE = 1024 * 1024 * 10; // 缓存大小为 10mb

        for (size_t i = 0; i < mergeOutVec.size(); i++)  // vector<future<tuple<string, map<uint32_t, uint32_t> > > > 
        {
            string chomosome;
            map<uint32_t, uint32_t> calOutMap;  // map<uint32_t, uint32_t>

            tie(chomosome, calOutMap) = move(mergeOutVec[i].get());  // 获取多线程结果

            for (const auto& iter1 : calOutMap)
            {
                outStream << chomosome << "\t" << iter1.first << "\t" << allSynNum << "\t" 
                        << iter1.second << "\t" << iter1.second/(float)allSynNum << "\n";

                if (outStream.tellp() >= CACHE_SIZE) {  // 缓存大小为 10mb
                    string outTxt = outStream.str();
                    SAVEClass.save(outTxt);
                    // 清空 stringstream
                    outStream.str(string());
                    outStream.clear();
                }
            }

            if (outStream.tellp() >= 0) {  // 最后写一次
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // 清空 stringstream
                outStream.str(string());
                outStream.clear();
            }
        }

        // 关闭线程池
        pool.shutdown();

        return 0;
    }
}

#endif