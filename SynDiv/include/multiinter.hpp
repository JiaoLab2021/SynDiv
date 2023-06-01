#ifndef MULTIINTER_hpp
#define MULTIINTER_hpp
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include "save.hpp"
#include "GzChunkReader.hpp"

using namespace std;

namespace MULTIINTER
{
    bool debug;
    // bool debug = false;
    // bool debug = true;

    pair<string, map<string, vector<pair<int, int> > > > build_syn_index(
        const string & inputFileName
    );

    tuple<int, int, vector<string> > loc_find(
        const vector<pair<int, int > > synVec, 
        const map<int, string> & idxLineMap
    );

    int save_result(
        const map<string, map<int, vector<tuple<int, vector<string> > > > > & outChrStartEndLineVecMap, 
        map<int, string> idxLineMap, 
        const string & outputName
    );



    /**
     * @brief build syn index for syri.out
     * 
     * @param inputFileName   the output of syri
     * 
     * @return pair<lineName, chrSynVecMap>        pair<lineName, map<chromosome, vector<pair<refStart, refEnd>>>>
    **/
    pair<string, map<string, vector<pair<int, int> > > > build_syn_index(
        const string & inputFileName
    )
    {
        // 品种的名字
        string lineName = inputFileName;

        vector<string> inputFileNameVecTmp = split(inputFileName, "/");  // 路径拆分
        lineName = inputFileNameVecTmp[inputFileNameVecTmp.size() - 1];  // 最后一个

        string prefix = ".syri.out";
        auto iter1 = lineName.find(prefix);
        if(iter1 != string::npos)
        {
            lineName.replace(lineName.begin() + iter1, lineName.begin() + iter1 + prefix.size(), "");
        }
        
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Building index: " << inputFileName << ".\n";  // print log

        map<string, vector<pair<int, int> > > chrSynVecMap;  // map<chromosome, vector<pair<refStart, refEnd>>>

        // open file
        GzChunkReader GzChunkReaderClass(inputFileName);

        // 上一个共线性的坐标
        string preChromosome;
        int PreRefStart = 0;
        int PreRefEnd = 0;

        // read line
        string line;

        while (GzChunkReaderClass.read_line(line))
        {
            // 跳过空行
            if (line.empty() || line.find("SYNAL") == string::npos)
            {
                continue;
            }

            // 拆分
            std::istringstream iss(line);
            vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

            string chromosome = lineVec[0];
            int refStart = stoi(lineVec[1]);
            int refEnd = stoi(lineVec[2]);

            if(preChromosome != chromosome && preChromosome.size() > 0)  // 如果换染色体了/第一条染色体略过。清零
            {
                chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));
                PreRefStart = refStart;
                PreRefEnd = refEnd;
            }
            else
            {
                if(PreRefStart == 0)  // 第一条染色体第一个syn
                {
                    PreRefStart = refStart;
                    PreRefEnd = refEnd;
                }
                else if(refStart <= PreRefEnd)  // 上一个包含该区间
                {
                    if(refEnd > PreRefEnd)  // 部分包含
                    {
                        PreRefStart = PreRefStart;
                        PreRefEnd = refEnd;
                    }
                    else  // 完全包含
                    {
                        PreRefStart = PreRefStart;
                        PreRefEnd = PreRefEnd;
                    }
                }
                else  // 不重叠
                {
                    chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));
                    PreRefStart = refStart;
                    PreRefEnd = refEnd;
                }
            }

            preChromosome = chromosome;
        }

        // 最后一个坐标添加
        chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));

        return make_pair(lineName, chrSynVecMap);
    }


    /**
        * @brief 找合集
        * 
        * @param chrLineSynVecMap          map<string, map<string, vector<pair<int, int> > > >, map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
        * 
        * @return pair(outChrStartEndLineVecMap, idxLineMap)   pair(map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>, 存储line的索引和名字)
    **/
    pair<map<string, map<int, vector<tuple<int, vector<string> > > > >, map<int, string> > syn_multiinter_find(
        const map<string, map<string, vector<pair<int, int> > > > & chrLineSynVecMap
    )
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Searching." << endl;
        // 保存结果
        map<string, map<int, vector<tuple<int, vector<string> > > > > outChrStartEndLineVecMap;  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>

        map<int, string> idxLineMap;  // 存储line的索引和名字

        // 遍历染色体
        for(auto it1 : chrLineSynVecMap)  // map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
        {
            // it1.first -> chromosome
            string chromosome = it1.first;

            // it1.second -> map<lineName, vector<pair<refStart, refEnd>>>

            int lineIdx = 0;  // lineName的索引
            
            // 将染色体对应的所有共线性列表和索引信息提取出来，并赋值第一个坐标为syn列表的第一个
            map<string, tuple<vector<pair<int, int> >, int, pair<int, int> > > LineSynVecIdxSynMap;  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
            for(auto it2 : it1.second)  // map<lineName, vector<pair<refStart, refEnd>>>
            {
                // it2.first -> lineName
                // it2.second -> vector<pair<refStart, refEnd>>

                LineSynVecIdxSynMap[it2.first] = make_tuple(it2.second, 0, it2.second[0]);

                idxLineMap[lineIdx] = it2.first;  // 存储line信息
                lineIdx++;  // 索引加1
            }


            // 寻找，先构造一个临时的共线性列表，存储每个样本第一个共线性区间
            vector<pair<int, int> > synVec;  // 存储共线性区间，找合集
            
            int lineNum = idxLineMap.size();  // sample数量，用于判断是否终止循环
            int endNum = 0;  // 用于判断是否终止while循环

            for(auto it1 : LineSynVecIdxSynMap)  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
            {
                // it1.first -> lineName
                string lineName = it1.first;

                // it1.second -> tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> >
                synVec.push_back(get<2>(it1.second));

                // 如果调试代码的话，打印log
                if (debug)
                {
                    cerr << chromosome << " " << lineName << ":" << get<2>(it1.second).first << "-" << get<2>(it1.second).second << endl;
                }
                
                if(get<1>(it1.second) == get<0>(it1.second).size())  // 判断是否循环完毕
                {
                    endNum++;
                }
            }
            
            while(endNum < lineNum)  // 当 endNum == lineNum 时结束循环
            {
                endNum = 0;  // 重置计数

                // 找交集
                int outRefStart = 0;
                int outRefEnd = 0;
                vector<string> outLineVec;
                tie(outRefStart, outRefEnd, outLineVec) = loc_find(
                    synVec, 
                    idxLineMap
                );

                // 如果调试代码的话，打印log
                if (debug)
                {
                    cerr << chromosome << " loc:" << outRefStart << "-" << outRefEnd << endl << chromosome << " outLine:";

                    for (auto outLine : outLineVec)
                    {
                        cerr << " " << outLine;
                    }
                    cerr << endl << endl;
                }

                // 存储交集结果
                outChrStartEndLineVecMap[chromosome][outRefStart].push_back(make_tuple(outRefEnd, outLineVec));

                // 清空坐标
                synVec.clear();
                vector<pair<int, int > >().swap(synVec);

                // 更新共线性索引和坐标
                for(auto it1 : LineSynVecIdxSynMap)  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
                {
                    // it1.first -> lineName
                    string lineName = it1.first;

                    // 品种对应的共线性坐标，用于判断是否更新坐标
                    vector<pair<int, int> > lineSynVec = get<0>(it1.second);
                    int lineSynIdx = get<1>(it1.second);
                    pair<int, int> lineSycLocTmp = get<2>(it1.second);

                    // 更新坐标
                    if(lineSycLocTmp.first >= outRefEnd)  // 如果没用到该共线性区间，直接下一个品种
                    {
                        synVec.push_back(lineSycLocTmp);  // 更新坐标，还用上一个的坐标
                    }
                    else if(outRefEnd >= lineSycLocTmp.second)  // 共线性坐标被用完了
                    {
                        lineSynIdx++;
                        if(lineSynIdx >= lineSynVec.size())  // 判断是否循环完毕
                        {
                            lineSycLocTmp = make_pair(0, 0);  // 全归零
                            endNum++;
                        }
                        else
                        {
                            lineSycLocTmp = lineSynVec[lineSynIdx];  // 共线性指向下一个坐标
                        }

                        synVec.push_back(lineSycLocTmp);
                        LineSynVecIdxSynMap[lineName] = make_tuple(lineSynVec, lineSynIdx, lineSycLocTmp);  // 更新哈希表的坐标
                    }
                    else  // 共线性坐标被用了一半
                    {
                        lineSycLocTmp.first = outRefEnd;

                        synVec.push_back(lineSycLocTmp);  // 更新坐标
                        LineSynVecIdxSynMap[lineName] = make_tuple(lineSynVec, lineSynIdx, lineSycLocTmp);  // 更新哈希表的坐标
                    }

                    // 如果调试代码的话，打印log
                    if (debug)
                    {
                        cerr << chromosome << " " << lineName << ":" << lineSycLocTmp.first << "-" << lineSycLocTmp.second << " endNum:" << endNum << endl;
                    }
                }
            }
        }

        return make_pair(outChrStartEndLineVecMap, idxLineMap);
    }


    /**
        * @brief 找集合_find
        * 
        * @param synVec      vector<pair<int, int > >, vector<pair<refStart, refEnd> >
        * 
        * @return tuple<int, int, vector<string> >  make_tuple(outStart, outEnd, lineName)
    **/
    tuple<int, int, vector<string> > loc_find(
        const vector<pair<int, int > > synVec, 
        const map<int, string> & idxLineMap
    )
    {
        // 共线性区间
        int outStart = INT_MAX;
        int outEnd = INT_MAX;

        vector<string> lineNameVec;  // 共线性区间对应的品种索引
        

        for(auto it1 : synVec)  // 找共线性最小的坐标，肯定在起始坐标中
        {
            if(it1.first == 0 && it1.second == 0) continue;  // 如果都为0代表这个line循环完毕

            outStart = min(outStart, it1.first);  // 所有共线性中的最小值，肯定在起始坐标中
        }

        for(auto it1 : synVec)  // 找共线性第二小的坐标
        {
            if(it1.first == 0 && it1.second == 0) continue;  // 如果都为0代表这个line循环完毕

            if (it1.first != outStart)  // 不能是最小值
            {
                outEnd = min(outEnd, it1.first);
            }
            outEnd = min(outEnd, it1.second);
        }

        if (outEnd == INT_MAX)  // 防止出现终止和起始一样的情况
        {
            outEnd = outStart;
        }
        

        int lineIdx = 0;
        for(auto it1 : synVec)  // 找最小区间对应的品种信息
        {
            if(it1.first == 0 && it1.second == 0)  // 如果都为0代表这个line循环完毕
            {
                lineIdx++;
                continue;
            }

            if(it1.first >= outStart && it1.first < outEnd)
            {
                string lineName;
                auto iter = idxLineMap.find(lineIdx);
                if(iter != idxLineMap.end())
                {
                    lineName = iter->second;
                }
                else
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "keyError: " << lineIdx << " not in idxLineMap." << endl;
                    exit(1);
                }
                lineNameVec.push_back(lineName);
            }
            lineIdx++;
        }
        
        return make_tuple(outStart, outEnd, lineNameVec);
    }


    /**
        * @brief 存储结果
        * 
        * @param outChrStartEndLineVecMap      const map<chr, map<refStart, vector<tuple<refEnd, vector<lineName> > > > >
        * @param idxLineMap                  map<index, lineName>
        * @param outputName                  输出文件名
        * 
        * @return tuple<int, int, vector<string> >  make_tuple(outStart, outEnd, lineName)
    **/
    int save_result(
        const map<string, map<int, vector<tuple<int, vector<string> > > > > & outChrStartEndLineVecMap, 
        map<int, string> idxLineMap, 
        const string & outputName
    )
    {
        SAVE::SAVE SAVEClass(outputName);

        stringstream outStream; // 使用 stringstream 代替字符串拼接
        static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // 缓存大小为 10mb
        outStream.str().reserve(CACHE_SIZE);

        // 将lineName添加到列表中
        vector<string> lineNameVec;
        for(auto it : idxLineMap)
        {
            lineNameVec.push_back(it.second); 
        }

        // 保存结果
        for(auto it1 : outChrStartEndLineVecMap)  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName> > > > >
        {
            for(auto it2 : it1.second)  // map<refStart, vector<tuple<refEnd, vector<lineName> > > >
            {
                for(auto it3 : it2.second)  // vector<tuple<refEnd, vector<lineName> > >
                {
                    outStream << it1.first << "\t" << to_string(it2.first) << "\t" + to_string(get<0>(it3)) << "\t" << to_string(get<1>(it3).size()) << "\t";

                    // 共线性含有的品种名字索引
                    int idxTmp=0;

                    // 输出(第五列及之后的内容)
                    vector<string> lineNameVecTmp;
                    vector<int> boolVecTmp;
                    for (size_t i = 0; i < lineNameVec.size(); i++)  // vector<lineName>
                    {
                        if(lineNameVec[i] == get<1>(it3)[idxTmp])  // 含有对应的品种
                        {
                            lineNameVecTmp.push_back(lineNameVec[i]);
                            boolVecTmp.push_back(1);
                            idxTmp++;
                        }
                        else
                        {
                            boolVecTmp.push_back(0);
                        }

                        // 如果到vector的结尾，跳出for循环
                        if (idxTmp == get<1>(it3).size())
                        {
                            break;
                        }  
                    }

                    // 如果结果为空，跳过该位点
                    if (lineNameVecTmp.size() == 0)
                    {
                        continue;
                    }
                    
                    outStream << join(lineNameVecTmp, ",") << "\t" + join(boolVecTmp, "\t") << "\n";

                    if (outStream.tellp() >= CACHE_SIZE)  // 缓存大小为 10mb
                    {
                        string outTxt = outStream.str();
                        SAVEClass.save(outTxt);
                        // 清空 stringstream
                        outStream.str(string());
                        outStream.clear();
                    }
                }
            }
        }

        if (outStream.tellp() >= 0)  // 最后写一次
        {
            string outTxt = outStream.str();
            SAVEClass.save(outTxt);
            // 清空 stringstream
            outStream.str(string());
            outStream.clear();
        }
        
        return 0;
    }
}

#endif