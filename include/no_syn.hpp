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
#include <getopt.h>

#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "save.hpp"
#include "GzChunkReader.hpp"
#include "ThreadPool.hpp"
#include "cal.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

void help_no_syn(char* argv[]);
int main_no_syn(int argc, char* argv[]);

// 调试代码
extern bool debugNoSyn;

namespace NOSYN
{
    // kseq.h 打开文件
    KSEQ_INIT(gzFile, gzread)

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


    // 将 CALNAME::SYNCOOR 作为基类
    class NOSYNCOOR:public CALNAME::SYNCOOR
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
}

#endif