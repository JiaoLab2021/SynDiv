#ifndef WINDOW_HPP
#define WINDOW_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include <numeric>
#include <tuple>
#include <string.h>
#include <getopt.h>

#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "save.hpp"
#include "GzChunkReader.hpp"
#include "cal.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

void help_window(char* argv[]);
int main_window(int argc, char* argv[]);

namespace Window
{
    // 分窗口统计
    class WINDOW
    {
    private:
        // 染色体长度
        map<string, uint32_t> refLenMap_;  // map<chr, length>

        // cal计算的得分
        string calFileName_;

        // 输出文件
        string outputFileName_;

        // 窗口和步长
        uint32_t windowSize_;
        uint32_t stepSize_;

        // cal的索引
        map<string, vector<double> > chrScoreVecMap_;  // map<chr, vector<score> >
        
    public:
        // 窗口划分以及平均得分
        map<string, vector<tuple<uint32_t, uint32_t, double> > > chrWinInfoTupVecMap_;  // map<chr, vector<tuple<start, end, score> > >

        WINDOW(
            const map<string, uint32_t>& refLenMap, 
            const string& calFileName, 
            string outputFileName = "",
            uint32_t windowSize = 5000,
            uint32_t stepSize = 1000
        )
        {
            refLenMap_ = refLenMap;
            calFileName_ = calFileName;
            outputFileName_ = outputFileName;
            windowSize_ = windowSize;
            stepSize_ = stepSize;
        }
        ~WINDOW() {}

        /**
         * @brief 根据染色体长度创建window
         * 
         * @return int
        **/
        int make_window()
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate the number of windows based on the input value ..." << endl;
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Window: " << windowSize_ << endl;
            cerr << "[" << __func__ << "::" << getTime() << "] " << "step: " << stepSize_ << endl;

            for (const auto& [chromosome, length] : refLenMap_)
            {
                // 变量绑定
                auto& WinInfoTupVec = chrWinInfoTupVecMap_[chromosome];

                // 步数
                double stepsNumber = double(length)/stepSize_;

                // 如果是浮点数，则向上取整
                if (stepsNumber - int(stepsNumber) != 0)
                {
                    stepsNumber = int(stepsNumber+1);
                }

                // 制作步长字典
                for (int i = 0; i < stepsNumber; i++)
                {
                    uint32_t stepStart = i * stepSize_ + 1;
                    uint32_t stepEnd = stepStart + windowSize_ - 1;
                    if (i < stepsNumber - 1)
                    {
                        WinInfoTupVec.push_back(make_tuple(stepStart, min(stepEnd, length), 0.0));
                    }
                    else
                    {
                        // 不能超过染色体长度
                        WinInfoTupVec.push_back(make_tuple(stepStart, min(stepEnd, length), 0.0));

                        // 最后一个窗口且位置小于染色体长度
                        if (stepEnd < length)
                        {
                            WinInfoTupVec.push_back(make_tuple(stepEnd + 1, length, 0.0));
                        }
                    }
                }
            }
            
            return 0;
        }


        /**
         * @brief 打开cal文件并构建索引
         * 
         * @return 0
        **/
        int cal_open()
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Build index for " << calFileName_ << endl;

            string chromosome;
            uint32_t length = 0;

            // 变量绑定
            vector<double>* ScoreVec;

            // open file
            GzChunkReader GzChunkReaderClass(calFileName_);
            // read line
            string line;
            while (GzChunkReaderClass.read_line(line))
            {
                // 跳过空行
                if (line.empty() || line.find("#") != string::npos)
                {
                    continue;
                }

                // 拆分
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // 检查数组是否越界
                if (lineVec.size() != 5)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of columns in the file is not 5." << endl;
                    exit(1);
                }

                string chromosomeTmp = lineVec[0];  // 临时染色体号
                uint32_t position = stoul(lineVec[1]);  // 位置
                string scoreS = lineVec[4];

                // 检查是否是数字
                if (!isdigit(scoreS[0]))
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Identify non-numeric data -> " << scoreS << endl;
                    exit(1);
                }

                if (chromosomeTmp != chromosome)
                {
                    chromosome = chromosomeTmp;

                    // 检查染色体是否在参考基因组中
                    if (refLenMap_.find(chromosome) == refLenMap_.end())
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The reference genome does not contain " << chromosome << "." << endl;
                        exit(1);
                    }

                    length = refLenMap_[chromosome];

                    // 初始化
                    chrScoreVecMap_[chromosomeTmp] = vector<double>(length, 0.0);
                    ScoreVec = &chrScoreVecMap_[chromosomeTmp];
                }
                
                double score = stod(scoreS);

                // 赋值
                (*ScoreVec)[position] = score;
            }

            return 0;
        }


        /**
         * @brief 计算给定索引下vector的平均值
         * 
         * @return 0
        **/
        double average(const vector<double>& v, uint32_t start_index, uint32_t end_index)
        {
            auto first = v.begin() + start_index;
            auto last = v.begin() + end_index + 1;
            double sum = accumulate(first, last, 0.0);
            return sum / (end_index - start_index + 1);
        }


        /**
         * @brief 计算平均值
         * 
         * @return 0
        **/
        int win_count()
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate the mean score for each window ..." << endl;

            for (auto& [chromosome, WinInfoTupVec] : chrWinInfoTupVecMap_)  // map<chr, vector<tuple<start, end, score> > >
            {
                // 染色体对应的score vector   map<chr, vector<score> >
                const auto& ScoreVec = chrScoreVecMap_[chromosome];

                for (auto& WinInfoTup : WinInfoTupVec)  // vector<tuple<start, end, score> >
                {
                    get<2>(WinInfoTup) = average(ScoreVec, get<0>(WinInfoTup) - 1, get<1>(WinInfoTup) - 1);
                }
            }
            
            return 0;
        }


        /**
         * @brief 保存结果
         * 
         * @return 0
        **/
        int save_result()
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName_ << "'" << endl;

            SAVE SAVEClass(outputFileName_);

            stringstream outStream; // 使用 stringstream 代替字符串拼接
            static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // 缓存大小为 10mb
            outStream.str().reserve(CACHE_SIZE);
            outStream << "#CHROM\tSTART\tEND\tSyntenic_Diversity\n";

            for (auto& [chromosome, WinInfoTupVec] : chrWinInfoTupVecMap_)  // map<chr, vector<tuple<start, end, score> > >
            {
                for (auto& WinInfoTup : WinInfoTupVec)  // vector<tuple<start, end, score> >
                {
                    outStream << chromosome << "\t" << get<0>(WinInfoTup) << "\t" << get<1>(WinInfoTup) << "\t" 
                            << get<2>(WinInfoTup) << "\n";
                    if (outStream.tellp() >= CACHE_SIZE)  // 缓存大小为 10mb
                    {
                        string outTxt = outStream.str();
                        SAVEClass.save(outTxt);
                        // 清空 stringstream
                        outStream.str(string());
                        outStream.clear();
                    }
                }

                if (outStream.tellp() > 0)  // 最后写一次
                {
                    string outTxt = outStream.str();
                    SAVEClass.save(outTxt);
                    // 清空 stringstream
                    outStream.str(string());
                    outStream.clear();
                }
            }
            return 0;
        }
    };
}

#endif