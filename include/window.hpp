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
    // �ִ���ͳ��
    class WINDOW
    {
    private:
        // Ⱦɫ�峤��
        map<string, uint32_t> refLenMap_;  // map<chr, length>

        // cal����ĵ÷�
        string calFileName_;

        // ����ļ�
        string outputFileName_;

        // ���ںͲ���
        uint32_t windowSize_;
        uint32_t stepSize_;

        // cal������
        map<string, vector<double> > chrScoreVecMap_;  // map<chr, vector<score> >
        
    public:
        // ���ڻ����Լ�ƽ���÷�
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
         * @brief ����Ⱦɫ�峤�ȴ���window
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
                // ������
                auto& WinInfoTupVec = chrWinInfoTupVecMap_[chromosome];

                // ����
                double stepsNumber = double(length)/stepSize_;

                // ����Ǹ�������������ȡ��
                if (stepsNumber - int(stepsNumber) != 0)
                {
                    stepsNumber = int(stepsNumber+1);
                }

                // ���������ֵ�
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
                        // ���ܳ���Ⱦɫ�峤��
                        WinInfoTupVec.push_back(make_tuple(stepStart, min(stepEnd, length), 0.0));

                        // ���һ��������λ��С��Ⱦɫ�峤��
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
         * @brief ��cal�ļ�����������
         * 
         * @return 0
        **/
        int cal_open()
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Build index for " << calFileName_ << endl;

            string chromosome;
            uint32_t length = 0;

            // ������
            vector<double>* ScoreVec;

            // open file
            GzChunkReader GzChunkReaderClass(calFileName_);
            // read line
            string line;
            while (GzChunkReaderClass.read_line(line))
            {
                // ��������
                if (line.empty() || line.find("#") != string::npos)
                {
                    continue;
                }

                // ���
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // ��������Ƿ�Խ��
                if (lineVec.size() != 5)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of columns in the file is not 5." << endl;
                    exit(1);
                }

                string chromosomeTmp = lineVec[0];  // ��ʱȾɫ���
                uint32_t position = stoul(lineVec[1]);  // λ��
                string scoreS = lineVec[4];

                // ����Ƿ�������
                if (!isdigit(scoreS[0]))
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Identify non-numeric data -> " << scoreS << endl;
                    exit(1);
                }

                if (chromosomeTmp != chromosome)
                {
                    chromosome = chromosomeTmp;

                    // ���Ⱦɫ���Ƿ��ڲο���������
                    if (refLenMap_.find(chromosome) == refLenMap_.end())
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The reference genome does not contain " << chromosome << "." << endl;
                        exit(1);
                    }

                    length = refLenMap_[chromosome];

                    // ��ʼ��
                    chrScoreVecMap_[chromosomeTmp] = vector<double>(length, 0.0);
                    ScoreVec = &chrScoreVecMap_[chromosomeTmp];
                }
                
                double score = stod(scoreS);

                // ��ֵ
                (*ScoreVec)[position] = score;
            }

            return 0;
        }


        /**
         * @brief �������������vector��ƽ��ֵ
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
         * @brief ����ƽ��ֵ
         * 
         * @return 0
        **/
        int win_count()
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate the mean score for each window ..." << endl;

            for (auto& [chromosome, WinInfoTupVec] : chrWinInfoTupVecMap_)  // map<chr, vector<tuple<start, end, score> > >
            {
                // Ⱦɫ���Ӧ��score vector   map<chr, vector<score> >
                const auto& ScoreVec = chrScoreVecMap_[chromosome];

                for (auto& WinInfoTup : WinInfoTupVec)  // vector<tuple<start, end, score> >
                {
                    get<2>(WinInfoTup) = average(ScoreVec, get<0>(WinInfoTup) - 1, get<1>(WinInfoTup) - 1);
                }
            }
            
            return 0;
        }


        /**
         * @brief ������
         * 
         * @return 0
        **/
        int save_result()
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName_ << "'" << endl;

            SAVE SAVEClass(outputFileName_);

            stringstream outStream; // ʹ�� stringstream �����ַ���ƴ��
            static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // �����СΪ 10mb
            outStream.str().reserve(CACHE_SIZE);
            outStream << "#CHROM\tSTART\tEND\tSyntenic_Diversity\n";

            for (auto& [chromosome, WinInfoTupVec] : chrWinInfoTupVecMap_)  // map<chr, vector<tuple<start, end, score> > >
            {
                for (auto& WinInfoTup : WinInfoTupVec)  // vector<tuple<start, end, score> >
                {
                    outStream << chromosome << "\t" << get<0>(WinInfoTup) << "\t" << get<1>(WinInfoTup) << "\t" 
                            << get<2>(WinInfoTup) << "\n";
                    if (outStream.tellp() >= CACHE_SIZE)  // �����СΪ 10mb
                    {
                        string outTxt = outStream.str();
                        SAVEClass.save(outTxt);
                        // ��� stringstream
                        outStream.str(string());
                        outStream.clear();
                    }
                }

                if (outStream.tellp() > 0)  // ���дһ��
                {
                    string outTxt = outStream.str();
                    SAVEClass.save(outTxt);
                    // ��� stringstream
                    outStream.str(string());
                    outStream.clear();
                }
            }
            return 0;
        }
    };
}

#endif