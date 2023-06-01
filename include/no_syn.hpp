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

// ���Դ���
extern bool debugNoSyn;

namespace NOSYN
{
    // kseq.h ���ļ�
    KSEQ_INIT(gzFile, gzread)

    // ���������������Ⱦɫ�峤��
    class LENGTH
    {
    private:
        // Ⱦɫ�峤���ļ�
        vector<string> lengthsVec_;
        vector<string> lengthsTitles_;

        // ������Ʒ�ĳ�����Ϣ
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
         * @brief ���ļ�
         * 
         * @return void
        **/
        void index_lengths()
        {
            for (size_t i = 0; i < lengthsVec_.size(); i++)
            {
                // Ⱦɫ�峤���ļ�
                string lengthFileName = lengthsVec_[i];

                // title
                string sampleName = "";


                // �������title�����б���ģ���������У������ļ�����ȡ
                if (lengthsTitles_.size() > 0)
                {
                    sampleName = lengthsTitles_[i];
                }
                else  // ������ȡ
                {
                    vector<string> lengthFileNameVec = split(lengthFileName, "/");  // ·�����
                    string baseName = lengthFileNameVec[lengthFileNameVec.size() - 1];  // �ļ���
                    vector<string> baseNameVec = split(baseName, "."); // �� '.' ���
                    sampleName = baseNameVec[0];  // ��Ʒ��Ϊ�б��е�һ��Ԫ��
                }

                // ��ʼ���ֵ�
                sampleChrLenMap[sampleName];

                // open file
                GzChunkReader GzChunkReaderClass(lengthFileName);
                // read line
                string line;
                while (GzChunkReaderClass.read_line(line))
                {
                    // ��������
                    if (line.empty())
                    {
                        continue;
                    }

                    // ���
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    // ��������Ƿ���Ϲ涨
                    if (lineVec.size() != 2)
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Error: the number of columns is not 2 -> " << line << endl;
                        exit(1);
                    }

                    // Ⱦɫ��źͳ�����Ϣ
                    string chromosome = lineVec[0];
                    uint32_t chrLen = stol(lineVec[1]);
                    
                    // �洢
                    sampleChrLenMap[sampleName][chromosome] = chrLen;
                }
            }
        }

        /**
         * @brief ��ȡȾɫ�峤��
         * 
         * @param sample ��Ʒ��
         * @param chr    Ⱦɫ���
         * 
         * @return uint32_t
        **/
        uint32_t get_length(
            const string & sample, 
            const string & chr
        )
        {
            // ����Ⱦɫ�峤��
            // ����Ʒ
            auto fintIter1 = sampleChrLenMap.find(sample);
            // ���û�У�����
            if (fintIter1 == sampleChrLenMap.end())
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: '" << sample << "' not present in chromosome length file." << endl;
                exit(1);
            }
            // ��Ⱦɫ��
            auto findIter2 = fintIter1->second.find(chr);
            if (findIter2 == fintIter1->second.end())
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: the length of '" << chr << "' was not found -> " << sample << endl;
                exit(1);
            }

            // Ⱦɫ�峤��
            return findIter2->second;
        }
    };


    // �� CALNAME::SYNCOOR ��Ϊ����
    class NOSYNCOOR:public CALNAME::SYNCOOR
    {
    private:

    public:
        // ��¼�ǹ����Ե�����
        map<string, map<uint32_t, map<string, vector<tuple<uint32_t, uint32_t> > > > > chrStartSampleLociVecMap; // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >

        // ����ַ���
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