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
    // ���Դ���
    bool debug = false;
    // bool debug = true;
    
    // kseq.h ���ļ�
    KSEQ_INIT(gzFile, gzread)



    /**
     * @brief ����ȫ�ֲ���
     * 
     * @param debugTmp     �Ƿ���Դ���
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


    // �� CAL::SYNCOOR ��Ϊ����
    class NOSYNCOOR:public CAL::SYNCOOR
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


    /**
     * @brief Ѱ�ҷǹ����Ե�����
     * 
     * @param lengthClass  ��ƷȾɫ�峤��class
     * 
     * @return void
    **/
    void NOSYNCOOR::find_no_syn_coor(
        LENGTH & lengthClass
    )
    {   
        // ��������������
        for (auto iter1 = coorChrLociSampleLociMap.begin(); iter1 != coorChrLociSampleLociMap.end(); iter1++)  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
        {
            string chrTmp = iter1->first;  // Ⱦɫ���

            // �ǹ����Ե��������ж��Ƿ���Ҫ���Ĳο�����������꣬��ֹ�ظ��������
            uint32_t preNoSynNum = 0;

            uint32_t refStart;  // �ο����������ʼ
            uint32_t refEnd;  // �ο����������ֹ


            // �����ж���Ʒ��С�����еķǹ����������Ƿ��ҵ�
            unordered_map<string, bool> sampleNoSynFindBoolMap;  // map<sample, bool>   true->�ҵ���   false->û�ҵ�
            // ��ʼ��ÿ����Ʒ�ķǹ����������Ƿ��ҵ�
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
                uint32_t refStartTmp;  // �ο����������ʼ
                uint32_t refEndTmp;  // �ο����������ֹ

                // �ڵ����ʼ����ֹλ��
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
                
                // �ǹ����Ե��������ж��Ƿ���Ҫ���Ĳο�����������꣬��ֹ�ظ��������
                uint32_t noSynNum = 0;

                // �����һ���ڵ�ȫ�ǹ����Եģ����²ο������������
                if (preNoSynNum  == 0)
                {
                    refStart = refStartTmp;
                    refEnd = refEndTmp;
                }

                for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++)  // map<sample, tuple(start, end)>
                {
                    string sampleTmp = iter3->first;  // ��Ʒ��
                    uint32_t StartTmp = get<0>(iter3->second);  // ��ʼ
                    uint32_t EndTmp = get<1>(iter3->second);  // ��ֹ

                    // �ο��������ֱ������ 
                    if (sampleTmp == "reference")
                    {
                        continue;
                    }

                    // ��¼Ϊ�ǹ����Ե�
                    if (StartTmp == 0 || EndTmp == 0)
                    {
                        // �ǹ�����������1
                        noSynNum++;

                        // �ҷǹ����Ե�λ��
                        // �ȳ�ʼ���ֵ�
                        if (chrStartSampleLociVecMap.find(chrTmp) == chrStartSampleLociVecMap.end())  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
                        {
                            chrStartSampleLociVecMap[chrTmp];
                        }
                        if (chrStartSampleLociVecMap[chrTmp].find(refStart) == chrStartSampleLociVecMap[chrTmp].end())  // map<refStart, map<sample, vector<tuple<start, end> > > >
                        {
                            chrStartSampleLociVecMap[chrTmp][refStart];
                        }

                        // �����sample���걻�ҵ��ˣ�����
                        if (sampleNoSynFindBoolMap[sampleTmp])  // map<sample, findBool>
                        {
                            continue;
                        }

                        // ��ʱ����
                        uint32_t noSynStartTmp = 1;
                        uint32_t noSynEndTmp = 1;

                        auto iter2Tmp = iter2;  // ��ʱ�ڵ�
                        
                        if (iter2 == iter1->second.begin())  // ����ǵ�һ���ڵ�
                        {
                            noSynStartTmp = 1;
                        }
                        else
                        {
                            iter2Tmp = iter2;
                            iter2Tmp--;  // �����ƶ�һλ

                            // �����һ���ڵ㻹��0���Ҳ�Ϊ��һ���ڵ㣬����ǰ��
                            while (iter2Tmp != iter1->second.begin() && get<1>(iter2Tmp->second[sampleTmp]) == 0)
                            {
                                iter2Tmp--;  // �����ƶ�һλ
                            }

                            // ��ֵǰһ������
                            if (get<1>(iter2Tmp->second[sampleTmp]) != 0)
                            {
                                noSynStartTmp = get<1>(iter2Tmp->second[sampleTmp]) + 1;
                            }
                        }

                        
                        iter2Tmp = iter2;  // ��ʼ��������
                        iter2Tmp++;  // �����ƶ�һλ

                        // �����һ���ڵ㻹��0���Ҳ�Ϊ���һ���ڵ㣬��������
                        while (iter2Tmp != iter1->second.end() && get<0>(iter2Tmp->second[sampleTmp]) == 0)
                        {
                            iter2Tmp++;  // �����ƶ�һλ
                        }

                        // ��������һ���ڵ㣬��ΪȾɫ�峤��
                        if (iter2Tmp == iter1->second.end())
                        {
                            // Ⱦɫ�峤��
                            noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);
                        }
                        else
                        {
                            // ��һ���ڵ����ʼ���� -1
                            noSynEndTmp = get<0>(iter2Tmp->second[sampleTmp]) - 1;
                        }

                        // ��ʼ����С����ֹ
                        if (noSynStartTmp < noSynEndTmp)
                        {
                            // ��ʼ���ֵ�
                            if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                            {
                                chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                            }
                            
                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(noSynStartTmp, noSynEndTmp));  // ��¼�ǹ���������

                            sampleNoSynFindBoolMap[sampleTmp] = true;  // ��¼����Ʒ���ҵ���
                        }
                    }
                    else
                    {
                        /* �ж��Ƿ�Ϊ���һ���ڵ� */
                        auto iter2Tmp = iter2;  // ��ʱ�ڵ�
                        iter2Tmp++;  // �����ƶ�һλ

                        // ����ǵ�һ���ڵ㣬�ж��Ƿ�� 1 ��ʼ�ıȶԡ�����Ͳο��������ǹ����ԵĻ���������Ⱦɫ��ǰ�벿�֣�
                        if (iter2 == iter1->second.begin() && StartTmp == 0)
                        {
                            // ��ʼ���ֵ�
                            if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                            {
                                chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                            }

                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(1, StartTmp - 1));  // ��ֵ
                        }
                        else if (iter2Tmp == iter1->second.end())  // ��������һ���ڵ㣬�ж��Ƿ񵽴�Ⱦɫ���β
                        {
                            // Ⱦɫ�峤��
                            uint32_t noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);

                            if (EndTmp < noSynEndTmp)
                            {
                                // ��ʼ���ֵ�
                                
                                if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                                {
                                    chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                                }

                                chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(EndTmp + 1, noSynEndTmp));  // ��ֵ
                            }
                        }
                        
                        // ���¸���Ʒ��û�ҵ�
                        sampleNoSynFindBoolMap[sampleTmp] = false;
                    }
                }
                // ���·ǹ����Ե�����
                preNoSynNum = noSynNum;
            }
        }


        // ���תΪ�ַ���
        for (auto iter1 : chrStartSampleLociVecMap)  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
        {
            for (auto iter2 : iter1.second)  // map<refStart, map<sample, vector<tuple<start, end> > > >
            {
                // �յĽڵ�����
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