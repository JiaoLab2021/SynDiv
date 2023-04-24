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
    // ���Դ���
    bool debug = false;
    // bool debug = true;

    // �߳���
    int threads = 30;

    // kseq.h ���ļ�
    KSEQ_INIT(gzFile, gzread)



    /**
     * @brief ����ȫ�ֲ���
     * 
     * @param debugTmp     �Ƿ���Դ���
     * @param threadsTmp   �߳���
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

        // ������Դ��룬�߳�������Ϊ1
        if (debug)
        {
            threads = 1;
        }
        
        return 0;
    }


    /**
     * @brief ��������
     * 
     * @param alignsTitles     aligns �ļ� title
     * @param alignsVec        aligns ·��Vec
     * @param syriConfigFileName   syri ����������ļ�
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

        uint32_t indexTmp = 0;  // ��¼ title ������
        for (auto iter1 : alignsVec)
        {
            string alignsPathTmp = iter1;
            string ailgnsTitleTmp;

            // ������� title
            if (alignsTitles.size() > 0)
            {
                ailgnsTitleTmp = alignsTitles[indexTmp];
            }
            else
            {
                vector<string> alignsPathVecTmp = split(alignsPathTmp, "/");  // ·�����
                ailgnsTitleTmp = alignsPathVecTmp[alignsPathVecTmp.size() - 1];  // ���һ��
                std::regex reg1(".aligns");  // �滻
                ailgnsTitleTmp = regex_replace(ailgnsTitleTmp, reg1, "");  // 'An-1.aligns' ɾ�� '.aligns'
            }

            // ��ֵ
            alignsMap[ailgnsTitleTmp] = alignsPathTmp;

            indexTmp++;  // ��������
        }

        /* ************************************* syri.out ************************************* */
        map<string, map<string, string>> sampleSampleSyriMap;  // map<sample1, map<sample2, syriOutPath> >

        // open file
        GzChunkReader GzChunkReaderClass(syriConfigFileName);

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

            // ����г��Ƿ�Ϊ3
            if (lineVec.size() != 3)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: The '" << syriConfigFileName << "' dataframe does not contain three columns." << endl;
                exit(1);
            }

            // ��Ʒ��·����Ϣ
            const string& sample1 = lineVec[0];
            const string& sample2 = lineVec[1];
            const string& syriPath = lineVec[2];

            // ��ʼ��
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

        // ��¼λ�㹲���Ե�����
        map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >

        // ��¼λ�㹲���Զ�Ӧ����Ʒ����λ��
        map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
    };




    // ����ת��
    class SYNCOOR
    {
    private:
        // �����������ļ�  'coor' ������
        string fileName_;

        // �ж� 'sampleNameVec' ���Ƿ������е�sample
        bool sampleNameVecBool = false;
    public:
        // ����������
        map<string, map<int, map<string, tuple<uint32_t, uint32_t> > > > coorChrLociSampleLociMap;  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >

        // ��Ʒ��
        vector<string> sampleNameVec;

        SYNCOOR() {}
        SYNCOOR(string fileName) {
            fileName_ = fileName;
        }
        ~SYNCOOR() {}

        // ��ʼ���ӿڣ����ڼ̳���ʹ��
        void init(string fileName) {
            fileName_ = fileName;
        }

        // ��������
        void open_coor()
        {
            // open file
            GzChunkReader GzChunkReaderClass(fileName_);

            // read line
            string line;
            while (GzChunkReaderClass.read_line(line))
            {
                // ��������
                if (line.empty())
                {
                    continue;
                }

                // �Ʊ�����
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // �ȼ�������ǲ������ı���
                if (lineVec.size() % 3 != 0)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "'Error: The number of columns is not a multiple of three -> " << fileName_ << endl;
                    exit(1);
                }

                // ÿ���м�¼һ��
                int idxTmp = lineVec.size() / 3;
                // ��ʱ�� 
                string chrTmp = "";
                uint32_t refStartTmp = 0;
                for (size_t i = 0; i < idxTmp; i++)
                {
                    string sampleNameTmp = lineVec[3*i + 0];
                    uint32_t qryStartTmp = 0;
                    uint32_t qryEndTmp = 0;

                    // �ж��Ƿ������֣�����һ��Ϊ0
                    if (isdigit(lineVec[3*i + 1][0]))
                    {
                        qryStartTmp = stol(lineVec[3*i + 1]);
                    }
                    if (isdigit(lineVec[3*i + 2][0]))
                    {
                        qryEndTmp = stol(lineVec[3*i + 2]);
                    }
                    
                    // Ⱦɫ���
                    if (i == 0)
                    {
                        chrTmp = sampleNameTmp;
                        refStartTmp = stoi(lineVec[3*i + 1]);

                        // ��ʼ�� map 
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

                    // ��¼ 'sampleName'��ֻ��¼��һ�е�
                    if (!sampleNameVecBool)
                    {
                        sampleNameVec.push_back(sampleNameTmp);
                    }
                    
                    // ��ӵ��ܵĹ�ϣ����
                    coorChrLociSampleLociMap[chrTmp][refStartTmp][sampleNameTmp] = make_tuple(qryStartTmp, qryEndTmp);
                }

                // ֻ��¼��һ�е� 'sampleName'
                sampleNameVecBool = true;
            }
        }
    };


    // ����ת��
    class COORTRANS
    {
    private:
        // show-align ������ 'xx.aligns'
        string aliFileName_;

        // sampleName
        string sampleName;

        // open file
        GzChunkReader* GzChunkReaderClass_;
        bool endBool = false;  // ��¼�ļ��Ƿ������

        // ��¼�Ƿ�Ҫѭ������������
        bool chrBool = true;  // ��¼��Ⱦɫ����Ƿ����Ҫ��Ŀǰ��Ⱦɫ����Ƿ�һ�£���һ������
        bool aliBool = true;  // ��¼��alignment����ʼ����ֹ�Ƿ�����Ҫ��
        bool findChrBool = true;  // ��¼Ҫ�ҵ�Ⱦɫ��ź͵�ǰȾɫ����Ƿ�һ��

        // �����ļ�ʱ���ĵ�ֵ
        // Ⱦɫ���
        string chr = "";
        // alignment�ķ�����ʼ����ֹ
        string aliRefStrand = "";
        uint32_t aliRefStart = 0;
        uint32_t aliRefEnd = 0;
        string aliQryStrand = "";
        uint32_t aliQryStart = 0;
        uint32_t aliQryEnd = 0;
        // �ļ���ָ�������еıȶ���Ϣ
        uint32_t refStart = 0;
        string refSeq = "";
        uint32_t refEnd = 0;
        uint32_t qryStart = 0;
        string qrySeq = "";
        uint32_t qryEnd = 0;
        // ��¼�������� '_ref_seq' �Լ���Ӧ������
        uint32_t refLoci = 0;
        int idxTmp = -1;
        uint32_t qryLoci = 0;

        tuple<string, uint32_t, uint32_t> get_alignment_loc(
            string infoTmp
        );
        int next_loci();  // ��һ���ȶԶε�����
    public:
        COORTRANS() {}
        COORTRANS(string aliFileName) {
            aliFileName_ = aliFileName;
            GzChunkReaderClass_ = new GzChunkReader(aliFileName_, 1024 * 1024 * 10);
        }
        ~COORTRANS() {
            // // �ͷ� GzChunkReader
            // if (GzChunkReaderClass_ != nullptr) {  
            //     delete GzChunkReaderClass_;  // �ͷ�ָ��
            //     GzChunkReaderClass_ = nullptr;  // �������free
            // }
        }
        
        // ������
        uint32_t find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };


    /**
     * @brief ���� '-- BEGIN alignment ' �ֶ�
     * 
     * @param infoTmp   ɾ�� '-- BEGIN alignment ' �� ' | ' �ָ���ֶΣ�'+1 278 - 1703'
     * 
     * @return tuple<strand, start, end>
    */
    tuple<string, uint32_t, uint32_t> COORTRANS::get_alignment_loc(
        string infoTmp
    )
    {
        // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ  ref
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
     * @brief ��һ���ȶԶε�����
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
                // ��������
                if (line.empty())
                {
                    return 0;
                }

                // �µ�Ⱦɫ���
                if (line.find("-- Alignments between") != string::npos)  // �µıȶ�Ⱦɫ��    -- Alignments between chr1 and chr1
                {
                    // ����ַ���
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    string refChr = lineVec[3];  // refȾɫ�����ʱֵ
                    string qryChr = lineVec[5];  // qryȾɫ�����ʱֵ

                    if (refChr != qryChr)
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Warning: Skipped alignment due to chromosomal mismatch -> " << refChr << " != " << qryChr << endl;
                        chrBool = false;  // ��¼Ⱦɫ����Ϲ涨
                    }
                    else
                    {
                        chrBool = true;  // ��¼Ⱦɫ�岻���Ϲ涨
                    }

                    // ���жϷ�ֹ��ʧ��һ�� "-- BEGIN alignment"
                    if (chr.length() == 0)
                    {
                        findChrBool = true;  // ��¼����Ҫ��Ⱦɫ��
                    }
                    
                    // �ȶ�Ⱦɫ����Ϣ����
                    chr = refChr;
                }
                // Ⱦɫ�������Ϲ涨
                else if (line.find("-- BEGIN alignment ") != string::npos && chrBool && findChrBool)  // �µıȶ�����   -- BEGIN alignment [ +1 1078 - 68996 | +1 1 - 67961 ]
                {
                    line = strip(line.erase(0, 21), '\n');  // ɾ�� '-- BEGIN alignment [ ' ��21���ַ�
                    line = line.erase(line.size()-2, 2);  // ɾ�� ' ]' ��2���ַ���λ�����������

                    vector<string> lineVec = split(line, " | ");  // �ָ� '+1 278 - 1703 -1 2148 - 751'

                    // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ  ref
                    tie(aliRefStrand, aliRefStart, aliRefEnd) = get_alignment_loc(
                        lineVec[0]
                    );

                    // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ  qry
                    tie(aliQryStrand, aliQryStart, aliQryEnd) = get_alignment_loc(
                        lineVec[1]
                    );
                }
                // ֻѭ�������ֿ�ͷ�Ͱ���syn����
                // Ⱦɫ�������Ϲ涨���ȶ���ʼ����ֹ������Ϲ涨
                else if (isdigit(line[0]) != 0 && chrBool && findChrBool && aliBool)
                {
                    /* ***************************** ref ***************************** */
                    // ����ַ���
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    // �������� (this)
                    refSeq = lineVec[1];
                    // ��ʱ��ref����
                    string refSeqTmp = refSeq;
                    // ɾ��  '.'  ���㳤��
                    refSeqTmp.erase(remove(refSeqTmp.begin(), refSeqTmp.end(), '.'), refSeqTmp.end());
                    if (aliRefStrand == "+")  // ����ȶ�
                    {
                        refStart = stoi(lineVec[0]);
                        refEnd = refStart + refSeqTmp.size() - 1;
                        refLoci = refStart - 1;  // ˢ������
                    }
                    else  // ����ȶ�
                    {
                        refEnd = stoi(lineVec[0]);
                        refStart = refEnd - refSeqTmp.size() + 1;
                        refLoci = refEnd + 1;  // ˢ������
                    }

                    /* ***************************** qry ***************************** */
                    // ��ȡ��һ��
                    GzChunkReaderClass_->read_line(line);

                    // ����ַ���
                    std::istringstream iss1(line);
                    vector<string> infoVecTmp(std::istream_iterator<std::string>{iss1}, std::istream_iterator<std::string>());
                    lineVec.assign(infoVecTmp.begin(), infoVecTmp.end());

                    // �������� (this)
                    qrySeq = lineVec[1];
                    // ��ʱ��ref����
                    string qrySeqTmp = qrySeq;
                    // ɾ��  '.'  ���㳤��
                    qrySeqTmp.erase(remove(qrySeqTmp.begin(), qrySeqTmp.end(), '.'), qrySeqTmp.end());
                    if (aliQryStrand == "+")  // ����ȶ�
                    {
                        qryStart = stoi(lineVec[0]);
                        qryEnd = qryStart + qrySeqTmp.size() - 1;
                        qryLoci = qryStart - 1;  // ˢ������
                    }
                    else  // ����ȶ�
                    {
                        qryEnd = stoi(lineVec[0]);
                        qryStart = qryEnd - qrySeqTmp.size() + 1;
                        qryLoci = qryEnd + 1;  // ˢ������
                    }

                    /* ***************************** ˢ������ ***************************** */
                    idxTmp = -1;
                }
            }
            else  // �ļ�������
            {
                // ȫ���������
                // �����ļ�ʱ���ĵ�ֵ
                // Ⱦɫ���
                chr.clear();
                // alignment�ķ�����ʼ����ֹ
                aliRefStrand.clear();
                aliRefStart = 0;
                aliRefEnd = 0;
                aliQryStrand.clear();
                aliQryStart = 0;
                aliQryEnd = 0;
                // �ļ���ָ�������еıȶ���Ϣ
                refStart = 0;
                refSeq.clear();
                refEnd = 0;
                qryStart = 0;
                qrySeq.clear();
                qryEnd = 0;
                // ��¼�������� '_ref_seq' �Լ���Ӧ������
                refLoci = 0;
                idxTmp = -1;
                qryLoci = 0;

                endBool = true;  // ��¼�ļ��Ƿ������
            }
        }

        return 0;
    }

    /**
     * @brief Ѱ��loci
     * 
     * @param chr_  Ҫ�ҵ�Ⱦɫ��
     * @param refLoci_  Ҫ�ҵ�����
     * @return uint32_t  0->������/û�ҵ�, >0->����
    **/
    uint32_t COORTRANS::find_loci(
        string chr_, 
        uint32_t refLoci_
    )
    {
        /* ***************************** �����꣬���� ��0�� ***************************** */
        if (endBool)
        {
            return 0;
        }
        
        /* ***************************** û�б����� ***************************** */
        // �ж�Ⱦɫ���ǲ�����Ҫ��
        while (chr != chr_ && !endBool)
        {
            findChrBool = false;  // ��¼������Ҫ��Ⱦɫ��
            next_loci();  // ��һ��
        }
        findChrBool = true;  // ��¼����Ҫ��Ⱦɫ��
        // ��������꣬���� '0'
        if (endBool)
        {
            return 0;
        }

        // �ж�alignment��ֹ�ǲ���С��  'refLoci_'
        while (aliRefEnd < refLoci_ && !endBool)  // �ϱ��Ѿ��ж���Ⱦɫ���Ƿ���������Դ˴������ж�
        {
            aliBool = false;  // ��¼����Ҫ�� align
            next_loci();  // ��һ��
        }
        aliBool = true;  // ��¼��Ҫ�� align
        // ��������꣬���� '0'
        if (endBool)
        {
            return 0;
        }

        // �ļ�ָ��ָ����Ҫ����
        while (refEnd < refLoci_ && !endBool)  // �ϱ��Ѿ��ж���Ⱦɫ���Ƿ���������Դ˴������ж�
        {
            next_loci();  // ��һ��
        }
        // ��������꣬���� '0'
        if (endBool)
        {
            return 0;
        }

        // �ҵ���Ҫ����
        if (refStart <= refLoci_ && refLoci_ <= refEnd and refLoci < refLoci_)
        {
            for (size_t idx = idxTmp + 1; idx < refSeq.length(); idx++)
            {
                // ref
                if (refSeq[idx] != '.')
                {
                    if (aliRefStrand == "+")  // ����
                    {
                        refLoci++;
                    }
                    else  // ����
                    {
                        refLoci--;
                    }
                }

                // qry
                if (qrySeq[idx] != '.')
                {
                    if (aliQryStrand == "+")  // ����
                    {
                        qryLoci++;
                    }
                    else  // ����
                    {
                        qryLoci--;
                    }
                }

                idxTmp = idx;  // ˢ������

                // �ж��ǲ���Ҫ�ҵ�����
                if (refLoci == refLoci_)
                {
                    break;  // ����ѭ��
                }
            }
        }
        else if (refStart <= refLoci_ && refLoci_ <= refEnd and refLoci == refLoci_)
        {
            return qryLoci;  // ��������
        }
        else  // û�ҵ�
        {
            return 0;
        }

        return qryLoci;
    }


    // ��ȡ����������
    class SYRIOUT
    {
    private:
        // syri.out
        string fileName_ = "";

        // open file
        GzChunkReader* GzChunkReaderClass_;
        bool endBool = false;  // ��¼�ļ��Ƿ������

        // ��¼��������Ϣ
        string chr = "";
        uint32_t refStart = 0;
        uint32_t refEnd = 0;
        uint32_t qryStart = 0;
        uint32_t qryEnd = 0;

        // ˽�к���
        int next_loci();  // ��һ���ȶԶε�����
    public:
        SYRIOUT() {}
        SYRIOUT(string fileName) {
            fileName_ = fileName;
            GzChunkReaderClass_ = new GzChunkReader(fileName_, 1024 * 1024 * 10);
        }
        ~SYRIOUT() {
            // // �ͷ� GzChunkReader
            // if (GzChunkReaderClass_ != nullptr) {  
            //     delete GzChunkReaderClass_;  // �ͷ�ָ��
            //     GzChunkReaderClass_ = nullptr;  // �������free
            // }
        }

        // void open_syri();

        int find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };

    /**
     * @brief ��һ������������
     * 
     * @return 0
    **/
    int SYRIOUT::next_loci()
    {
        // û�б�����
        if (!endBool)
        {
            // read line
            string line;
            if (GzChunkReaderClass_->read_line(line))
            {
                // �������л��߲��� SYNAL ����
                if (line.empty() || line.find("SYNAL") == string::npos)
                {
                    return 0;
                }

                // ���
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // ��������Ƿ�����
                if (lineVec.size() != 12)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "'" << fileName_ << "': does not contain 12 columns. -> " << line << endl;
                    exit(1);
                }

                // ��ֵ
                chr = lineVec[0];
                refStart = stoi(lineVec[1]);
                refEnd = stoi(lineVec[2]);
                qryStart = stoi(lineVec[6]);
                qryEnd = stoi(lineVec[7]);
            }
            else  // �ļ�������
            {
                // ȫ���������
                chr.clear();
                refStart = 0;
                refEnd = 0;
                qryStart = 0;
                qryEnd = 0;

                endBool = true;  // ��¼�ļ�������
            }
        }
        return 0;
    }

    /**
     * @brief Ѱ��loci
     * 
     * @param chr_  Ҫ�ҵ�Ⱦɫ��
     * @param refLoci_  Ҫ�ҵ�����
     * 
     * @return int  -1->�����꣬0->û�ҵ�, 1->����
    **/
    int SYRIOUT::find_loci(
        string chr_, 
        uint32_t refLoci_
    )
    {
        /* ***************************** �����꣬���� ��-1�� ***************************** */
        if (endBool)
        {
            return -1;
        }
        
        /* ***************************** û�б����� ***************************** */
        // �ж�Ⱦɫ���ǲ�����Ҫ��
        while (chr != chr_ && !endBool)
        {
            next_loci();  // ��һ��
        }
        // ��������꣬���� '-1'
        if (endBool)
        {
            return -1;
        }

        // �ж�refEnd�ǲ���С��  'refLoci_'
        while (refEnd < refLoci_ && !endBool)  // �ϱ��Ѿ��ж���Ⱦɫ���Ƿ���������Դ˴������ж�
        {
            next_loci();  // ��һ��
        }
        // ��������꣬���� '-1'
        if (endBool)
        {
            return -1;
        }

        // �ڹ���������
        if (refStart <= refLoci_ && refLoci_ <= refEnd)
        {
            return 1;
        }
        else  // ����������
        {
            return 0;
        }

        return 0;
    }



    /**
     * @brief �����ο�����������
     * 
     * @param referenceFileName -> refgenome
     * 
     * @return refLenMap  map<chr, length>
    **/
    map<string, uint32_t> build_reference_index(
        const string & referenceFileName
    )
    {
        map<string, uint32_t> refLenMap;  // �ο������鳤����Ϣ

        // �����ļ���
        gzFile gzfp = gzopen(referenceFileName.c_str(), "rb");

        // ���ļ�
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

            // �ͷ��ڴ棬�ر��ļ�
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        return refLenMap;
    }




    /* ************************************** calculate syntenic diversity ************************************** */

    /**
     * @brief ����
     * 
     * @param sampleName               ��Ʒ��
     * @param refLenMap                �ο������鳤��
     * @param SynCoorTmp               ����������
     * @param sampleSampleSyriMap      syri.out ���·���ֵ�
     * @param alignsMap                show-aligns ���·���ֵ�
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
        // ������
        CALSTRUCTURE CALSTRUCTURETMP;

        // ��Ʒ��
        CALSTRUCTURETMP.sampleName = sampleName;

        /* ************************************************ Find the index of sampleName ************************************************ */
        // �� sampleName ��һ������
        vector<string>::const_iterator findResult = find(SynCoorTmp.sampleNameVec.begin(), SynCoorTmp.sampleNameVec.end(), sampleName);
        int idxTmp = distance(SynCoorTmp.sampleNameVec.begin(), findResult) + 1;

        // ��������һ��������ֱ������
        if (idxTmp > SynCoorTmp.sampleNameVec.size() - 1)
        {
            return CALSTRUCTURETMP;
        }
        
        /* ************************************************ Syntenic Coordinate ************************************************ */
        // ����������  class��map
        map<string, SYRIOUT> SyriOutMap;  // map<sampleName, SYRIOUT>
        if (sampleName != "reference")  // �ο����������ıȽϲ���Ҫ���������β�ˣ�Ҳ����
        {
            // ��sample�ĵ�����
            map<string, map<string, string> >::const_iterator findIter1 = sampleSampleSyriMap.find(sampleName);
            if (findIter1 == sampleSampleSyriMap.end()) // ���û���ύ����Ʒ��Ӧ�� syri.out �ļ�������
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: '" << sampleName << "' is not present in sampleSampleSyriMap." << endl;
                exit(1);
            }

            // ��������������ʼ�ж�
            for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
            {
                string sampleName2 = SynCoorTmp.sampleNameVec[i];

                if (sampleName2 == "reference")  // �ο�������û�� syri.out
                {
                    continue;
                }

                // ��sample2�ĵ�����
                map<string, string>::const_iterator findIter2 = findIter1->second.find(sampleName2);
                if (findIter2 == findIter1->second.end()) // ���û���ύ����Ʒ��Ӧ�� syri.out �ļ�������
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Both '" << sampleName << "' and '" << sampleName2 << "' are not present in sampleSampleSyriMap." << endl;
                    exit(1);
                }
                string syriOutPath = findIter2->second;

                // ��ʱ�ṹ��
                SyriOutMap[sampleName2] = SYRIOUT(syriOutPath);
            }
        }

        /* ************************************************ Coordinate transformation ************************************************ */
        // ����ת�� class��map
        map<string, COORTRANS> CoorTransMap;  // map<sampleName, COORTRANS>
        for (auto iter1 : SynCoorTmp.sampleNameVec)  // vector<sampleName>
        {
            if (iter1 == "reference")  // �ο����������ıȽϲ���Ҫת��
            {
                continue;
            }

            // ���� aligns ·��
            map<string, string>::const_iterator findIter1 = alignsMap.find(iter1);
            if (findIter1 == alignsMap.end()) // ���û���ύ����Ʒ��Ӧ�� aligns �ļ�������
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: '" << iter1 << "' cannot be found in alignsMap." << endl;
                exit(1);
            }
            
            // ��ʱ�ṹ��
            CoorTransMap[iter1] = COORTRANS(findIter1->second);
        }

        /* ************************************************ Calculate syntenic diversity ************************************************ */
        // ��¼ÿ��Ⱦɫ�����һ�����꣬�����ж��Ƿ񵽴�Ⱦɫ��ĩβ
        map<string, uint32_t> refChrMapTmp;  // map<chr, length>

        // �����ܵ������ϣ��
        for (const auto& iter1 : SynCoorTmp.coorChrLociSampleLociMap)  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
        {
            // Ⱦɫ���
            string chrTmp = iter1.first;

            // ��ʼ���ֵ�
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

            // ��¼��һ������
            uint32_t refEndTmp = 0;

            // ��¼��һ�� refLoci �Ĺ����Ե÷֣����һ�����򲻴洢�������ڴ�����
            uint32_t preSynNum = 0;
            map<string, uint32_t> preSynNumMap;
            for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
            {
                string sampleName2 = SynCoorTmp.sampleNameVec[i];
                preSynNumMap[sampleName2] = 0;
            }

            for (const auto& iter2 : iter1.second)  // map<refStart, map<sample, tuple(start, end)> >
            {
                // ���е���ʼ����ֹ
                uint32_t refStart = get<0>(iter2.second.at("reference"));
                uint32_t refEnd = get<1>(iter2.second.at("reference"));

                // ��¼�����Ե���ֹ����
                refChrMapTmp[chrTmp] = refEnd;

                // ��ʱ����map������sample�ģ����ٲ�ѯ��ʱ��
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

                // �������һ�������м��
                if (refStart - refEndTmp > 1)
                {
                    for (size_t refLoci = refEndTmp + 1; refLoci < refStart; refLoci++)
                    {
                        uint32_t synNumTmp = 0;  // ��������������ʱ

                        for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
                        {
                            string sampleName2 = SynCoorTmp.sampleNameVec[i];

                            // �����reference����û���꣬�Ǿ��Ƕ�û����
                            if (sampleName == "reference")
                            {
                                continue;
                            }
                            
                            // ���������
                            uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // �ο������������ת��sampleName
                            uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // �ο������������ת��sampleName2

                            if (lociA == 0 && lociB == 0) // both del
                            {
                                synNumTmp++;  // ��Ϊ syn

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
                            else if (lociA != 0 && lociB == 0)  // һ�������꣬��һ��û�У�one del
                            {
                                continue;
                            }
                            else if (lociA == 0 && lociB != 0)  // һ�������꣬��һ��û�У�one del
                            {
                                continue;
                            }
                            else  // �������꣬���ж������A��B�ǲ��ǹ�����
                            {
                                if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A��B�ǹ�����
                                {
                                    synNumTmp++;  // ��Ϊ syn

                                    // debug
                                    if (debug)
                                    {
                                        sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";

                                        cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                    }
                                }
                            }
                        }

                        // �������һ���� synNum ��һ�������
                        // if (synNumTmp != preSynNum || (synNumTmp == 0 && preSynNum == 0))
                        if (synNumTmp != preSynNum)
                        {
                            sampleSynOutMap[refLoci] += synNumTmp;
                            preSynNum = synNumTmp;
                        }
                    }
                }

                // ����ȶ�������ĵ÷�
                for (size_t refLoci = refStart; refLoci < refEnd + 1; refLoci++)
                {
                    uint32_t synNumTmp = 0;  // ��������������ʱ

                    for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
                    {
                        string sampleName2 = SynCoorTmp.sampleNameVec[i];

                        // ����� 'coor.txt' �������꣬�����ǹ����ԣ�ֱ�Ӽ�¼
                        // SampleLociMapTmp  ->  // map<sample, tuple<start, end> >
                        if (
                            get<0>(SampleLociMapTmp[sampleName]) > 0 && 
                            get<0>(SampleLociMapTmp[sampleName2]) > 0
                        )  // both syntenic
                        {
                            synNumTmp++;  // ��Ϊ syn

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
                            // �����reference����û���꣬�Ǿ��Ƕ�û����
                            if (sampleName == "reference")
                            {
                                continue;
                            }

                            // ���������
                            uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // �ο������������ת��sampleName
                            uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // �ο������������ת��sampleName2

                            if (lociA == 0 && lociB == 0) // both del
                            {
                                synNumTmp++;  // ��Ϊ syn
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";

                                // debug
                                if (debug)
                                {
                                    cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                }
                            }
                            else if (lociA != 0 && lociB == 0)  // һ�������꣬��һ��û�У�one del
                            {
                                continue;
                            }
                            else if (lociA == 0 && lociB != 0)  // һ�������꣬��һ��û�У�one del
                            {
                                continue;
                            }
                            else  // �������꣬���ж������Ӧ��A��B�ǲ��ǹ�����
                            {
                                if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A��B�ǹ�����
                                {
                                    synNumTmp++;  // ��Ϊ syn

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

                    // �������һ���� synNum ��һ�������
                    // if (synNumTmp != preSynNum || (synNumTmp == 0 && preSynNum == 0))
                    if (synNumTmp != preSynNum)
                    {
                        sampleSynOutMap[refLoci] += synNumTmp;
                        preSynNum = synNumTmp;
                    }
                }

                // ��������
                refEndTmp = refEnd;
            }
        }

        // �ж��Ƿ񵽴�Ⱦɫ��ĩβ
        for (const auto& iter1 : refChrMapTmp)  // map<chr, length>
        {
            string chrTmp = iter1.first;  // Ⱦɫ��
            uint32_t refPosTmp = iter1.second;  // ���һ������

            // ��ʼ���ֵ�
            auto& sampleSynOutMap = CALSTRUCTURETMP.sampleSynOutMap[chrTmp];
            auto& sampleSynOutMapTmp =  CALSTRUCTURETMP.sampleSynOutMapTmp[chrTmp];
            
            // ���ο����������Ƿ��ж�Ӧ��Ⱦɫ��
            map<string, uint32_t>::const_iterator findIter1 = refLenMap.find(chrTmp);
            if (findIter1 == refLenMap.end())
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The reference genome does not include '" << chrTmp << "'." << endl;
                exit(1);
            }

            uint32_t refLen = findIter1->second;  // Ⱦɫ�峤��

            // �������ĩβ������ loop
            if (refPosTmp >= refLen)
            {
                break;
            }
            
            // ���û����Ⱦɫ��ĩβ
            uint32_t preSynNum = 0;
            for (uint32_t refLoci = refPosTmp + 1; refLoci < refLen + 1; refLoci++)
            {
                uint32_t synNumTmp = 0;  // ��������������ʱ

                for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
                {
                    string sampleName2 = SynCoorTmp.sampleNameVec[i];

                    // �����reference����û���꣬�Ǿ��Ƕ�û����
                    if (sampleName == "reference")
                    {
                        continue;
                    }
                    
                    // ���������
                    uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // �ο������������ת��sampleName
                    uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // �ο������������ת��sampleName2

                    if (lociA == 0 && lociB == 0) // both del
                    {
                        synNumTmp++;  // ��Ϊ syn

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
                    else if (lociA != 0 && lociB == 0)  // һ�������꣬��һ��û�У�one del
                    {
                        continue;
                    }
                    else if (lociA == 0 && lociB != 0)  // һ�������꣬��һ��û�У�one del
                    {
                        continue;
                    }
                    else  // �������꣬���ж������A��B�ǲ��ǹ�����
                    {
                        if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A��B�ǹ�����
                        {
                            synNumTmp++;  // ��Ϊ syn

                            // debug
                            if (debug)
                            {
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                cerr << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                            }
                        }
                    }
                }

                // �������һ���� synNum ��һ�������
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
     * @brief �ϲ����
     * 
     * @param chromosome                         Ⱦɫ���
     * @param chrLen                             Ⱦɫ�峤��
     * @param CALSTRUCTUREVec                    ���߳�������
     * 
     * @return tuple<chromosome, synOutMapTmp>   map<refLoci, synNum>
    **/
    tuple<string, map<uint32_t, uint32_t> > merge(
        string chromosome, 
        uint32_t chrLen, 
        const vector<CALSTRUCTURE> & CALSTRUCTUREVec
    )
    {
        // ������
        map<uint32_t, uint32_t> synOutMapTmp;  //  map<refLoci, synNum>
        
        for (const auto& CALSTRUCTURETmp : CALSTRUCTUREVec) // vector<CALSTRUCTURE>
        {
            // Ⱦɫ���Ӧ�ĵ�����
            map<string, map<uint32_t, uint32_t> >::const_iterator iter0 = CALSTRUCTURETmp.sampleSynOutMap.find(chromosome);
            map<string, map<uint32_t, string> >::const_iterator iter3 = CALSTRUCTURETmp.sampleSynOutMapTmp.find(chromosome);

            // ��ʱ����������¼��ǰ����һ��λ��ĵ÷�
            map<uint32_t, uint32_t>::const_iterator iter1;
            map<uint32_t, uint32_t>::const_iterator iter2;
            
            if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end())  // ������и�Ⱦɫ��
            {
                iter1 = iter0->second.begin();
                iter2 = iter0->second.begin();
                iter2++;
            }
            
            // ��¼��һ���ڵ�� synNum
            uint32_t preSynNumTmp = 0;

            // ѭ��Ⱦɫ�峤�ȣ�ÿ��λ�õ������
            for (uint32_t refLoci = 1; refLoci < chrLen + 1; refLoci++)
            {
                // �ж���ǰλ���Ƿ�ָ���һ��������
                if (refLoci == iter1->first && iter0 != CALSTRUCTURETmp.sampleSynOutMap.end())  // ����ǣ���ֵ�õ�����������
                {
                    preSynNumTmp = iter1->second;
                }
                else if (refLoci == iter2->first && iter0 != CALSTRUCTURETmp.sampleSynOutMap.end() && iter2 != iter0->second.end())  // �ߵ���һ���ڵ㣬���µ�����
                {
                    iter1++;  // ���µ�����
                    iter2++;  // ���µ�����

                    preSynNumTmp = iter1->second;  // ��������
                }

                // λ��Ƶ�ʵ���
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
     * @brief ����
     * 
     * @param refLenMap                �ο������鳤��
     * @param SynCoorTmp               ����������
     * @param sampleSampleSyriMap      syri.out ���·���ֵ�
     * @param alignsMap                show-aligns ���·���ֵ�
     * @param allSynNum                �����������
     * @param outputFileName           ����ļ���
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
        // ���̳�
        ThreadPool pool(threads);

        // ��ʼ���̳߳�
        pool.init();

        // ������̵߳Ľ��
        /*
            map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >
            map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
        */
        vector<future<CALSTRUCTURE> > calOutStrVecTmp;
        
        // ��Ʒ����
        for (const auto& iter1 : SynCoorTmp.sampleNameVec)  // vector<sampleName>
        {
            string sampleNameTmp = iter1;

            // ���һ����������
            if (sampleNameTmp == SynCoorTmp.sampleNameVec[SynCoorTmp.sampleNameVec.size() - 1])
            {
                continue;
            }

            cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate: " << sampleNameTmp << endl;

            // ���߳��ύ��������
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

            sleep(0.1);  // �̼߳��
        }

        // ���߳̽������
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

            // ���߳��ύ��������
            mergeOutVec.push_back(
                pool.submit(
                    merge, 
                    chromosome, 
                    chrLen, 
                    ref(CALSTRUCTUREVecTmp)
                )
            );

            sleep(0.1);  // �̼߳��
        }
        

        /* ********************************************** save the result ********************************************** */
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Save the result." << endl;
        SAVE::SAVE SAVEClass(outputFileName);

        stringstream outStream; // ʹ�� stringstream �����ַ���ƴ��
        constexpr uint64_t CACHE_SIZE = 1024 * 1024 * 10; // �����СΪ 10mb

        for (size_t i = 0; i < mergeOutVec.size(); i++)  // vector<future<tuple<string, map<uint32_t, uint32_t> > > > 
        {
            string chomosome;
            map<uint32_t, uint32_t> calOutMap;  // map<uint32_t, uint32_t>

            tie(chomosome, calOutMap) = move(mergeOutVec[i].get());  // ��ȡ���߳̽��

            for (const auto& iter1 : calOutMap)
            {
                outStream << chomosome << "\t" << iter1.first << "\t" << allSynNum << "\t" 
                        << iter1.second << "\t" << iter1.second/(float)allSynNum << "\n";

                if (outStream.tellp() >= CACHE_SIZE) {  // �����СΪ 10mb
                    string outTxt = outStream.str();
                    SAVEClass.save(outTxt);
                    // ��� stringstream
                    outStream.str(string());
                    outStream.clear();
                }
            }

            if (outStream.tellp() >= 0) {  // ���дһ��
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // ��� stringstream
                outStream.str(string());
                outStream.clear();
            }
        }

        // �ر��̳߳�
        pool.shutdown();

        return 0;
    }
}

#endif