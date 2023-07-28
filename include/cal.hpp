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
#include <memory>
#include <string.h>
#include <malloc.h>
#include <getopt.h>

#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "GzChunkReader.hpp"
#include "save.hpp"
#include "ThreadPool.hpp"


using namespace std;


// ���Դ���
extern bool debugCal;

// �߳���
extern int threadsCal;

// ��ȡ�ļ��Ļ����С
extern int32_t readBuffer;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

void help_cal(char* argv[]);
int main_cal(int argc, char* argv[]);


namespace CALNAME
{
    // kseq.h ���ļ�
    KSEQ_INIT(gzFile, gzread)

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
    );



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
        unique_ptr<GzChunkReader> GzChunkReaderClass_;
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

        explicit COORTRANS(string aliFileName) 
            : aliFileName_(aliFileName), GzChunkReaderClass_(std::make_unique<GzChunkReader>(aliFileName_, 1024 * 1024 * readBuffer)) {}

        COORTRANS(const COORTRANS& other) = delete;  // ��ֹ�������캯��
        COORTRANS& operator=(const COORTRANS& other) = delete;  // ��ֹ������ֵ�����

        // �ƶ����캯��
        COORTRANS(COORTRANS&& other) noexcept : 
            aliFileName_(std::move(other.aliFileName_)), 
            sampleName(std::move(other.sampleName)),
            GzChunkReaderClass_(std::move(other.GzChunkReaderClass_)),
            endBool(other.endBool),
            chrBool(other.chrBool),
            aliBool(other.aliBool),
            findChrBool(other.findChrBool),
            chr(std::move(other.chr)),
            aliRefStrand(std::move(other.aliRefStrand)),
            aliRefStart(other.aliRefStart),
            aliRefEnd(other.aliRefEnd),
            aliQryStrand(std::move(other.aliQryStrand)),
            aliQryStart(other.aliQryStart),
            aliQryEnd(other.aliQryEnd),
            refStart(other.refStart),
            refSeq(std::move(other.refSeq)),
            refEnd(other.refEnd),
            qryStart(other.qryStart),
            qrySeq(std::move(other.qrySeq)),
            qryEnd(other.qryEnd),
            refLoci(other.refLoci),
            idxTmp(other.idxTmp),
            qryLoci(other.qryLoci)
        {
            // ��ԭ����״̬��Ϊ��Ч
            other.GzChunkReaderClass_ = nullptr;
        }

        // �ƶ���ֵ�����
        COORTRANS& operator=(COORTRANS&& other) {
            if (this != &other) {
                aliFileName_ = std::move(other.aliFileName_);
                sampleName = std::move(other.sampleName);
                GzChunkReaderClass_ = std::move(other.GzChunkReaderClass_);
                endBool = other.endBool;
                chrBool = other.chrBool;
                aliBool = other.aliBool;
                findChrBool = other.findChrBool;
                chr = std::move(other.chr);
                aliRefStrand = std::move(other.aliRefStrand);
                aliRefStart = other.aliRefStart;
                aliRefEnd = other.aliRefEnd;
                aliQryStrand = std::move(other.aliQryStrand);
                aliQryStart = other.aliQryStart;
                aliQryEnd = other.aliQryEnd;
                refStart = other.refStart;
                refSeq = std::move(other.refSeq);
                refEnd = other.refEnd;
                qryStart = other.qryStart;
                qrySeq = std::move(other.qrySeq);
                qryEnd = other.qryEnd;
                refLoci = other.refLoci;
                idxTmp = other.idxTmp;
                qryLoci = other.qryLoci;
            }
            return *this;
        }
        
        // ������
        uint32_t find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };


    // ��ȡ����������
    class SYRIOUT
    {
    private:
        // syri.out
        string fileName_ = "";

        // open file
        unique_ptr<GzChunkReader> GzChunkReaderClass_;
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
        SYRIOUT() : endBool(true) {}
        explicit SYRIOUT(string fileName) : 
            fileName_(fileName), GzChunkReaderClass_(std::make_unique<GzChunkReader>(fileName_, 1024 * 1024 * readBuffer)), endBool(false) {}

        SYRIOUT(const SYRIOUT& other) = delete;  // ��ֹ�������캯��
        SYRIOUT& operator=(const SYRIOUT& other) = delete;  // ��ֹ������ֵ�����

        // �ƶ����캯��
        SYRIOUT(SYRIOUT&& other) noexcept : 
            fileName_(std::move(other.fileName_)), 
            GzChunkReaderClass_(std::move(other.GzChunkReaderClass_)),
            endBool(other.endBool),
            chr(std::move(other.chr)),
            refStart(other.refStart),
            refEnd(other.refEnd),
            qryStart(other.qryStart),
            qryEnd(other.qryEnd) 
        {
            // ��ԭ����״̬��Ϊ��Ч
            other.GzChunkReaderClass_ = nullptr;
        }

        // �ƶ���ֵ�����
        SYRIOUT& operator=(SYRIOUT&& other) {
            if (this != &other) {
                fileName_ = std::move(other.fileName_);
                GzChunkReaderClass_ = std::move(other.GzChunkReaderClass_);
                endBool = other.endBool;
                chr = std::move(other.chr);
                refStart = other.refStart;
                refEnd = other.refEnd;
                qryStart = other.qryStart;
                qryEnd = other.qryEnd;
            }
            return *this;
        }

        int find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };


    /**
     * @brief �����ο�����������
     * 
     * @param referenceFileName -> refgenome
     * 
     * @return refLenMap  map<chr, length>
    **/
    map<string, uint32_t> build_fasta_index(
        const string & referenceFileName
    );


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
    );


    /* ********************************************* memory-saving mode ********************************************* */

    /**
     * @brief �ϲ����
     * 
     * @param calOutStr                 ĳһ��Ʒ�����м�����
     * @param refLenMap                 Ⱦɫ�峤����Ϣ  map<string, length>
     * @param chrLociSynNumMap          �������յĽ��map<chr, vector<synNum> >
     * 
     * @return 0
    **/
    int merge(
        const CALSTRUCTURE& calOutStr,
        const map<string, uint32_t> & refLenMap,
        map<string, vector<uint32_t> >& chrSynNumVecMap
    );


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
    );



    /* ********************************************* quick mode ********************************************* */

    /**
     * @brief �ϲ����
     * 
     * @param chromosome                         Ⱦɫ���
     * @param chrLen                             Ⱦɫ�峤��
     * @param CALSTRUCTUREVec                    ���߳�������
     * 
     * @return tuple<chromosome, synOutVecTmp>   tuple<chr, vector<synNum> >
    **/
    tuple<string, vector<uint32_t> > merge_fast(
        string chromosome, 
        uint32_t chrLen, 
        const vector<CALSTRUCTURE> & CALSTRUCTUREVec
    );


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
    int calculate_fast(
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap, 
        const uint32_t & allSynNum, 
        const string & outputFileName
    );
}

#endif