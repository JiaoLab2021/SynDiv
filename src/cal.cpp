// g++ -c cal.cpp -o cal -lz -lpthread -O3 -std=c++17
#include "../include/cal.hpp"

using namespace std;

// debug code
bool debugCal = false;

// Threads
int threadsCal = 30;

// Cache size for reading files
int32_t readBuffer = 1;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

void help_cal(char* argv[]);

int main_cal(int argc, char* argv[])
{
    // reference genome
    string referenceFileName;

    // collinear coordinates
    string coorFileName;

    // Configuration file output by syri
    string syriConfigFileName;

    // aligns file with the reference genome
    vector<string> alignsVec;
    vector<string> alignsTitles;
    // Determine whether there is a titlename
    bool haveTitles = false;

    // output filename
    string outputFileName;

    // Whether the software runs in express mode
    bool fastBool = false;

    //Parse command line options
    if(argc <= 2)
    {
        help_cal(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_cal(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-r", 2, parameterLength)) || 
        (PARAMETER_CHECK("--reference", 11, parameterLength))) {
            if ((i+1) < argc) {
                referenceFileName = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("--coor", 6, parameterLength)) {
            if ((i+1) < argc) {
                coorFileName = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("--syri_outs", 11, parameterLength)) {
            if ((i+1) < argc) {
                syriConfigFileName = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("--aligns", 8, parameterLength)) {
            if ((i+1) < argc) {
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    alignsVec.push_back(file);
                    i++;
                    if (i < argc)
                        file = argv[i];
                }
                i--;
            }
        }
        else if((PARAMETER_CHECK("-n", 2, parameterLength)) || 
        (PARAMETER_CHECK("--names", 7, parameterLength))) {
            if ((i+1) < argc) {
                haveTitles = true;
                i = i+1;
                string title = argv[i];
                while (title[0] != '-' && i < argc) {
                    alignsTitles.push_back(title);
                    i++;
                    if (i < argc)
                        title = argv[i];
                }
                i--;
            }
        }
        else if(PARAMETER_CHECK("-o", 2, parameterLength) ||
        (PARAMETER_CHECK("--output", 8, parameterLength))) {
            if ((i+1) < argc) {
                outputFileName = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-t", 2, parameterLength) ||
        (PARAMETER_CHECK("--threads", 9, parameterLength))) {
            if ((i+1) < argc) {
                threadsCal = stoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("--buffer", 8, parameterLength)) {
            if ((i+1) < argc) {
                readBuffer = stoul(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("--fast", 6, parameterLength)) {
            if ((i) < argc) {
                fastBool = true;
            }
        }
        else if(PARAMETER_CHECK("--debug", 7, parameterLength)) {
            if ((i) < argc) {
                debugCal = true;
            }
        }
    }

    if (argc <= 2) {
        help_cal(argv);
        exit(1);
    }
    if (threadsCal <= 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Threads must be greater than 0 (-t)." << endl;
        help_cal(argv);
        exit(1);
    }
    if (readBuffer <= 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: ReadBuffer must be greater than 0 (--buffer)." << endl;
        help_cal(argv);
        exit(1);
    }
    if (referenceFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '-r'." << endl;
        help_cal(argv);
        exit(1);
    }
    if (coorFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '--coor'." << endl;
        help_cal(argv);
        exit(1);
    }
    if (syriConfigFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '--syri_outs'." << endl;
        help_cal(argv);
        exit(1);
    }
    if (alignsVec.size() <= 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Only " << alignsVec.size() << " aligns file was specified. Nothing to extract, exiting." << endl;
        help_cal(argv);
        exit(1);
    }
    if ((haveTitles == true) && (alignsVec.size() != alignsTitles.size())) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of file titles (-n) does not match the number of files (--aligns)." << endl;
        help_cal(argv);
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;


    /* ************************************ Change Global Parameters ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Change Global Parameters ..." << endl;
    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "debug: " << debugCal << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "threads: " << threadsCal << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "readBuffer: " << readBuffer << " MB" << endl;

    // ������Դ��룬�߳�������Ϊ1
    if (debugCal)
    {
        threadsCal = 1;
    }

    /* ************************************ Parse Parameters ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Parse Parameters ..." << endl;
    // ��ȡ aligns/syri.out �ļ��ֵ�
    map<string, string> alignsMap;  // map<sampleName, alignsPath>
    map<string, map<string, string> > sampleSampleSyriMap;  // map<sample1, map<sample2, syriOutPath> >
    tie(alignsMap, sampleSampleSyriMap) = CALNAME::aligns_parameter(
         alignsTitles, 
         alignsVec, 
         syriConfigFileName
    );

    /* ************************************ Building Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building Syntenic Coordinates Index ..." << endl;
    CALNAME::SYNCOOR SynCoorTmp(coorFileName);
    SynCoorTmp.open_coor();

    // �����������
    uint32_t allSynNum = SynCoorTmp.sampleNameVec.size() * (SynCoorTmp.sampleNameVec.size() - 1) / 2;  // n*(n-1)/2

    /* ************************************ Building Reference Genome Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building Reference Genome Index ..." << endl;
    map<string, uint32_t> refLenMap = CALNAME::build_fasta_index(
        referenceFileName
    );

    /* ************************************ Calculating Syntenic Diversity ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating Syntenic Diversity ..." << endl;
    if (fastBool)
    {
        CALNAME::calculate_fast(
            refLenMap, 
            SynCoorTmp, 
            sampleSampleSyriMap, 
            alignsMap, 
            allSynNum, 
            outputFileName
        );
    }
    else
    {
        CALNAME::calculate(
            refLenMap, 
            SynCoorTmp, 
            sampleSampleSyriMap, 
            alignsMap, 
            allSynNum, 
            outputFileName
        );
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ..." << endl;
    
    return 0;
}


void help_cal(char* argv[])
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -r FILE --coor FILE --syri_outs FILE --aligns FILE1 FILE2 .. FILEn [options]" << endl
       << "compute syntenic diversity" << endl
       << endl
       << "required arguments:" << endl
       << "    -r, --reference     FILE      input FASTA reference" << endl
       << "    --coor              FILE      syntenic coordinates, output file of coor" << endl
       << "    --syri_outs         FILE      config file for syri output (format: sample1\tsample2\tsyri.out)" << endl
       << "    --aligns            FILE      list of output files of show-aligns (.aligns), one for multiple mate" << endl
       << endl
       << "optional arguments:" << endl
       << "    -n, --names         STRING    list of names to describe each file in --aligns, one for multiple mate" << endl
       << "    -o, --output        FILE      output results to FILE [stdout]" << endl
       << "    -t, --threads       INT       number of compute threads to use [30]" << endl
       << "    --buffer            INT       buffer size for file reading, measured in MB [1]" << endl
       << "    --fast                        enabling quick mode will increase memory consumption (false by default)" << endl
       << "    --debug                       debug code (false by default)" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
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
tuple<map<string, string>, map<string, map<string, string> > > CALNAME::aligns_parameter(
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


/**
 * @brief ���� '-- BEGIN alignment ' �ֶ�
 * 
 * @param infoTmp   ɾ�� '-- BEGIN alignment ' �� ' | ' �ָ���ֶΣ�'+1 278 - 1703'
 * 
 * @return tuple<strand, start, end>
*/
tuple<string, uint32_t, uint32_t> CALNAME::COORTRANS::get_alignment_loc(
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
int CALNAME::COORTRANS::next_loci()
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
uint32_t CALNAME::COORTRANS::find_loci(
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
        next_loci();  // next line
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
        next_loci();  // next line
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
        next_loci();  // next line
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


/**
 * @brief ��һ������������
 * 
 * @return 0
**/
int CALNAME::SYRIOUT::next_loci()
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
int CALNAME::SYRIOUT::find_loci(
    string chr_, 
    uint32_t refLoci_
)
{
    /* ***************************** After traversing, return '-1' ***************************** */
    if (endBool) {
        return -1;
    }
    
    /* ***************************** not traversed ***************************** */
    // Determine whether the chromosome is desired
    while (chr != chr_ && !endBool) {
        next_loci();  // next line
    }
    // If the traversal is complete, return '-1'
    if (endBool) {
        return -1;
    }

    // Determine whether refEnd is less than 'refLoci_'
    while (refEnd < refLoci_ && !endBool) {  // The above has already judged whether the chromosomes match, so there is no need to judge here
        next_loci();  // next line
    }
    // If the traversal is complete, return '-1'
    if (endBool) {
        return -1;
    }

    // in the collinear interval
    if (refStart <= refLoci_ && refLoci_ <= refEnd) {
        return 1;
    } else {  // not in range
        return 0;
    }

    return 0;
}


/**
 * @brief Build a Reference Genome Index
 * 
 * @param referenceFileName -> refgenome
 * 
 * @return refLenMap  map<chr, length>
**/
map<string, uint32_t> CALNAME::build_fasta_index(
    const string & referenceFileName
)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Compute the length of each chromosome in the reference genome." << endl;

    map<string, uint32_t> refLenMap;  // Reference Genome Length Information

    // input file stream
    gzFile gzfp = gzopen(referenceFileName.c_str(), "rb");

    // open file
    if(!gzfp) {

        cerr << "[" << __func__ << "::" << getTime() << "] "
                << "'"
                << referenceFileName 
                << "': No such file or directory or possibly reached the maximum open file limit. You can set 'ulimit -n' to a larger value to continue." 
                << endl;
        exit(1);
    } else {

        kseq_t *ks;
        ks = kseq_init(gzfp);
    
        while( kseq_read(ks) > 0 ) {
            
            string chromosome = ks->name.s;
            uint32_t chrLen = ks->seq.l;
            // string sequence = ks->seq.s;

            refLenMap[chromosome] = chrLen;
        }

        // free memory, close file
        kseq_destroy(ks);
        gzclose(gzfp);
    }

    return refLenMap;
}




/* ************************************** calculate syntenic diversity ************************************** */

/**
 * @brief ����
 * 
 * @param sampleName               sample name
 * @param refLenMap                chromosome length
 * @param SynCoorTmp               collinear coordinates
 * @param sampleSampleSyriMap      syri.out output path dictionary
 * @param alignsMap                show-aligns output path dictionary
 * 
 * @return CALSTRUCTURE
**/
CALNAME::CALSTRUCTURE CALNAME::calculate_run(
    const string sampleName, 
    const map<string, uint32_t> & refLenMap, 
    const SYNCOOR & SynCoorTmp, 
    const map<string, map<string, string> > & sampleSampleSyriMap, 
    const map<string, string> & alignsMap
)
{
    // save result
    CALSTRUCTURE CALSTRUCTURETMP;

    // Sample name
    CALSTRUCTURETMP.sampleName = sampleName;

    /* ************************************************ Find the index of sampleName ************************************************ */
    // Find the next index after sampleName
    vector<string>::const_iterator findResult = find(SynCoorTmp.sampleNameVec.begin(), SynCoorTmp.sampleNameVec.end(), sampleName);
    int idxTmp = distance(SynCoorTmp.sampleNameVec.begin(), findResult) + 1;

    // If it is the last sample, just skip it
    if (idxTmp > SynCoorTmp.sampleNameVec.size() - 1) {
        return CALSTRUCTURETMP;
    }
    
    /* ************************************************ Syntenic Coordinate ************************************************ */
    // The map of the collinear coordinate class
    map<string, SYRIOUT> SyriOutMap;  // map<sampleName, SYRIOUT>

    if (sampleName != "reference") {  // It is not necessary to compare the reference genome with others. If it is at the end, skip it
        // Find the sample iterator
        map<string, map<string, string> >::const_iterator findIter1 = sampleSampleSyriMap.find(sampleName);
        if (findIter1 == sampleSampleSyriMap.end()) {  // Warn if no syri.out file for this sample is submitted

            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: '" << sampleName << "' is not present in sampleSampleSyriMap." << endl;
        }

        // Judging from the post-index sample
        for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++) {
            string sampleName2 = SynCoorTmp.sampleNameVec[i];

            if (sampleName2 == "reference") {  // Reference genome without syri.out
                continue;
            }

            // Due to too few coordinates or other reasons, there is no alignment result for this sample, and an empty class is constructed
            if (findIter1 == sampleSampleSyriMap.end()) {
                // SYRIOUT
                SyriOutMap[sampleName2] = SYRIOUT();
                continue;
            }
            
            // Find the iterator of sample2
            map<string, string>::const_iterator findIter2 = findIter1->second.find(sampleName2);

            if (findIter2 == findIter1->second.end()) {  // Warn if no syri.out file for this sample is submitted

                cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: '" << sampleName << "' and '" << sampleName2 << "' are not present in sampleSampleSyriMap." << endl;

                // Due to too few coordinates or other reasons, there is no alignment result for this sample, and an empty class is constructed
                SyriOutMap[sampleName2] = SYRIOUT();  // SYRIOUT
                continue;
            }
            string syriOutPath = findIter2->second;

            // SYRIOUT
            SyriOutMap[sampleName2] = SYRIOUT(syriOutPath);
        }
    }

    /* ************************************************ Coordinate transformation ************************************************ */
    // Coordinate conversion class map
    map<string, COORTRANS> CoorTransMap;  // map<sampleName, COORTRANS>

    // Judging from the post-index sample
    for (size_t i = idxTmp - 1; i < SynCoorTmp.sampleNameVec.size(); i++) {
    
        string sampleName2 = SynCoorTmp.sampleNameVec[i];

        if (sampleName2 == "reference") {  // The reference genome has no .aligns
        
            continue;
        }

        // Find the aligns path
        map<string, string>::const_iterator findIter1 = alignsMap.find(sampleName2);
        if (findIter1 == alignsMap.end()) {  // If the aligns file corresponding to the sample is not submitted, an error will be reported
        
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: '" << sampleName2 << "' cannot be found in alignsMap." << endl;
            exit(1);
        }
        
        // ��ʱ�ṹ��
        CoorTransMap[sampleName2] = COORTRANS(findIter1->second);
    }

    /* ************************************************ Calculate syntenic diversity ************************************************ */
    // Record the last coordinate of each chromosome to judge whether it has reached the end of the chromosome
    map<string, uint32_t> refChrMapTmp;  // map<chr, length>

    // Traverse the total coordinate hash table
    for (const auto& iter1 : SynCoorTmp.coorChrLociSampleLociMap) {  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
    
        // chromosome
        string chrTmp = iter1.first;

        // initialize dictionary
        auto& sampleSynOutMap = CALSTRUCTURETMP.sampleSynOutMap[chrTmp];
        auto& sampleSynOutMapTmp =  CALSTRUCTURETMP.sampleSynOutMapTmp[chrTmp];

        // record last coordinate
        uint32_t refEndTmp = 0;

        // Record the collinearity score of the last refLoci, if it is the same, it will not be stored. Reduce memory consumption
        uint32_t preSynNum = 0;

        for(const auto& [refStart, SampleLociMap] : iter1.second) {  // map<refStart, map<sample, tuple(start, end)> >
        
            // start and end of the line
            uint32_t refEnd = get<1>(SampleLociMap.at("reference"));

            // End coordinates for record collinearity
            refChrMapTmp[chrTmp] = refEnd;

            // If there is a distance from the previous coordinate
            if (refStart - refEndTmp > 1) {
            
                for (size_t refLoci = refEndTmp + 1; refLoci < refStart; refLoci++) {
                
                    uint32_t synNumTmp = 0;  // Collinearity Quantity, Temporary

                    for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++) {
                    
                        string sampleName2 = SynCoorTmp.sampleNameVec[i];

                        // If the reference has no coordinates with you, then there are no coordinates
                        if (sampleName == "reference") {
                        
                            continue;
                        }
                        
                        // if other
                        uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // The coordinates of the reference genome go to sampleName
                        uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // The coordinates of the reference genome go to sampleName2

                        if (lociA == 0 && lociB == 0) { // both del
                        
                            synNumTmp++;  // denoted as syn

                            // debug
                            if (debugCal) {
                            
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                            }
                        }
                        else if (lociA != 0 && lociB == 0) {  // One has coordinates, the other doesn't, one del
                        
                            continue;
                        }
                        else if (lociA == 0 && lociB != 0) {  // One has coordinates, the other doesn't, one del
                        
                            continue;
                        }
                        else {  // There are coordinates, and then judge whether A and B of the coordinates are collinear
                        
                            if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A and B are collinear
                            {
                                synNumTmp++;  // denoted as syn

                                // debug
                                if (debugCal) {

                                    sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                    cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                }
                            }
                        }
                    }

                    // �������һ���� synNum ��һ�������
                    // if (synNumTmp != preSynNum || (synNumTmp == 0 && preSynNum == 0))
                    if (synNumTmp != preSynNum)
                    {
                        sampleSynOutMap[refLoci] = synNumTmp;
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
                    // SampleLociMap  ->  // map<sample, tuple<start, end> >
                    if (
                        get<0>(SampleLociMap.at(sampleName)) > 0 && 
                        get<0>(SampleLociMap.at(sampleName2)) > 0
                    )  // both syntenic
                    {
                        synNumTmp++;  // ��Ϊ syn

                        // debug
                        if (debugCal)
                        {
                            sampleSynOutMapTmp[refLoci] += sampleName2 + ";";
                            cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " syn:ture" << endl;
                        }
                    }
                    else if (
                        get<0>(SampleLociMap.at(sampleName)) == 0 && 
                        get<0>(SampleLociMap.at(sampleName2)) > 0
                    )  // one syntenic
                    {
                        continue;
                    }
                    else if (
                        get<0>(SampleLociMap.at(sampleName)) > 0 && 
                        get<0>(SampleLociMap.at(sampleName2)) == 0
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
                            
                            // debug
                            if (debugCal)
                            {
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
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
                                if (debugCal)
                                {
                                    sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                    cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                }
                            }
                        }
                    }
                }

                // �������һ���� synNum ��һ�������
                if (synNumTmp != preSynNum)
                {
                    sampleSynOutMap[refLoci] = synNumTmp;
                    preSynNum = synNumTmp;
                }
            }

            // ��������
            refEndTmp = refEnd;
        }
    }

    // �ж��Ƿ񵽴�Ⱦɫ��ĩβ
    for (const auto& [chrTmp, refPosTmp] : refChrMapTmp)  // map<chr, length>
    {
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
        if (refPosTmp >= refLen) {

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
                    if (debugCal)
                    {
                        sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                        cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
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
                        if (debugCal)
                        {
                            sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                            cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                        }
                    }
                }
            }

            // �������һ���� synNum ��һ�������
            if (synNumTmp != preSynNum)
            {
                sampleSynOutMap[refLoci] = synNumTmp;
                preSynNum = synNumTmp;
            }
        }
    }

    return CALSTRUCTURETMP;
}




/* ********************************************* memory-saving mode ********************************************* */

/**
 * @brief �ϲ����
 * 
 * @param calOutStr                 ĳһ��Ʒ�����м�����
 * @param refLenMap                 Ⱦɫ�峤����Ϣ  map<string, length>
 * @param chrSynNumVecMap           �������յĽ��map<chr, vector<synNum> >
 * 
 * @return 0
**/
int CALNAME::merge(
    const CALSTRUCTURE& calOutStr,
    const map<string, uint32_t> & refLenMap,
    map<string, vector<uint32_t> >& chrSynNumVecMap
)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Combine the computation results of " << calOutStr.sampleName << "." << endl;
    
    for(const auto& [chromosome, chrLen] : refLenMap)  // map<string, length>
    {
        // ������
        auto& SynNumVec = chrSynNumVecMap[chromosome];

        // �����һ�κϲ��Ļ����ȳ�ʼ��SynNumVec�ĳ���ΪȾɫ�峤��
        if (SynNumVec.empty())
        {
            SynNumVec = vector<uint32_t>(chrLen, 0);
        }
        
        // Ⱦɫ���Ӧ�ĵ�����
        map<string, map<uint32_t, uint32_t> >::const_iterator iter0 = calOutStr.sampleSynOutMap.find(chromosome);
        map<string, map<uint32_t, string> >::const_iterator iter3 = calOutStr.sampleSynOutMapTmp.find(chromosome);

        // ��ʱ����������¼��ǰ����һ��λ��ĵ÷�
        map<uint32_t, uint32_t>::const_iterator iter1;
        map<uint32_t, uint32_t>::const_iterator iter2;
        
        if (iter0 != calOutStr.sampleSynOutMap.end())  // ������и�Ⱦɫ��
        {
            iter1 = iter0->second.begin();
            iter2 = iter0->second.begin();
            iter2++;
        }
        else
        {
            continue;
        }
        
        // ��¼��һ���ڵ�� synNum
        uint32_t preSynNumTmp = 0;

        // ѭ��Ⱦɫ�峤�ȣ�ÿ��λ�õ������
        for (uint32_t refLoci = 1; refLoci < chrLen + 1; refLoci++)
        {
            // �ж���ǰλ���Ƿ�ָ���һ��������
            if (iter0 != calOutStr.sampleSynOutMap.end() && refLoci == iter1->first)  // ����ǣ���ֵ�õ�������Ӧ������
            {
                preSynNumTmp = iter1->second;
            }
            else if (iter0 != calOutStr.sampleSynOutMap.end() && iter2 != iter0->second.end() && refLoci == iter2->first)  // �ߵ���һ���ڵ㣬���µ�����
            {
                iter1++;  // ���µ�����
                iter2++;  // ���µ�����

                preSynNumTmp = iter1->second;  // ��������
            }

            // λ��Ƶ�ʵ���
            SynNumVec[refLoci - 1] += preSynNumTmp;

            // debug
            if (debugCal)
            {
                if (preSynNumTmp > 0)
                {
                    cerr << chromosome << " " << refLoci << " " << calOutStr.sampleName << " " << preSynNumTmp << " " << iter3->second.at(refLoci) << endl;
                }
            }
        }
    }

    return 0;
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
int CALNAME::calculate(
    const map<string, uint32_t> & refLenMap, 
    const SYNCOOR & SynCoorTmp, 
    const map<string, map<string, string> > & sampleSampleSyriMap, 
    const map<string, string> & alignsMap, 
    const uint32_t & allSynNum, 
    const string & outputFileName
)
{
    // final result
    map<string, vector<uint32_t> > chrSynNumVecMap;  // map<chr, vector<synNum> >

    /* ********************************************** calculate syntenic diversity ********************************************** */
    // process pool
    ThreadPool pool(threadsCal);
    const int MAX_THREADS_NUM = threadsCal*2;  // The most stored element in a multi-threaded vector, write out when it is greater than

    // init
    pool.init();

    // Save the results of multiple threads
    /*
        map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >
        map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
    */
    vector<future<CALSTRUCTURE> > calOutStrFutureVec;
    
    // sample name
    for (const auto& iter1 : SynCoorTmp.sampleNameVec) {  // vector<sampleName>

        string sampleNameTmp = iter1;

        // last sample skipped
        if (sampleNameTmp == SynCoorTmp.sampleNameVec[SynCoorTmp.sampleNameVec.size() - 1]) {

            continue;
        }

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate: " << sampleNameTmp << endl;

        // Submit and save results in multiple threads
        calOutStrFutureVec.push_back(
            pool.submit(
                calculate_run, 
                sampleNameTmp, 
                ref(refLenMap), 
                ref(SynCoorTmp), 
                ref(sampleSampleSyriMap), 
                ref(alignsMap)
            )
        );

        // If the length of calOutStrVecTmp is greater than the threshold, write it first and then submit the task
        if (calOutStrFutureVec.size() >= MAX_THREADS_NUM) {

            // Traversing the results returned by multiple threads
            for (auto&& calOutStrFuture : calOutStrFutureVec) {  // vector<future<CALSTRUCTURE> >
            
                CALSTRUCTURE calOutStr = move(calOutStrFuture.get());  // future<CALSTRUCTURE>
                // merged result
                merge(calOutStr, refLenMap, chrSynNumVecMap);
            }

            calOutStrFutureVec.clear();
            vector<future<CALSTRUCTURE> >().swap(calOutStrFutureVec);

            malloc_trim(0); // 0 is for heap memory
        }
    }

    // The last combined results
    if (calOutStrFutureVec.size() > 0) {

        // Traversing the results returned by multiple threads
        for (auto&& calOutStrFuture : calOutStrFutureVec) {  // vector<future<CALSTRUCTURE> >
        
            CALSTRUCTURE calOutStr = move(calOutStrFuture.get());  // future<CALSTRUCTURE>
            // merged result
            merge(calOutStr, refLenMap, chrSynNumVecMap);
        }
        calOutStrFutureVec.clear();
        vector<future<CALSTRUCTURE> >().swap(calOutStrFutureVec);

        malloc_trim(0); // 0 is for heap memory
    }

    /* ********************************************** save the result ********************************************** */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;
    // total number of samples
    const uint32_t sampleNum = SynCoorTmp.sampleNameVec.size();
    // Correction factor
    // const double correctionFactor = static_cast<double> (sampleNum - 1) / sampleNum;  // (n-1)/n

    SAVE SAVEClass(outputFileName);

    stringstream outStream;  // Use stringstream instead of string concatenation
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);
    outStream << "#CHROM\tPOS\tAll_States\tSyntenic_States\tSyntenic_Diversity\n";

    for (const auto& [chromosome, SynNumVec] : chrSynNumVecMap)  // map<chr, vector<synNum> >
    {
        uint32_t loci = 1;
        for (const auto& synNum : SynNumVec)  // vector<synNum>
        {
            // outStream << chromosome << "\t" << loci << "\t" << allSynNum << "\t" 
            //         << synNum << "\t" << 1 - (synNum/(double)allSynNum)*(double)correctionFactor << "\n";

            outStream << chromosome << "\t" << loci << "\t" << allSynNum << "\t" 
                    << synNum << "\t" << 1 - synNum/(double)allSynNum << "\n";

            if (outStream.tellp() >= CACHE_SIZE)  // Cache size is 10mb
            {
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // empty stringstream
                outStream.str(string());
                outStream.clear();
            }
            ++loci;
        }

        if (outStream.tellp() > 0)  // write one last time
        {
            string outTxt = outStream.str();
            SAVEClass.save(outTxt);
            // empty stringstream
            outStream.str(string());
            outStream.clear();
        }
    }

    // Close the thread pool
    pool.shutdown();

    return 0;
}



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
tuple<string, vector<uint32_t> > CALNAME::merge_fast(
    string chromosome, 
    uint32_t chrLen, 
    const vector<CALSTRUCTURE> & CALSTRUCTUREVec
)
{
    // ������
    vector<uint32_t> synOutVecTmp(chrLen, 0);  //  vector<synNum>
    
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
            if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end() && refLoci == iter1->first)  // ����ǣ���ֵ�õ�������Ӧ������
            {
                preSynNumTmp = iter1->second;
            }
            else if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end() && iter2 != iter0->second.end() && refLoci == iter2->first)  // �ߵ���һ���ڵ㣬���µ�����
            {
                iter1++;  // ���µ�����
                iter2++;  // ���µ�����

                preSynNumTmp = iter1->second;  // ��������
            }

            // λ��Ƶ�ʵ���
            synOutVecTmp[refLoci - 1] += preSynNumTmp;

            // debug
            if (debugCal)
            {
                if (preSynNumTmp > 0)
                {
                    cerr << chromosome << " " << refLoci << " " << CALSTRUCTURETmp.sampleName << " " << preSynNumTmp << " " << iter3->second.at(refLoci) << endl;
                }
            }
        }
    }

    return make_tuple(chromosome, synOutVecTmp);
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
int CALNAME::calculate_fast(
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
    ThreadPool pool(threadsCal);

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
    vector<future<tuple<string, vector<uint32_t> > > > mergeOutVec;
    for (const auto& iter1 : refLenMap)  // map<chromosome, length>
    {
        string chromosome = iter1.first;
        uint32_t chrLen = iter1.second;

        // ���߳��ύ��������
        mergeOutVec.push_back(
            pool.submit(
                merge_fast, 
                chromosome, 
                chrLen, 
                ref(CALSTRUCTUREVecTmp)
            )
        );

        sleep(0.1);  // �̼߳��
    }
    

    /* ********************************************** save the result ********************************************** */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;
    // �ܵ�������
    const uint32_t sampleNum = SynCoorTmp.sampleNameVec.size();
    // // ����ϵ��
    // const double correctionFactor = static_cast<double> (sampleNum - 1) / sampleNum;  // (n-1)/n

    SAVE SAVEClass(outputFileName);

    stringstream outStream; // ʹ�� stringstream �����ַ���ƴ��
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // �����СΪ 10mb
    outStream.str().reserve(CACHE_SIZE);
    outStream << "#CHROM\tPOS\tAll_States\tSyntenic_States\tSyntenic_Diversity\n";

    for (size_t i = 0; i < mergeOutVec.size(); ++i)  // vector<future<tuple<string, vector<uint32_t> > > > 
    {
        string chromosome;
        vector<uint32_t> calOutVec;  // vector<uint32_t>

        tie(chromosome, calOutVec) = move(mergeOutVec[i].get());  // ��ȡ���߳̽��
        uint32_t loci = 1;
        for (const auto& iter1 : calOutVec)
        {
            // outStream << chromosome << "\t" << iter1.first << "\t" << allSynNum << "\t" 
            //         << iter1.second << "\t" << 1.0 - (iter1.second/(double)allSynNum)*(double)correctionFactor << "\n";

            outStream << chromosome << "\t" << loci << "\t" << allSynNum << "\t" 
                    << iter1 << "\t" << 1.0 - iter1/(double)allSynNum << "\n";

            if (outStream.tellp() >= CACHE_SIZE)  // �����СΪ 10mb
            {
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // ��� stringstream
                outStream.str(string());
                outStream.clear();
            }
            ++loci;
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

    // �ر��̳߳�
    pool.shutdown();

    return 0;
}