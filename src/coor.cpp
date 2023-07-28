// g++ -c src/coor.cpp -O3 -std=c++17
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <future>
#include <malloc.h>
#include "../include/coor.hpp"
#include "../include/ThreadPool.hpp"

using namespace std;

// ȫ�ֱ���
int thresholdLength = 1000;  // ����������ref��qry���ȱ���ֵ

// �Ƿ���Դ���
bool debugCoor = false;


int main_coor(int argc, char* argv[])
{
    // loc ����ļ�
    string synLocFileName;

    // �����ļ��б������  show-aligns
    vector<string> inputFiles;
    vector<string> inputTitles;
    // �ж��Ƿ��� titlename
    bool haveTitles = false;

    // ����ļ���
    string outputFileName;

    // �߳���
    int threads = 30;

    //Parse command line options
    if(argc <= 2)
    {
        help_coor(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_coor(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-s", 2, parameterLength)) || 
        (PARAMETER_CHECK("--syn", 5, parameterLength))) {
            if ((i+1) < argc) {
                synLocFileName = argv[i + 1];
                i++;
            }
        }
        else if((PARAMETER_CHECK("-i", 2, parameterLength)) || 
        (PARAMETER_CHECK("--inputs", 8, parameterLength))) {
            if ((i+1) < argc) {
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    inputFiles.push_back(file);
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
                    inputTitles.push_back(title);
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
                threads = stoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("--debug", 7, parameterLength)) {
            if ((i) < argc) {
                debugCoor = true;
            }
        }
    }

    if (argc <= 2 || inputFiles.size() <= 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Only " << inputFiles.size() << " aligns file was specified. Nothing to extract, exiting." << endl;
        help_coor(argv);
        exit(1);
    }
    if ((haveTitles == true) && (inputFiles.size() != inputTitles.size())) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of file titles (-n) does not match the number of files (-i)." << endl;
        help_coor(argv);
        exit(1);
    }
    if (synLocFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name (-s)." << endl;
        help_coor(argv);
        exit(1);
    }


    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;

    
    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index ..." << endl;
    map<string, vector<tuple<int, int, vector<string> > > > synLocSampleVecMap;  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    synLocSampleVecMap = COOR::build_syn_idx(
        synLocFileName
    );


    /* ************************************ Find Syntenic Coordinates on query genomes ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find Syntenic Coordinates on query genomes ..." << endl;
    ThreadPool pool(threads);  // ���̳�

    // ��ʼ���̳߳�
    pool.init();

    // ������̵߳Ľ��
    vector<future<COOR::synAllStructure> > chrStartSynQryLocMapVec;

    // ��ȡ��������qry�ϵ�����
    int indexTmp = 0;  // ��¼��Ʒ��������
    for (auto it : inputFiles)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "get_syn_coor: " << it << endl;

        string sampleName;
        if (haveTitles)
        {
            sampleName = inputTitles[indexTmp];
        }
        
        // string alignsFileName = it;
        // COOR::synAllStructure chrStartSynQryLocMap;
        // ���߳��ύ��������
        chrStartSynQryLocMapVec.push_back(
            pool.submit(
                COOR::get_syn_coor, 
                sampleName, 
                it, 
                ref(synLocSampleVecMap), 
                false
            )
        );

        indexTmp++;  // ��������
    }

    // ���߳̽������
    map<string, unordered_map<string, unordered_map<int, tuple<int, int> > > > sampleChrStartSynQryLocMap;
    for (size_t i = 0; i < chrStartSynQryLocMapVec.size(); i++)
    {
        COOR::synAllStructure synAllStructureTmp = move(chrStartSynQryLocMapVec[i].get());
        sampleChrStartSynQryLocMap[move(synAllStructureTmp.sampleName)] = move(synAllStructureTmp.chrStartSynQryLocMap);
    }
    chrStartSynQryLocMapVec.clear();
    vector<future<COOR::synAllStructure> >().swap(chrStartSynQryLocMapVec);
    malloc_trim(0);	// 0 is for heap memory

    // �ر��̳߳�
    pool.shutdown();
    
    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;
    COOR::save_result(
        synLocSampleVecMap, 
        sampleChrStartSynQryLocMap, 
        outputFileName
    );

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ..." << endl;

    return 0;
}

void help_coor(char* argv[])
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -s FILE -i FILE1 FILE2 .. FILEn [options]" << endl
       << "retrieve syntenic coordinates on the query genome" << endl
       << endl
       << "required arguments:" << endl
       << "    -s, --syn           FILE      syntenic intersection, output by multiinter" << endl
       << "    -i, --inputs        FILE      list of output files of show-aligns (.aligns), one for multiple mate" << endl
       << endl
       << "optional arguments:" << endl
       << "    -n, --names         STRING    list of names to describe each file in -i, one for multiple mate" << endl
       << "    -o, --output        FILE      output coordinates to FILE [stdout]" << endl
       << "    -t, --threads       INT       number of compute threads to use [30]" << endl
       << "    --debug                       debug code (false by default)" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}



/**
 * @brief ���� '-- BEGIN alignment ' �ֶ�
 * 
 * @param informationTmp    '+1 278 - 1703'
 * 
 * @return tuple<strand, start, end>
*/
tuple<string, int, int> COOR::get_alignment_loc(string informationTmp)
{
    // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ
    istringstream iss(informationTmp);  // +1 278 - 1703
    string strandTmp;
    int startTmp, endTmp;
    char sep;

    iss >> strandTmp >> startTmp >> sep >> endTmp;

    // ���Ϊ����ȶԣ�������ʼ����ֹ����
    if (strandTmp == "-1") {
        swap(startTmp, endTmp);
    }

    return make_tuple(strandTmp, startTmp, endTmp);
}


/**
 * @brief ��qry�Ϲ����Ե�����
 * 
 * @param synLoc        �����Ե�����
 * @param refStart      ref����ʼ
 * @param refEnd        ref����ֹ
 * @param refSeq        ref������
 * @param qryStrand     qry�ȶԷ���
 * @param qryStart      qry����ʼ
 * @param qryEnd        qry����ֹ
 * @param qrySeq        qry������
 * 
 * @return qrySynLoc    qry��syn�����꣬ 0-û�ҵ�
*/
int COOR::find_qry_syn(
    const int & synLoc, 
    int refStart, 
    int refEnd, 
    const string & refSeq, 
    const string & AliQryStrand, 
    int qryStart, 
    int qryEnd, 
    const string & qrySeq
)
{
    int qrySynLoc = 0;

    if (refStart <= synLoc && synLoc <= refEnd)  // ��������˹���������
    {
        int iIdx = 0;  // ��¼synStart��Ӧ��λ��

        --refStart;  // �����ȼ�1
        for (size_t i = 0; i < refSeq.size(); ++i)  // ѭ��ref����
        {
            if (refSeq[i] != '.')
            {
                ++refStart;  // ref�������
            }

            if (refStart == synLoc)  // ��syn���������
            {
                iIdx = i;
                break;  // �ҵ�������˳���ѭ��
            }
        }

        if (AliQryStrand == "+1")  // qry������ȶ�ʱ
        {
            --qryStart;  // �����ȼ�1
            for (size_t i = 0; i < qrySeq.size(); ++i)
            {
                if (qrySeq[i] != '.')
                {
                    ++qryStart;  // ref�������
                }
                
                if (i == iIdx)  // ѭ������λ�ú��˳�
                {
                    break;  // �ҵ�������˳���ѭ��
                }
            }
            qrySynLoc = qryStart;
        }
        else  // qry�Ƿ���ȶ�ʱ
        {
            ++qryEnd;  // �����ȼ�1
            for (size_t i = 0; i < qrySeq.size(); ++i)
            {
                if (qrySeq[i] != '.')
                {
                    --qryEnd;  // ref����ݼ�
                }
                
                if (i == iIdx)  // ѭ������λ�ú��˳�
                {
                    break;  // �ҵ�������˳���ѭ��
                }
            }
            qrySynLoc = qryEnd;
        }
    }

    return qrySynLoc;
}


/**
 * @brief ��qry�Ϲ����Ե�����
 * 
 * @param chrStartSynQryLocMap  synAllStructure
 * @param refChr                Ⱦɫ���
 * @param qrySynStartTmp        find_qry_syn�ҵ�qrySynStart
 * @param qrySynEndTmp          find_qry_syn�ҵ�qrySynEnd
 * @param refStart              �ñȶ���ref����ʼ
 * @param refEnd                �ñȶ���ref����ֹ
 * @param qryStart              �ñȶ���qry����ʼ
 * @param qryEnd                ���ȶ���qry����ֹ
 * @param synStart              ĿǰҪ�ҵ�syn��ref��ʼ
 * @param synEnd                ĿǰҪ�ҵ�syn��ref��ֹ
 * @param qrySynStart           syn��qry��ʼ
 * @param qrySynEnd             syn��qry��ֹ
 * @param whileBool             �ж��Ƿ���Ҫwhileѭ��
 * @param aliRowStartNum        syn��ʼλ���ڸ�alignment�Ĵ�������
 * @param aliRowEndNum          syn��ֹλ���ڸ�alignment�Ĵ�������
 * 
 * @return qrySynLoc    qry��syn�����꣬ 0-û�ҵ�
*/
int COOR::syn_all_loc_push(
    synAllStructure & chrStartSynQryLocMap, 
    const string & refChr, 
    const int & qrySynStartTmp, 
    const int & qrySynEndTmp, 
    const int & refStart, 
    const int & refEnd, 
    const int & qryStart, 
    const int & qryEnd, 
    const int & synStart, 
    const int & synEnd, 
    int & qrySynStart, 
    int & qrySynEnd, 
    bool & whileBool, 
    int & aliRowStartNum, 
    int & aliRowEndNum
)
{
    whileBool = false;

    if (qrySynStartTmp != 0) {
        if (qrySynStart > 0) {  // ֮ǰ�ҵ���qrySynStart��ѡ����preQrySynStart���������
            if (abs(qrySynStart - synStart) > abs(qrySynStartTmp - synStart))
            {
                qrySynStart = qrySynStartTmp;
            } else if (debugCoor) {
                cerr << "skip_start:" << qrySynStartTmp << endl;
            }
        }
        else  // û���ҵ��Ļ�ֱ�Ӹ�ֵ
        {
            qrySynStart = qrySynStartTmp;
        }
        
        aliRowStartNum = aliRowEndNum;  // ����aliRowStartNum������aliRowEndNum����ֹ����
    }
    if (qrySynEndTmp != 0 && qrySynStart != 0)  // ��ʼ�������ҵ�
    {
        if (qrySynEnd > 0)  // ֮ǰ�ҵ���qrySynEnd��ѡ����preQrySynEnd���������
        {
            if (abs(qrySynEnd - synEnd) > abs(qrySynEndTmp - synEnd))
            {
                qrySynEnd = qrySynEndTmp;
            }
            else if (debugCoor)
            {
                cerr << "skip_end:" << qrySynEndTmp << endl;
            }
        }
        else  // û���ҵ��Ļ�ֱ�Ӹ�ֵ
        {
            qrySynEnd = qrySynEndTmp;
        }

        // ����0���ٸ�ֵ���ȶԳ�������С��thresholdLength���ٸ�ֵ
        int refSynLen = abs(synEnd - synStart);
        int qrySynLen = abs(qrySynEnd - qrySynStart);
        if (qrySynEnd > 0 && max(refSynLen, qrySynLen)/(float)min(refSynLen, qrySynLen) < thresholdLength) {
            chrStartSynQryLocMap.chrStartSynQryLocMap[refChr][synStart] = make_tuple(
                min(qrySynStart, qrySynEnd), 
                max(qrySynStart, qrySynEnd)
            );
        }
    }

    // debug
    if (qrySynStartTmp != 0 || qrySynEndTmp != 0)
    {
        if (debugCoor)
        {
            cerr << "ref:" << refChr << "-" << refStart << "-" << refEnd << "\t" << 
                "qry:" << qryStart << "-" << qryEnd << "\t" << 
                "refSyn:" << synStart << "-" << synEnd << "\t" <<
                "qrySyn:" << qrySynStart << "-" << qrySynEnd << "\t" <<
                "synOut:" << qrySynStartTmp << "-" << qrySynEndTmp << endl;
        }
    }

    return 0;
}


/**
 * @brief ������ʱ��syn����
 * 
 * @param sampleName            ��Ʒ���������жϸ����Ƿ�Ҫ�����ж�
 * @param synLocSampleVecMap    map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param synChr                Ⱦɫ��
 * @param synIdx                syn�ڸ�Ⱦɫ���vector������
 * @param synStart              syn����ʼ
 * @param synEnd                syn����ֹ
 * @param qrySynStart           qry��syn����ʼ
 * @param qrySynEnd             qry��syn����ֹ
 * @param whileBool             �ж��Ƿ�����whileѭ��
 * 
 * @return 0
*/
int COOR::renew_syn_loc(
    const string& sampleName, 
    const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap, 
    const string& synChr, 
    int& synIdx, 
    int& synStart, 
    int& synEnd, 
    int& qrySynStart, 
    int& qrySynEnd, 
    bool& whileBool
)
{
    map<string, vector<tuple<int, int, vector<string> > > >::const_iterator findIter = synLocSampleVecMap.find(synChr);  //  �ڹ�����map����Ⱦɫ��
    if (findIter != synLocSampleVecMap.end())  // �ҵ���
    {
        auto& synLocSampleVec = findIter->second;  // vector<tuple<int, int, vector<string> > >

        if (synIdx < synLocSampleVec.size())  // ��ֹԽ��
        {
            // ����ù�����û�ж�Ӧ����Ʒ����������һ������
            while (find(get<2>(synLocSampleVec[synIdx]).begin(), get<2>(synLocSampleVec[synIdx]).end(), sampleName) == get<2>(synLocSampleVec[synIdx]).end())
            {
                synIdx++;  // ��������
                // ���Խ�磬�������겢��������
                if (synIdx >= synLocSampleVec.size())
                {
                    synIdx = 0;
                    synStart = 0;
                    synEnd = 0;
                    qrySynStart = 0;
                    qrySynEnd = 0;
                    whileBool = false;

                    return 0;
                }
            }
            
            synStart = get<0>(synLocSampleVec[synIdx]);
            synEnd = get<1>(synLocSampleVec[synIdx]);
            synIdx++;  // ��������
            qrySynStart = 0;
            qrySynEnd = 0;
            whileBool = true;
        }
        else
        {
            synIdx = 0;
            synStart = 0;
            synEnd = 0;
            qrySynStart = 0;
            qrySynEnd = 0;
            whileBool = false;
        }
    }
    
    return 0;
}


/**
 * @brief ��������sample���е�syn��������
 * 
 * @param inputFileName         syn_loc����ļ�
 * 
 * @return synLocSampleVecMap   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
*/
map<string, vector<tuple<int, int, vector<string> > > > COOR::build_syn_idx(const string & inputFileName)
{
    // ���湲��������
    map<string, vector<tuple<int, int, vector<string> > > > synLocSampleVecMap;  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >

    // open syntenic coordinations file
    GzChunkReader GzChunkReaderClass(inputFileName);

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

        synLocSampleVecMap[lineVec[0]].push_back(make_tuple(stoi(lineVec[1]), stoi(lineVec[2]), split(lineVec[4], ",")));  // �洢��������ref�ϵ������Լ���λ�����Ʒ����
    }

    return synLocSampleVecMap;
}


/**
 * @brief �õ�syn��qry�ϵ�����
 * 
 * @param sampleName                 ��Ʒ��
 * @param inputFileName              show-aligns ����ļ�
 * @param synLocSampleVecMap         map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param findRevBool                �Ƿ��������ȶ�����
 * 
 * @return chrStartSynQryLocMap      synAllStructure
*/
COOR::synAllStructure COOR::get_syn_coor(
    string sampleName, 
    const string& inputFileName, 
    const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap, 
    const bool& findRevBool
)
{
    // ��ȡ��������qry�ϵ�����
    // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
    synAllStructure chrStartSynQryLocMap;

    // ���û���ύ��Ʒ���������ļ���ȡ
    if (sampleName.empty())
    {
        vector<string> inputFileNameVecTmp = split(inputFileName, "/");  // ·�����
        sampleName = inputFileNameVecTmp[inputFileNameVecTmp.size() - 1];
        std::regex reg1(".aligns");  // �滻
        sampleName = regex_replace(sampleName, reg1, "");  // 'An-1.aligns' ɾ�� '.aligns'
    }

    // ��¼��Ʒ����
    chrStartSynQryLocMap.sampleName = sampleName;

    // �ȳ�ʼ��sylLocVecOutMap
    for (const auto& [chromosome, locationSampleVec] : synLocSampleVecMap)  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    {
        for (const auto& [synStart, synEnd, sampleVec] : locationSampleVec)  // vector<tuple<refStart, refEnd, vector<sample> > >
        {
            chrStartSynQryLocMap.chrStartSynQryLocMap[chromosome].emplace(synStart, make_tuple(0, 0));
        }
    }
    
    // ��ʱ��syn����
    string synChr = "";
    int synIdx = 0;  // ������ȡ�����Ե�����
    int synStart = 0;
    int synEnd = 0;
    // �洢qry��syn������
    int qrySynStart;
    int qrySynEnd;

    // ��ʱ��Ⱦɫ��ź�����
    string refChr = "";
    string qryChr = "";
    // alignmenmt����Ϣ
    string AliRefStrand = "";
    string AliQryStrand = "";
    int aliRefStart = 0;
    int aliRefEnd = 0;
    int aliQryStart = 0;
    int aliQryEnd = 0;
    // �����ж��Ƿ�ѭ������alignment�Ĳ���ֵ  (false)
    bool aliBool = false;
    // ���ڼ�¼syn�������ڸ�alignment�ĵڼ���
    int aliRowStartNum = 0;
    int aliRowEndNum = INT32_MAX;

    // �����е���Ϣ
    int refStart = 0;
    int refEnd = 0;
    string refSeq = "";
    int qryStart = 0;
    int qryEnd = 0;
    string qrySeq = "";

    // �������ʱboolֵ����ʹ�ã�ֻ����Ϊ�����ύ
    bool whileBoolTmp = true;

    // open aligns file
    GzChunkReader GzChunkReaderClass(inputFileName);

    // ��¼ѭ��������  ż��->ref  ����->qry
    int forIdx = -1;

    // read line
    string line;
    
    while (GzChunkReaderClass.read_line(line))
    {
        // ��������
        if (line.empty())
        {
            continue;
        }
        
        if (line == "\n" || 
            line[0] == ' ' || 
            line[0] == '=' || 
            line.find("--   END alignment ") != string::npos)  // �ո�/=/\n/END alignment��ֱ������
        {
            continue;
        }
        else  // ��ȡsown-aligns����
        {
            if (line.find("-- Alignments between") != string::npos)  // �µıȶ�Ⱦɫ��    -- Alignments between chr1 and chr1
            {
                // ����ַ���
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                refChr = lineVec[3];  // refȾɫ�����ʱֵ
                qryChr = lineVec[5];  // qryȾɫ�����ʱֵ
                
                // ��������Ϣ����
                synChr = refChr;
                synIdx = 0;
                synStart = 0;
                synEnd = 0;
                renew_syn_loc(
                    sampleName, 
                    synLocSampleVecMap, 
                    synChr, 
                    synIdx, 
                    synStart, 
                    synEnd, 
                    qrySynStart, 
                    qrySynEnd, 
                    whileBoolTmp
                );

                // ��������
                AliRefStrand = "";
                AliQryStrand = "";

                // ����ַ���
                line.clear();
                string().swap(line);
                continue;
            }
            else if (line.find("-- BEGIN alignment ") != string::npos)  // �µıȶ�����    -- BEGIN alignment [ +1 1078 - 68996 | +1 1 - 67961 ]
            {
                line = strip(line.erase(0, 21), '\n');  // ɾ�� '-- BEGIN alignment [ ' ��21���ַ�
                line = line.erase(line.size()-2, 2);  // ɾ�� ' ]' ��2���ַ���λ�����������

                vector<string> lineVec = split(line, " | ");  // �ָ� '+1 278 - 1703 -1 2148 - 751'

                // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ  ref
                string AliRefStrandTmp;
                int aliRefStartTmp;
                int aliRefEndTmp;
                tie(AliRefStrandTmp, aliRefStartTmp, aliRefEndTmp) = get_alignment_loc(
                    lineVec[0]
                );

                // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ  qry
                string AliQryStrandTmp;
                int aliQryStartTmp;
                int aliQryEndTmp;
                tie(AliQryStrandTmp, aliQryStartTmp, aliQryEndTmp) = get_alignment_loc(
                    lineVec[1]
                );

                // �������е�����
                forIdx = -1;

                // �����ڹ�����֮ǰ����������read��
                if (aliRefEndTmp < synStart)
                {
                    aliBool = false;
                }
                // �ж��Ƿ��������򻥲����У���������read��
                else if (!findRevBool && (AliRefStrandTmp == "-1" || AliQryStrandTmp == "-1"))
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Reverse alignment skipped -> " << synChr << " " << line << endl;
                    aliBool = false;
                }
                // ref��qry��Ⱦɫ�岻һ������������read��
                else if (refChr != qryChr)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Skipped alignment due to chromosomal mismatch -> " << refChr << " != " << qryChr << endl;
                    aliBool = false;
                }
                else
                {
                    // ��������
                    AliRefStrand = AliRefStrandTmp;
                    aliRefStart = aliRefStartTmp;
                    aliRefEnd = aliRefEndTmp;
                    AliQryStrand = AliQryStrandTmp;
                    aliQryStart = aliQryStartTmp;
                    aliQryEnd = aliQryEndTmp;
                    aliBool = true;

                    // ����syn���ڵĴ����У�
                    aliRowStartNum = ((synStart - aliRefStart)/49)*2-10;
                    aliRowEndNum = ((synEnd - aliRefStart)/49)*2-10;
                }
            }
            // ֻѭ�������ֿ�ͷ�Ͱ���syn����
            // ���map�������ˣ�ͬ�������������������  'synEnd==0'
            else if (isdigit(line[0]) != 0 && aliBool && synEnd > 0)   
            {
                forIdx++;  // ��¼����  ż��ref ����qry

                // ���С����ֵ��������
                if (forIdx < aliRowStartNum)
                {
                    continue;
                }

                // ����ַ���
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                if (( forIdx & 1 ) == 0)  // ż��
                {
                    // �������� (this)
                    refSeq = lineVec[1];
                    // ��ʱ��ref����
                    string refSeqTmp = refSeq;
                    // ɾ��  '.'  ���㳤��
                    refSeqTmp.erase(remove(refSeqTmp.begin(), refSeqTmp.end(), '.'), refSeqTmp.end());
                    if (AliRefStrand == "+1")  // ����ȶ�
                    {
                        refStart = stoi(lineVec[0]);
                        refEnd = refStart + refSeqTmp.size() - 1;
                    }
                    else  // ����ȶ�
                    {
                        refEnd = stoi(lineVec[0]);
                        refStart = refEnd - refSeqTmp.size() + 1;
                    }
                }
                else
                {
                    // �������� (this)
                    qrySeq = lineVec[1];
                    // ��ʱ��ref����
                    string qrySeqTmp = qrySeq;
                    // ɾ��  '.'  ���㳤��
                    qrySeqTmp.erase(remove(qrySeqTmp.begin(), qrySeqTmp.end(), '.'), qrySeqTmp.end());
                    if (AliQryStrand == "+1")  // ����ȶ�
                    {
                        qryStart = stoi(lineVec[0]);
                        qryEnd = qryStart + qrySeqTmp.size() - 1;
                    }
                    else  // ����ȶ�
                    {
                        qryEnd = stoi(lineVec[0]);
                        qryStart = qryEnd - qrySeqTmp.size() + 1;
                    }
                    
                    // ��������������ڸ�����֮ǰ�ˣ����¹����Ե����꣬����Ҫ��ֹ��ѭ�������Ե�'synEnd==0'ʱ����������ֵ��ˣ�����
                    while (synEnd > 0 && synEnd < refStart)
                    {
                        renew_syn_loc(
                            sampleName, 
                            synLocSampleVecMap, 
                            synChr, 
                            synIdx, 
                            synStart, 
                            synEnd, 
                            qrySynStart, 
                            qrySynEnd, 
                            whileBoolTmp
                        );

                        // ����syn���ڵĴ����У�
                        aliRowStartNum = ((synStart - aliRefStart)/49)*2-10;
                        aliRowEndNum = ((synEnd - aliRefStart)/49)*2-10;

                        // �жϸ����ȶ��Ƿ������syn���������±ߵ��ж�����
                        if (aliRefEnd < synStart)  // �����ڹ�����֮ǰ����������read
                        {
                            aliBool = false;
                        }
                    }

                    // �ж��Ƿ����whileѭ��
                    bool whileBool = true;

                    // �����꣬push����ͼ��
                    while (((refStart <= synStart && synStart <= refEnd) || 
                            (refStart <= synEnd && synEnd <= refEnd)) &&
                            whileBool)
                    {
                        int qrySynStartTmp = find_qry_syn(
                            synStart, 
                            refStart, 
                            refEnd, 
                            refSeq, 
                            AliQryStrand, 
                            qryStart, 
                            qryEnd, 
                            qrySeq
                        );

                        int qrySynEndTmp = find_qry_syn(
                            synEnd, 
                            refStart, 
                            refEnd, 
                            refSeq, 
                            AliQryStrand, 
                            qryStart, 
                            qryEnd, 
                            qrySeq
                        );

                        // ���Ƿ��ҵ����ҵ������ܹ�ϣ�������
                        syn_all_loc_push(
                            chrStartSynQryLocMap, 
                            refChr, 
                            qrySynStartTmp, 
                            qrySynEndTmp, 
                            refStart, 
                            refEnd, 
                            qryStart, 
                            qryEnd, 
                            synStart, 
                            synEnd, 
                            qrySynStart, 
                            qrySynEnd, 
                            whileBool, 
                            aliRowStartNum, 
                            aliRowEndNum
                        );

                        // ����������궼�ҵ��������Ѿ�����synEnd�����У�����syn������
                        if ((qrySynStart > 0 && qrySynEnd > 0) || 
                            (refStart <= synEnd && synEnd <= refEnd) || 
                            synEnd < refStart)
                        {
                            renew_syn_loc(
                                sampleName, 
                                synLocSampleVecMap, 
                                synChr, 
                                synIdx, 
                                synStart, 
                                synEnd, 
                                qrySynStart, 
                                qrySynEnd, 
                                whileBool
                            );
                            // �жϸ����ȶ��Ƿ������syn���������±ߵ��ж�����
                            if (aliRefEnd < synStart)  // �����ڹ�����֮ǰ����������read
                            {
                                aliBool = false;
                            }
                            // ����syn���ڵĴ����У�
                            aliRowStartNum = ((synStart - aliRefStart)/49)*2-10;
                            aliRowEndNum = ((synEnd - aliRefStart)/49)*2-10;
                        }
                    }
                }
            }
        }
    }

    return chrStartSynQryLocMap;
}


/**
 * @brief ������
 * 
 * @param synLocSampleVecMap            build_syn_idx����������   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param sampleChrStartSynQryLocMap    get_syn_coor����ֵ   map<sampleName, chrStartSynQryLocMap>   map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
 * @param outputFileName                ����ļ���
 * 
 * @return chrStartSynQryLocMap      synAllStructure
*/
int COOR::save_result(
    const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap,
    const map<string, unordered_map<string, unordered_map<int, tuple<int, int> > > >& sampleChrStartSynQryLocMap,
    const string& outputFileName
)
{
    SAVE SAVEClass(outputFileName);

    stringstream outStream; // ʹ�� stringstream �����ַ���ƴ��
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // �����СΪ 10mb
    outStream.str().reserve(CACHE_SIZE);

    // ������
    for (const auto& [chromosome, locationSampleVec] : synLocSampleVecMap)  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    {
        // ��¼��һ�����������꣬��ֹ������֮���м��
        uint64_t preSynEnd = 0;

        for (const auto& [synStart, synEnd, sampleVec] : locationSampleVec)  // vector<tuple<refStart, refEnd, vector<sample> > >
        {
            // �����ж��Ƿ��м��
            if(preSynEnd != 0 && synStart - preSynEnd > 1)
            {
                outStream << chromosome << '\t' << preSynEnd + 1 << '\t' << synStart - 1;
                
                for (const auto& it3 : sampleChrStartSynQryLocMap)  // map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
                {
                    outStream << '\t' << it3.first << '\t' << 0 << '\t' << 0;
                }
                outStream << '\n';
            }

            preSynEnd = synEnd;  // ��¼��һ�������Ե�����

            // ��ǰ�ڵ����Ϣ
            outStream << chromosome << '\t' << synStart << '\t' << synEnd;
            for (const auto& it3 : sampleChrStartSynQryLocMap)  // map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
            {
                const auto& chrStartSynQryLoc = it3.second.at(chromosome).at(synStart);
                outStream << '\t' << it3.first << '\t' << get<0>(chrStartSynQryLoc) << '\t' << get<1>(chrStartSynQryLoc);
            }
            outStream << '\n';

            if (outStream.tellp() >= CACHE_SIZE)  // �����СΪ 10mb
            {
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // ��� stringstream
                outStream.str(string());
                outStream.clear();
            }
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

    return 0;
}