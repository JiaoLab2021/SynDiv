// g++ -c src/coor.cpp -O3 -std=c++17
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <future>
#include <malloc.h>
#include "../include/coor.hpp"
#include "../include/ThreadPool.hpp"

using namespace std;

// 全局变量
int thresholdLength = 1000;  // 共线性坐标ref和qry长度比阈值

// 是否调试代码
bool debugCoor = false;


int main_coor(int argc, char* argv[])
{
    // loc 输出文件
    string synLocFileName;

    // 输入文件列表和名称  show-aligns
    vector<string> inputFiles;
    vector<string> inputTitles;
    // 判断是否有 titlename
    bool haveTitles = false;

    // 输出文件名
    string outputFileName;

    // 线程数
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
    ThreadPool pool(threads);  // 进程池

    // 初始化线程池
    pool.init();

    // 保存多线程的结果
    vector<future<COOR::synAllStructure> > chrStartSynQryLocMapVec;

    // 获取共线性在qry上的坐标
    int indexTmp = 0;  // 记录样品名的索引
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
        // 多线程提交并保存结果
        chrStartSynQryLocMapVec.push_back(
            pool.submit(
                COOR::get_syn_coor, 
                sampleName, 
                it, 
                ref(synLocSampleVecMap), 
                false
            )
        );

        indexTmp++;  // 索引叠加
    }

    // 多线程结果保存
    map<string, unordered_map<string, unordered_map<int, tuple<int, int> > > > sampleChrStartSynQryLocMap;
    for (size_t i = 0; i < chrStartSynQryLocMapVec.size(); i++)
    {
        COOR::synAllStructure synAllStructureTmp = move(chrStartSynQryLocMapVec[i].get());
        sampleChrStartSynQryLocMap[move(synAllStructureTmp.sampleName)] = move(synAllStructureTmp.chrStartSynQryLocMap);
    }
    chrStartSynQryLocMapVec.clear();
    vector<future<COOR::synAllStructure> >().swap(chrStartSynQryLocMapVec);
    malloc_trim(0);	// 0 is for heap memory

    // 关闭线程池
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
 * @brief 解析 '-- BEGIN alignment ' 字段
 * 
 * @param informationTmp    '+1 278 - 1703'
 * 
 * @return tuple<strand, start, end>
*/
tuple<string, int, int> COOR::get_alignment_loc(string informationTmp)
{
    // 获取比对方向，起始和终止的信息
    istringstream iss(informationTmp);  // +1 278 - 1703
    string strandTmp;
    int startTmp, endTmp;
    char sep;

    iss >> strandTmp >> startTmp >> sep >> endTmp;

    // 如果为反向比对，交换起始和终止坐标
    if (strandTmp == "-1") {
        swap(startTmp, endTmp);
    }

    return make_tuple(strandTmp, startTmp, endTmp);
}


/**
 * @brief 找qry上共线性的坐标
 * 
 * @param synLoc        共线性的坐标
 * @param refStart      ref的起始
 * @param refEnd        ref的终止
 * @param refSeq        ref的序列
 * @param qryStrand     qry比对方向
 * @param qryStart      qry的起始
 * @param qryEnd        qry的终止
 * @param qrySeq        qry的序列
 * 
 * @return qrySynLoc    qry上syn的坐标， 0-没找到
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

    if (refStart <= synLoc && synLoc <= refEnd)  // 如果包含了共线性坐标
    {
        int iIdx = 0;  // 记录synStart对应的位置

        --refStart;  // 坐标先减1
        for (size_t i = 0; i < refSeq.size(); ++i)  // 循环ref序列
        {
            if (refSeq[i] != '.')
            {
                ++refStart;  // ref坐标递增
            }

            if (refStart == synLoc)  // 找syn坐标的索引
            {
                iIdx = i;
                break;  // 找到坐标后退出该循环
            }
        }

        if (AliQryStrand == "+1")  // qry是正向比对时
        {
            --qryStart;  // 坐标先减1
            for (size_t i = 0; i < qrySeq.size(); ++i)
            {
                if (qrySeq[i] != '.')
                {
                    ++qryStart;  // ref坐标递增
                }
                
                if (i == iIdx)  // 循环到该位置后退出
                {
                    break;  // 找到坐标后退出该循环
                }
            }
            qrySynLoc = qryStart;
        }
        else  // qry是反向比对时
        {
            ++qryEnd;  // 坐标先加1
            for (size_t i = 0; i < qrySeq.size(); ++i)
            {
                if (qrySeq[i] != '.')
                {
                    --qryEnd;  // ref坐标递减
                }
                
                if (i == iIdx)  // 循环到该位置后退出
                {
                    break;  // 找到坐标后退出该循环
                }
            }
            qrySynLoc = qryEnd;
        }
    }

    return qrySynLoc;
}


/**
 * @brief 找qry上共线性的坐标
 * 
 * @param chrStartSynQryLocMap  synAllStructure
 * @param refChr                染色体号
 * @param qrySynStartTmp        find_qry_syn找的qrySynStart
 * @param qrySynEndTmp          find_qry_syn找的qrySynEnd
 * @param refStart              该比对行ref的起始
 * @param refEnd                该比对行ref的终止
 * @param qryStart              该比对行qry的起始
 * @param qryEnd                带比对行qry的终止
 * @param synStart              目前要找的syn的ref起始
 * @param synEnd                目前要找的syn的ref终止
 * @param qrySynStart           syn的qry起始
 * @param qrySynEnd             syn的qry终止
 * @param whileBool             判断是否还需要while循环
 * @param aliRowStartNum        syn起始位置在该alignment的大致行数
 * @param aliRowEndNum          syn终止位置在该alignment的大致行数
 * 
 * @return qrySynLoc    qry上syn的坐标， 0-没找到
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
        if (qrySynStart > 0) {  // 之前找到过qrySynStart，选择离preQrySynStart最近的序列
            if (abs(qrySynStart - synStart) > abs(qrySynStartTmp - synStart))
            {
                qrySynStart = qrySynStartTmp;
            } else if (debugCoor) {
                cerr << "skip_start:" << qrySynStartTmp << endl;
            }
        }
        else  // 没有找到的话直接赋值
        {
            qrySynStart = qrySynStartTmp;
        }
        
        aliRowStartNum = aliRowEndNum;  // 重置aliRowStartNum，跳至aliRowEndNum找终止坐标
    }
    if (qrySynEndTmp != 0 && qrySynStart != 0)  // 起始必须先找到
    {
        if (qrySynEnd > 0)  // 之前找到过qrySynEnd，选择离preQrySynEnd最近的序列
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
        else  // 没有找到的话直接赋值
        {
            qrySynEnd = qrySynEndTmp;
        }

        // 大于0了再赋值，比对长度相差倍数小于thresholdLength了再赋值
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
 * @brief 更新临时的syn坐标
 * 
 * @param sampleName            样品名，用于判断该行是否要进行判断
 * @param synLocSampleVecMap    map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param synChr                染色体
 * @param synIdx                syn在该染色体的vector中索引
 * @param synStart              syn的起始
 * @param synEnd                syn的终止
 * @param qrySynStart           qry上syn的起始
 * @param qrySynEnd             qry上syn的终止
 * @param whileBool             判断是否跳出while循环
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
    map<string, vector<tuple<int, int, vector<string> > > >::const_iterator findIter = synLocSampleVecMap.find(synChr);  //  在共线性map中找染色体
    if (findIter != synLocSampleVecMap.end())  // 找到了
    {
        auto& synLocSampleVec = findIter->second;  // vector<tuple<int, int, vector<string> > >

        if (synIdx < synLocSampleVec.size())  // 防止越界
        {
            // 如果该共线性没有对应的样品名，继续下一个坐标
            while (find(get<2>(synLocSampleVec[synIdx]).begin(), get<2>(synLocSampleVec[synIdx]).end(), sampleName) == get<2>(synLocSampleVec[synIdx]).end())
            {
                synIdx++;  // 更新坐标
                // 如果越界，归零坐标并跳出函数
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
            synIdx++;  // 更新坐标
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
 * @brief 构件所有sample共有的syn坐标索引
 * 
 * @param inputFileName         syn_loc输出文件
 * 
 * @return synLocSampleVecMap   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
*/
map<string, vector<tuple<int, int, vector<string> > > > COOR::build_syn_idx(const string & inputFileName)
{
    // 保存共线性坐标
    map<string, vector<tuple<int, int, vector<string> > > > synLocSampleVecMap;  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >

    // open syntenic coordinations file
    GzChunkReader GzChunkReaderClass(inputFileName);

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

        synLocSampleVecMap[lineVec[0]].push_back(make_tuple(stoi(lineVec[1]), stoi(lineVec[2]), split(lineVec[4], ",")));  // 存储共线性在ref上的坐标以及该位点的样品名称
    }

    return synLocSampleVecMap;
}


/**
 * @brief 得到syn在qry上的坐标
 * 
 * @param sampleName                 样品名
 * @param inputFileName              show-aligns 输出文件
 * @param synLocSampleVecMap         map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param findRevBool                是否遍历反向比对序列
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
    // 获取共线性在qry上的坐标
    // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
    synAllStructure chrStartSynQryLocMap;

    // 如果没有提交样品名，根据文件提取
    if (sampleName.empty())
    {
        vector<string> inputFileNameVecTmp = split(inputFileName, "/");  // 路径拆分
        sampleName = inputFileNameVecTmp[inputFileNameVecTmp.size() - 1];
        std::regex reg1(".aligns");  // 替换
        sampleName = regex_replace(sampleName, reg1, "");  // 'An-1.aligns' 删除 '.aligns'
    }

    // 记录样品名称
    chrStartSynQryLocMap.sampleName = sampleName;

    // 先初始化sylLocVecOutMap
    for (const auto& [chromosome, locationSampleVec] : synLocSampleVecMap)  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    {
        for (const auto& [synStart, synEnd, sampleVec] : locationSampleVec)  // vector<tuple<refStart, refEnd, vector<sample> > >
        {
            chrStartSynQryLocMap.chrStartSynQryLocMap[chromosome].emplace(synStart, make_tuple(0, 0));
        }
    }
    
    // 临时的syn坐标
    string synChr = "";
    int synIdx = 0;  // 用于提取共线性的坐标
    int synStart = 0;
    int synEnd = 0;
    // 存储qry上syn的坐标
    int qrySynStart;
    int qrySynEnd;

    // 临时的染色体号和坐标
    string refChr = "";
    string qryChr = "";
    // alignmenmt的信息
    string AliRefStrand = "";
    string AliQryStrand = "";
    int aliRefStart = 0;
    int aliRefEnd = 0;
    int aliQryStart = 0;
    int aliQryEnd = 0;
    // 用于判断是否循环该条alignment的布尔值  (false)
    bool aliBool = false;
    // 用于记录syn的坐标在该alignment的第几行
    int aliRowStartNum = 0;
    int aliRowEndNum = INT32_MAX;

    // 具体行的信息
    int refStart = 0;
    int refEnd = 0;
    string refSeq = "";
    int qryStart = 0;
    int qryEnd = 0;
    string qrySeq = "";

    // 构造的临时bool值，不使用，只是作为参数提交
    bool whileBoolTmp = true;

    // open aligns file
    GzChunkReader GzChunkReaderClass(inputFileName);

    // 记录循环的索引  偶数->ref  奇数->qry
    int forIdx = -1;

    // read line
    string line;
    
    while (GzChunkReaderClass.read_line(line))
    {
        // 跳过空行
        if (line.empty())
        {
            continue;
        }
        
        if (line == "\n" || 
            line[0] == ' ' || 
            line[0] == '=' || 
            line.find("--   END alignment ") != string::npos)  // 空格/=/\n/END alignment，直接跳过
        {
            continue;
        }
        else  // 获取sown-aligns坐标
        {
            if (line.find("-- Alignments between") != string::npos)  // 新的比对染色体    -- Alignments between chr1 and chr1
            {
                // 拆分字符串
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                refChr = lineVec[3];  // ref染色体的临时值
                qryChr = lineVec[5];  // qry染色体的临时值
                
                // 共线性信息更新
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

                // 重置坐标
                AliRefStrand = "";
                AliQryStrand = "";

                // 清空字符串
                line.clear();
                string().swap(line);
                continue;
            }
            else if (line.find("-- BEGIN alignment ") != string::npos)  // 新的比对坐标    -- BEGIN alignment [ +1 1078 - 68996 | +1 1 - 67961 ]
            {
                line = strip(line.erase(0, 21), '\n');  // 删除 '-- BEGIN alignment [ ' 共21个字符
                line = line.erase(line.size()-2, 2);  // 删除 ' ]' 共2个字符，位置在最后两个

                vector<string> lineVec = split(line, " | ");  // 分割 '+1 278 - 1703 -1 2148 - 751'

                // 获取比对方向，起始和终止的信息  ref
                string AliRefStrandTmp;
                int aliRefStartTmp;
                int aliRefEndTmp;
                tie(AliRefStrandTmp, aliRefStartTmp, aliRefEndTmp) = get_alignment_loc(
                    lineVec[0]
                );

                // 获取比对方向，起始和终止的信息  qry
                string AliQryStrandTmp;
                int aliQryStartTmp;
                int aliQryEndTmp;
                tie(AliQryStrandTmp, aliQryStartTmp, aliQryEndTmp) = get_alignment_loc(
                    lineVec[1]
                );

                // 重置序列的索引
                forIdx = -1;

                // 坐标在共线性之前，跳过该条read。
                if (aliRefEndTmp < synStart)
                {
                    aliBool = false;
                }
                // 判断是否跳过反向互补序列，跳过该条read。
                else if (!findRevBool && (AliRefStrandTmp == "-1" || AliQryStrandTmp == "-1"))
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Reverse alignment skipped -> " << synChr << " " << line << endl;
                    aliBool = false;
                }
                // ref和qry的染色体不一样，跳过该条read。
                else if (refChr != qryChr)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Skipped alignment due to chromosomal mismatch -> " << refChr << " != " << qryChr << endl;
                    aliBool = false;
                }
                else
                {
                    // 重置坐标
                    AliRefStrand = AliRefStrandTmp;
                    aliRefStart = aliRefStartTmp;
                    aliRefEnd = aliRefEndTmp;
                    AliQryStrand = AliQryStrandTmp;
                    aliQryStart = aliQryStartTmp;
                    aliQryEnd = aliQryEndTmp;
                    aliBool = true;

                    // 计算syn所在的大致行，
                    aliRowStartNum = ((synStart - aliRefStart)/49)*2-10;
                    aliRowEndNum = ((synEnd - aliRefStart)/49)*2-10;
                }
            }
            // 只循环是数字开头和包含syn的行
            // 如果map遍历完了，同样跳过，遍历完代表是  'synEnd==0'
            else if (isdigit(line[0]) != 0 && aliBool && synEnd > 0)   
            {
                forIdx++;  // 记录行数  偶数ref 奇数qry

                // 如果小于阈值的行跳过
                if (forIdx < aliRowStartNum)
                {
                    continue;
                }

                // 拆分字符串
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                if (( forIdx & 1 ) == 0)  // 偶数
                {
                    // 重置坐标 (this)
                    refSeq = lineVec[1];
                    // 临时的ref序列
                    string refSeqTmp = refSeq;
                    // 删除  '.'  计算长度
                    refSeqTmp.erase(remove(refSeqTmp.begin(), refSeqTmp.end(), '.'), refSeqTmp.end());
                    if (AliRefStrand == "+1")  // 正向比对
                    {
                        refStart = stoi(lineVec[0]);
                        refEnd = refStart + refSeqTmp.size() - 1;
                    }
                    else  // 反向比对
                    {
                        refEnd = stoi(lineVec[0]);
                        refStart = refEnd - refSeqTmp.size() + 1;
                    }
                }
                else
                {
                    // 重置坐标 (this)
                    qrySeq = lineVec[1];
                    // 临时的ref序列
                    string qrySeqTmp = qrySeq;
                    // 删除  '.'  计算长度
                    qrySeqTmp.erase(remove(qrySeqTmp.begin(), qrySeqTmp.end(), '.'), qrySeqTmp.end());
                    if (AliQryStrand == "+1")  // 正向比对
                    {
                        qryStart = stoi(lineVec[0]);
                        qryEnd = qryStart + qrySeqTmp.size() - 1;
                    }
                    else  // 反向比对
                    {
                        qryEnd = stoi(lineVec[0]);
                        qryStart = qryEnd - qrySeqTmp.size() + 1;
                    }
                    
                    // 如果共线性区间在该区间之前了，更新共线性的坐标，但是要防止死循环，所以当'synEnd==0'时代表遍历完字典了，跳过
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

                        // 计算syn所在的大致行，
                        aliRowStartNum = ((synStart - aliRefStart)/49)*2-10;
                        aliRowEndNum = ((synEnd - aliRefStart)/49)*2-10;

                        // 判断该条比对是否包含了syn，不包含下边的行都跳过
                        if (aliRefEnd < synStart)  // 坐标在共线性之前，跳过该条read
                        {
                            aliBool = false;
                        }
                    }

                    // 判断是否进入while循环
                    bool whileBool = true;

                    // 找坐标，push到总图中
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

                        // 看是否找到，找到了往总哈希表中添加
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

                        // 如果两个坐标都找到，或者已经到达synEnd所在行，更新syn的坐标
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
                            // 判断该条比对是否包含了syn，不包含下边的行都跳过
                            if (aliRefEnd < synStart)  // 坐标在共线性之前，跳过该条read
                            {
                                aliBool = false;
                            }
                            // 计算syn所在的大致行，
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
 * @brief 保存结果
 * 
 * @param synLocSampleVecMap            build_syn_idx构建的索引   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param sampleChrStartSynQryLocMap    get_syn_coor返回值   map<sampleName, chrStartSynQryLocMap>   map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
 * @param outputFileName                输出文件名
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

    stringstream outStream; // 使用 stringstream 代替字符串拼接
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // 缓存大小为 10mb
    outStream.str().reserve(CACHE_SIZE);

    // 保存结果
    for (const auto& [chromosome, locationSampleVec] : synLocSampleVecMap)  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    {
        // 记录上一个共线性坐标，防止共线性之间有间隔
        uint64_t preSynEnd = 0;

        for (const auto& [synStart, synEnd, sampleVec] : locationSampleVec)  // vector<tuple<refStart, refEnd, vector<sample> > >
        {
            // 否则判断是否有间隔
            if(preSynEnd != 0 && synStart - preSynEnd > 1)
            {
                outStream << chromosome << '\t' << preSynEnd + 1 << '\t' << synStart - 1;
                
                for (const auto& it3 : sampleChrStartSynQryLocMap)  // map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
                {
                    outStream << '\t' << it3.first << '\t' << 0 << '\t' << 0;
                }
                outStream << '\n';
            }

            preSynEnd = synEnd;  // 记录上一个共线性的坐标

            // 当前节点的信息
            outStream << chromosome << '\t' << synStart << '\t' << synEnd;
            for (const auto& it3 : sampleChrStartSynQryLocMap)  // map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
            {
                const auto& chrStartSynQryLoc = it3.second.at(chromosome).at(synStart);
                outStream << '\t' << it3.first << '\t' << get<0>(chrStartSynQryLoc) << '\t' << get<1>(chrStartSynQryLoc);
            }
            outStream << '\n';

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

    if (outStream.tellp() > 0)  // 最后写一次
    {
        string outTxt = outStream.str();
        SAVEClass.save(outTxt);
        // 清空 stringstream
        outStream.str(string());
        outStream.clear();
    }

    return 0;
}