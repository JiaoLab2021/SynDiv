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

    // 如果调试代码，线程数设置为1
    if (debugCal)
    {
        threadsCal = 1;
    }

    /* ************************************ Parse Parameters ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Parse Parameters ..." << endl;
    // 获取 aligns/syri.out 文件字典
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

    // 所有组合数量
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
 * @brief 解析参数
 * 
 * @param alignsTitles     aligns 文件 title
 * @param alignsVec        aligns 路径Vec
 * @param syriConfigFileName   syri 输出的配置文件
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

    uint32_t indexTmp = 0;  // 记录 title 的索引
    for (auto iter1 : alignsVec)
    {
        string alignsPathTmp = iter1;
        string ailgnsTitleTmp;

        // 如果含有 title
        if (alignsTitles.size() > 0)
        {
            ailgnsTitleTmp = alignsTitles[indexTmp];
        }
        else
        {
            vector<string> alignsPathVecTmp = split(alignsPathTmp, "/");  // 路径拆分
            ailgnsTitleTmp = alignsPathVecTmp[alignsPathVecTmp.size() - 1];  // 最后一个
            std::regex reg1(".aligns");  // 替换
            ailgnsTitleTmp = regex_replace(ailgnsTitleTmp, reg1, "");  // 'An-1.aligns' 删除 '.aligns'
        }

        // 赋值
        alignsMap[ailgnsTitleTmp] = alignsPathTmp;

        indexTmp++;  // 索引迭代
    }

    /* ************************************* syri.out ************************************* */
    map<string, map<string, string>> sampleSampleSyriMap;  // map<sample1, map<sample2, syriOutPath> >

    // open file
    GzChunkReader GzChunkReaderClass(syriConfigFileName);

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

        // 检查列出是否为3
        if (lineVec.size() != 3)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: The '" << syriConfigFileName << "' dataframe does not contain three columns." << endl;
            exit(1);
        }

        // 样品和路径信息
        const string& sample1 = lineVec[0];
        const string& sample2 = lineVec[1];
        const string& syriPath = lineVec[2];

        // 初始化
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
 * @brief 解析 '-- BEGIN alignment ' 字段
 * 
 * @param infoTmp   删除 '-- BEGIN alignment ' 后按 ' | ' 分割的字段，'+1 278 - 1703'
 * 
 * @return tuple<strand, start, end>
*/
tuple<string, uint32_t, uint32_t> CALNAME::COORTRANS::get_alignment_loc(
    string infoTmp
)
{
    // 获取比对方向，起始和终止的信息  ref
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
 * @brief 下一个比对段的坐标
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
            // 跳过空行
            if (line.empty())
            {
                return 0;
            }

            // 新的染色体对
            if (line.find("-- Alignments between") != string::npos)  // 新的比对染色体    -- Alignments between chr1 and chr1
            {
                // 拆分字符串
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                string refChr = lineVec[3];  // ref染色体的临时值
                string qryChr = lineVec[5];  // qry染色体的临时值

                if (refChr != qryChr)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Skipped alignment due to chromosomal mismatch -> " << refChr << " != " << qryChr << endl;
                    chrBool = false;  // 记录染色体符合规定
                }
                else
                {
                    chrBool = true;  // 记录染色体不符合规定
                }

                // 该判断防止丢失第一个 "-- BEGIN alignment"
                if (chr.length() == 0)
                {
                    findChrBool = true;  // 记录是需要的染色体
                }
                
                // 比对染色体信息更新
                chr = refChr;
            }
            // 染色体必须符合规定
            else if (line.find("-- BEGIN alignment ") != string::npos && chrBool && findChrBool)  // 新的比对坐标   -- BEGIN alignment [ +1 1078 - 68996 | +1 1 - 67961 ]
            {
                line = strip(line.erase(0, 21), '\n');  // 删除 '-- BEGIN alignment [ ' 共21个字符
                line = line.erase(line.size()-2, 2);  // 删除 ' ]' 共2个字符，位置在最后两个

                vector<string> lineVec = split(line, " | ");  // 分割 '+1 278 - 1703 -1 2148 - 751'

                // 获取比对方向，起始和终止的信息  ref
                tie(aliRefStrand, aliRefStart, aliRefEnd) = get_alignment_loc(
                    lineVec[0]
                );

                // 获取比对方向，起始和终止的信息  qry
                tie(aliQryStrand, aliQryStart, aliQryEnd) = get_alignment_loc(
                    lineVec[1]
                );
            }
            // 只循环是数字开头和包含syn的行
            // 染色体必须符合规定，比对起始和终止必须符合规定
            else if (isdigit(line[0]) != 0 && chrBool && findChrBool && aliBool)
            {
                /* ***************************** ref ***************************** */
                // 拆分字符串
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // 重置坐标 (this)
                refSeq = lineVec[1];
                // 临时的ref序列
                string refSeqTmp = refSeq;
                // 删除  '.'  计算长度
                refSeqTmp.erase(remove(refSeqTmp.begin(), refSeqTmp.end(), '.'), refSeqTmp.end());
                if (aliRefStrand == "+")  // 正向比对
                {
                    refStart = stoi(lineVec[0]);
                    refEnd = refStart + refSeqTmp.size() - 1;
                    refLoci = refStart - 1;  // 刷新坐标
                }
                else  // 反向比对
                {
                    refEnd = stoi(lineVec[0]);
                    refStart = refEnd - refSeqTmp.size() + 1;
                    refLoci = refEnd + 1;  // 刷新坐标
                }

                /* ***************************** qry ***************************** */
                // 读取下一行
                GzChunkReaderClass_->read_line(line);

                // 拆分字符串
                std::istringstream iss1(line);
                vector<string> infoVecTmp(std::istream_iterator<std::string>{iss1}, std::istream_iterator<std::string>());
                lineVec.assign(infoVecTmp.begin(), infoVecTmp.end());

                // 重置坐标 (this)
                qrySeq = lineVec[1];
                // 临时的ref序列
                string qrySeqTmp = qrySeq;
                // 删除  '.'  计算长度
                qrySeqTmp.erase(remove(qrySeqTmp.begin(), qrySeqTmp.end(), '.'), qrySeqTmp.end());
                if (aliQryStrand == "+")  // 正向比对
                {
                    qryStart = stoi(lineVec[0]);
                    qryEnd = qryStart + qrySeqTmp.size() - 1;
                    qryLoci = qryStart - 1;  // 刷新坐标
                }
                else  // 反向比对
                {
                    qryEnd = stoi(lineVec[0]);
                    qryStart = qryEnd - qrySeqTmp.size() + 1;
                    qryLoci = qryEnd + 1;  // 刷新坐标
                }

                /* ***************************** 刷新坐标 ***************************** */
                idxTmp = -1;
            }
        }
        else  // 文件遍历完
        {
            // 全部坐标归零
            // 遍历文件时更改的值
            // 染色体号
            chr.clear();
            // alignment的方向、起始和终止
            aliRefStrand.clear();
            aliRefStart = 0;
            aliRefEnd = 0;
            aliQryStrand.clear();
            aliQryStart = 0;
            aliQryEnd = 0;
            // 文件流指针所在行的比对信息
            refStart = 0;
            refSeq.clear();
            refEnd = 0;
            qryStart = 0;
            qrySeq.clear();
            qryEnd = 0;
            // 记录遍历到了 '_ref_seq' 以及对应的坐标
            refLoci = 0;
            idxTmp = -1;
            qryLoci = 0;

            endBool = true;  // 记录文件是否遍历完
        }
    }

    return 0;
}

/**
 * @brief 寻找loci
 * 
 * @param chr_  要找的染色体
 * @param refLoci_  要找的坐标
 * @return uint32_t  0->遍历完/没找到, >0->坐标
**/
uint32_t CALNAME::COORTRANS::find_loci(
    string chr_, 
    uint32_t refLoci_
)
{
    /* ***************************** 遍历完，返回 ‘0’ ***************************** */
    if (endBool)
    {
        return 0;
    }
    
    /* ***************************** 没有遍历完 ***************************** */
    // 判断染色体是不是想要的
    while (chr != chr_ && !endBool)
    {
        findChrBool = false;  // 记录不是需要的染色体
        next_loci();  // next line
    }
    findChrBool = true;  // 记录是需要的染色体
    // 如果遍历完，返回 '0'
    if (endBool)
    {
        return 0;
    }

    // 判断alignment终止是不是小于  'refLoci_'
    while (aliRefEnd < refLoci_ && !endBool)  // 上边已经判断了染色体是否相符，所以此处不用判断
    {
        aliBool = false;  // 记录不是要的 align
        next_loci();  // next line
    }
    aliBool = true;  // 记录是要的 align
    // 如果遍历完，返回 '0'
    if (endBool)
    {
        return 0;
    }

    // 文件指针指向想要的行
    while (refEnd < refLoci_ && !endBool)  // 上边已经判断了染色体是否相符，所以此处不用判断
    {
        next_loci();  // next line
    }
    // 如果遍历完，返回 '0'
    if (endBool)
    {
        return 0;
    }

    // 找到想要的行
    if (refStart <= refLoci_ && refLoci_ <= refEnd and refLoci < refLoci_)
    {
        for (size_t idx = idxTmp + 1; idx < refSeq.length(); idx++)
        {
            // ref
            if (refSeq[idx] != '.')
            {
                if (aliRefStrand == "+")  // 正向
                {
                    refLoci++;
                }
                else  // 反向
                {
                    refLoci--;
                }
            }

            // qry
            if (qrySeq[idx] != '.')
            {
                if (aliQryStrand == "+")  // 正向
                {
                    qryLoci++;
                }
                else  // 反向
                {
                    qryLoci--;
                }
            }

            idxTmp = idx;  // 刷新索引

            // 判断是不是要找的坐标
            if (refLoci == refLoci_)
            {
                break;  // 跳出循环
            }
        }
    }
    else if (refStart <= refLoci_ && refLoci_ <= refEnd and refLoci == refLoci_)
    {
        return qryLoci;  // 返回坐标
    }
    else  // 没找到
    {
        return 0;
    }

    return qryLoci;
}


/**
 * @brief 下一个共线性坐标
 * 
 * @return 0
**/
int CALNAME::SYRIOUT::next_loci()
{
    // 没有遍历完
    if (!endBool)
    {
        // read line
        string line;
        if (GzChunkReaderClass_->read_line(line))
        {
            // 跳过空行或者不含 SYNAL 的行
            if (line.empty() || line.find("SYNAL") == string::npos)
            {
                return 0;
            }

            // 拆分
            std::istringstream iss(line);
            vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

            // 检查列数是否正常
            if (lineVec.size() != 12)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "'" << fileName_ << "': does not contain 12 columns. -> " << line << endl;
                exit(1);
            }

            // 赋值
            chr = lineVec[0];
            refStart = stoi(lineVec[1]);
            refEnd = stoi(lineVec[2]);
            qryStart = stoi(lineVec[6]);
            qryEnd = stoi(lineVec[7]);
        }
        else  // 文件遍历完
        {
            // 全部变量清空
            chr.clear();
            refStart = 0;
            refEnd = 0;
            qryStart = 0;
            qryEnd = 0;

            endBool = true;  // 记录文件遍历完
        }
    }
    return 0;
}

/**
 * @brief 寻找loci
 * 
 * @param chr_  要找的染色体
 * @param refLoci_  要找的坐标
 * 
 * @return int  -1->遍历完，0->没找到, 1->坐标
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
 * @brief 计算
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
        
        // 临时结构体
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

                    // 如果和上一个的 synNum 不一样则添加
                    // if (synNumTmp != preSynNum || (synNumTmp == 0 && preSynNum == 0))
                    if (synNumTmp != preSynNum)
                    {
                        sampleSynOutMap[refLoci] = synNumTmp;
                        preSynNum = synNumTmp;
                    }
                }
            }

            // 计算比对上区域的得分
            for (size_t refLoci = refStart; refLoci < refEnd + 1; refLoci++)
            {
                uint32_t synNumTmp = 0;  // 共线性数量，临时

                for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
                {
                    string sampleName2 = SynCoorTmp.sampleNameVec[i];

                    // 如果在 'coor.txt' 中有坐标，代表是共线性，直接记录
                    // SampleLociMap  ->  // map<sample, tuple<start, end> >
                    if (
                        get<0>(SampleLociMap.at(sampleName)) > 0 && 
                        get<0>(SampleLociMap.at(sampleName2)) > 0
                    )  // both syntenic
                    {
                        synNumTmp++;  // 记为 syn

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
                        // 如果是reference与大家没坐标，那就是都没坐标
                        if (sampleName == "reference")
                        {
                            continue;
                        }

                        // 如果是其它
                        uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName
                        uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName2

                        if (lociA == 0 && lociB == 0) // both del
                        {
                            synNumTmp++;  // 记为 syn
                            
                            // debug
                            if (debugCal)
                            {
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                            }
                        }
                        else if (lociA != 0 && lociB == 0)  // 一个有坐标，另一个没有，one del
                        {
                            continue;
                        }
                        else if (lociA == 0 && lociB != 0)  // 一个有坐标，另一个没有，one del
                        {
                            continue;
                        }
                        else  // 都有坐标，再判断坐标对应的A和B是不是共线性
                        {
                            if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A和B是共线性
                            {
                                synNumTmp++;  // 记为 syn

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

                // 如果和上一个的 synNum 不一样则添加
                if (synNumTmp != preSynNum)
                {
                    sampleSynOutMap[refLoci] = synNumTmp;
                    preSynNum = synNumTmp;
                }
            }

            // 更新坐标
            refEndTmp = refEnd;
        }
    }

    // 判断是否到达染色体末尾
    for (const auto& [chrTmp, refPosTmp] : refChrMapTmp)  // map<chr, length>
    {
        // 初始化字典
        auto& sampleSynOutMap = CALSTRUCTURETMP.sampleSynOutMap[chrTmp];
        auto& sampleSynOutMapTmp =  CALSTRUCTURETMP.sampleSynOutMapTmp[chrTmp];
        
        // 检查参考基因组中是否有对应的染色体
        map<string, uint32_t>::const_iterator findIter1 = refLenMap.find(chrTmp);
        if (findIter1 == refLenMap.end())
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The reference genome does not include '" << chrTmp << "'." << endl;
            exit(1);
        }

        uint32_t refLen = findIter1->second;  // 染色体长度

        // 如果到达末尾，跳出 loop
        if (refPosTmp >= refLen) {

            break;
        }
        
        // 如果没到达染色体末尾
        uint32_t preSynNum = 0;
        for (uint32_t refLoci = refPosTmp + 1; refLoci < refLen + 1; refLoci++)
        {
            uint32_t synNumTmp = 0;  // 共线性数量，临时

            for (size_t i = idxTmp; i < SynCoorTmp.sampleNameVec.size(); i++)
            {
                string sampleName2 = SynCoorTmp.sampleNameVec[i];

                // 如果是reference与大家没坐标，那就是都没坐标
                if (sampleName == "reference")
                {
                    continue;
                }
                
                // 如果是其它
                uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName
                uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // 参考基因组的坐标转到sampleName2

                if (lociA == 0 && lociB == 0) // both del
                {
                    synNumTmp++;  // 记为 syn

                    // debug
                    if (debugCal)
                    {
                        sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                        cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                    }
                }
                else if (lociA != 0 && lociB == 0)  // 一个有坐标，另一个没有，one del
                {
                    continue;
                }
                else if (lociA == 0 && lociB != 0)  // 一个有坐标，另一个没有，one del
                {
                    continue;
                }
                else  // 都有坐标，再判断坐标的A和B是不是共线性
                {
                    if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0)  // A和B是共线性
                    {
                        synNumTmp++;  // 记为 syn

                        // debug
                        if (debugCal)
                        {
                            sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                            cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                        }
                    }
                }
            }

            // 如果和上一个的 synNum 不一样则添加
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
 * @brief 合并结果
 * 
 * @param calOutStr                 某一样品的所有计算结果
 * @param refLenMap                 染色体长度信息  map<string, length>
 * @param chrSynNumVecMap           保存最终的结果map<chr, vector<synNum> >
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
        // 变量绑定
        auto& SynNumVec = chrSynNumVecMap[chromosome];

        // 如果第一次合并的话，先初始化SynNumVec的长度为染色体长度
        if (SynNumVec.empty())
        {
            SynNumVec = vector<uint32_t>(chrLen, 0);
        }
        
        // 染色体对应的迭代器
        map<string, map<uint32_t, uint32_t> >::const_iterator iter0 = calOutStr.sampleSynOutMap.find(chromosome);
        map<string, map<uint32_t, string> >::const_iterator iter3 = calOutStr.sampleSynOutMapTmp.find(chromosome);

        // 临时迭代器，记录当前和下一个位点的得分
        map<uint32_t, uint32_t>::const_iterator iter1;
        map<uint32_t, uint32_t>::const_iterator iter2;
        
        if (iter0 != calOutStr.sampleSynOutMap.end())  // 如果含有该染色体
        {
            iter1 = iter0->second.begin();
            iter2 = iter0->second.begin();
            iter2++;
        }
        else
        {
            continue;
        }
        
        // 记录上一个节点的 synNum
        uint32_t preSynNumTmp = 0;

        // 循环染色体长度，每个位置单独添加
        for (uint32_t refLoci = 1; refLoci < chrLen + 1; refLoci++)
        {
            // 判读当前位置是否指向第一个迭代器
            if (iter0 != calOutStr.sampleSynOutMap.end() && refLoci == iter1->first)  // 如果是，赋值该迭代器对应的数量
            {
                preSynNumTmp = iter1->second;
            }
            else if (iter0 != calOutStr.sampleSynOutMap.end() && iter2 != iter0->second.end() && refLoci == iter2->first)  // 走到下一个节点，更新迭代器
            {
                iter1++;  // 更新迭代器
                iter2++;  // 更新迭代器

                preSynNumTmp = iter1->second;  // 更新坐标
            }

            // 位点频率叠加
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
 * @brief 计算
 * 
 * @param refLenMap                参考基因组长度
 * @param SynCoorTmp               共线性坐标
 * @param sampleSampleSyriMap      syri.out 输出路径字典
 * @param alignsMap                show-aligns 输出路径字典
 * @param allSynNum                所有组合数量
 * @param outputFileName           输出文件名
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
 * @brief 合并结果
 * 
 * @param chromosome                         染色体号
 * @param chrLen                             染色体长度
 * @param CALSTRUCTUREVec                    多线程输出结果
 * 
 * @return tuple<chromosome, synOutVecTmp>   tuple<chr, vector<synNum> >
**/
tuple<string, vector<uint32_t> > CALNAME::merge_fast(
    string chromosome, 
    uint32_t chrLen, 
    const vector<CALSTRUCTURE> & CALSTRUCTUREVec
)
{
    // 保存结果
    vector<uint32_t> synOutVecTmp(chrLen, 0);  //  vector<synNum>
    
    for (const auto& CALSTRUCTURETmp : CALSTRUCTUREVec) // vector<CALSTRUCTURE>
    {
        // 染色体对应的迭代器
        map<string, map<uint32_t, uint32_t> >::const_iterator iter0 = CALSTRUCTURETmp.sampleSynOutMap.find(chromosome);
        map<string, map<uint32_t, string> >::const_iterator iter3 = CALSTRUCTURETmp.sampleSynOutMapTmp.find(chromosome);

        // 临时迭代器，记录当前和下一个位点的得分
        map<uint32_t, uint32_t>::const_iterator iter1;
        map<uint32_t, uint32_t>::const_iterator iter2;
        
        if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end())  // 如果含有该染色体
        {
            iter1 = iter0->second.begin();
            iter2 = iter0->second.begin();
            iter2++;
        }
        
        // 记录上一个节点的 synNum
        uint32_t preSynNumTmp = 0;

        // 循环染色体长度，每个位置单独添加
        for (uint32_t refLoci = 1; refLoci < chrLen + 1; refLoci++)
        {
            // 判读当前位置是否指向第一个迭代器
            if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end() && refLoci == iter1->first)  // 如果是，赋值该迭代器对应的数量
            {
                preSynNumTmp = iter1->second;
            }
            else if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end() && iter2 != iter0->second.end() && refLoci == iter2->first)  // 走到下一个节点，更新迭代器
            {
                iter1++;  // 更新迭代器
                iter2++;  // 更新迭代器

                preSynNumTmp = iter1->second;  // 更新坐标
            }

            // 位点频率叠加
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
 * @brief 计算
 * 
 * @param refLenMap                参考基因组长度
 * @param SynCoorTmp               共线性坐标
 * @param sampleSampleSyriMap      syri.out 输出路径字典
 * @param alignsMap                show-aligns 输出路径字典
 * @param allSynNum                所有组合数量
 * @param outputFileName           输出文件名
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
    // 进程池
    ThreadPool pool(threadsCal);

    // 初始化线程池
    pool.init();

    // 保存多线程的结果
    /*
        map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >
        map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
    */
    vector<future<CALSTRUCTURE> > calOutStrVecTmp;
    
    // 样品名称
    for (const auto& iter1 : SynCoorTmp.sampleNameVec)  // vector<sampleName>
    {
        string sampleNameTmp = iter1;

        // 最后一个样本跳过
        if (sampleNameTmp == SynCoorTmp.sampleNameVec[SynCoorTmp.sampleNameVec.size() - 1])
        {
            continue;
        }

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate: " << sampleNameTmp << endl;

        // 多线程提交并保存结果
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

        sleep(0.1);  // 线程间隔
    }

    // 多线程结果保存
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

        // 多线程提交并保存结果
        mergeOutVec.push_back(
            pool.submit(
                merge_fast, 
                chromosome, 
                chrLen, 
                ref(CALSTRUCTUREVecTmp)
            )
        );

        sleep(0.1);  // 线程间隔
    }
    

    /* ********************************************** save the result ********************************************** */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;
    // 总的样本数
    const uint32_t sampleNum = SynCoorTmp.sampleNameVec.size();
    // // 矫正系数
    // const double correctionFactor = static_cast<double> (sampleNum - 1) / sampleNum;  // (n-1)/n

    SAVE SAVEClass(outputFileName);

    stringstream outStream; // 使用 stringstream 代替字符串拼接
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // 缓存大小为 10mb
    outStream.str().reserve(CACHE_SIZE);
    outStream << "#CHROM\tPOS\tAll_States\tSyntenic_States\tSyntenic_Diversity\n";

    for (size_t i = 0; i < mergeOutVec.size(); ++i)  // vector<future<tuple<string, vector<uint32_t> > > > 
    {
        string chromosome;
        vector<uint32_t> calOutVec;  // vector<uint32_t>

        tie(chromosome, calOutVec) = move(mergeOutVec[i].get());  // 获取多线程结果
        uint32_t loci = 1;
        for (const auto& iter1 : calOutVec)
        {
            // outStream << chromosome << "\t" << iter1.first << "\t" << allSynNum << "\t" 
            //         << iter1.second << "\t" << 1.0 - (iter1.second/(double)allSynNum)*(double)correctionFactor << "\n";

            outStream << chromosome << "\t" << loci << "\t" << allSynNum << "\t" 
                    << iter1 << "\t" << 1.0 - iter1/(double)allSynNum << "\n";

            if (outStream.tellp() >= CACHE_SIZE)  // 缓存大小为 10mb
            {
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // 清空 stringstream
                outStream.str(string());
                outStream.clear();
            }
            ++loci;
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

    // 关闭线程池
    pool.shutdown();

    return 0;
}