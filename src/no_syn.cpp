// g++ no_syn.cpp -o no_syn -lz -lpthread
#include "../include/no_syn.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// 调试代码
bool debugNoSyn = false;

int main_no_syn(int argc, char* argv[])
{
    // 共线性坐标
    string coorFileName;

    // 染色体长度文件
    vector<string> lengthsVec;
    vector<string> lengthsTitles;
    // 判断是否有 titlename
    bool haveTitles = false;

    // 输出文件名
    string outputFileName;

    //Parse command line options
    if(argc <= 2)
    {
        help_no_syn(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_no_syn(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("--coor", 6, parameterLength)) {
            if ((i+1) < argc) {
                coorFileName = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("--lengths", 9, parameterLength)) {
            if ((i+1) < argc) {
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    lengthsVec.push_back(file);
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
                    lengthsTitles.push_back(title);
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
        else if(PARAMETER_CHECK("--debug", 7, parameterLength)) {
            if ((i) < argc) {
                debugNoSyn = true;
            }
        }
    }

    if (argc <= 2) {
        help_no_syn(argv);
        exit(1);
    }
    if (coorFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name (--coor)." << endl;
        help_no_syn(argv);
        exit(1);
    }
    if (lengthsVec.size() <= 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Only " << lengthsVec.size() << " length file was specified. Nothing to extract, exiting." << endl;
        help_no_syn(argv);
        exit(1);
    }
    if ((haveTitles == true) && (lengthsVec.size() != lengthsTitles.size())) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of file titles (-n) does not match the number of files (--lengths)." << endl;
        help_no_syn(argv);
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running." << endl;


   /* ************************************ Build Chromosome Length Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Chromosome Length Index." << endl;
    NOSYN::LENGTH lengthClass(lengthsVec, lengthsTitles);
    lengthClass.index_lengths();

    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index." << endl;
    NOSYN::NOSYNCOOR NoSynCoor(coorFileName);
    NoSynCoor.open_coor();

    /* ************************************ Find no-syntenic Coordinates on query genomes ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find No-syntenic Coordinates on query genomes." << endl;
    NoSynCoor.find_no_syn_coor(
        lengthClass
    );

    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Save the Result." << endl;
    SAVE SaveClass(outputFileName);
    SaveClass.save(NoSynCoor.outTxt);  // 保存

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done." << endl;
    
    return 0;
}


void help_no_syn(char* argv[])
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " --coor FILE --lengths FILE1 FILE2 .. FILEn [options]" << endl
       << "retrieve non-syntenic coordinates" << endl
       << endl
       << "required arguments:" << endl
       << "    --coor              FILE      syntenic coordinates, output file of coor" << endl
       << "    --lengths           FILE      chromosome length files for query genomes (format: chr\tlength), one for multiple mate" << endl
       << endl
       << "optional arguments:" << endl
       << "    -n, --names         STRING    list of names to describe each file in --lengths, one for multiple mate" << endl
       << "    -o, --output        FILE      output coordinates to FILE [stdout]" << endl
       << "    --debug                       debug code (false by default)" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}


/**
 * @brief 寻找非共线性的坐标
 * 
 * @param lengthClass  样品染色体长度class
 * 
 * @return void
**/
void NOSYN::NOSYNCOOR::find_no_syn_coor(
    LENGTH & lengthClass
)
{   
    // 遍历共线性坐标
    for (auto iter1 = coorChrLociSampleLociMap.begin(); iter1 != coorChrLociSampleLociMap.end(); iter1++)  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
    {
        string chrTmp = iter1->first;  // 染色体号

        // 非共线性的数量，判断是否需要更改参考基因组的坐标，防止重复添加坐标
        uint32_t preNoSynNum = 0;

        uint32_t refStart;  // 参考基因组的起始
        uint32_t refEnd;  // 参考基因组的终止


        // 用于判断样品在小区域中的非共线性坐标是否被找到
        unordered_map<string, bool> sampleNoSynFindBoolMap;  // map<sample, bool>   true->找到了   false->没找到
        // 初始化每个样品的非共线性坐标是否被找到
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
            uint32_t refStartTmp;  // 参考基因组的起始
            uint32_t refEndTmp;  // 参考基因组的终止

            // 节点的起始和终止位置
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
            
            // 非共线性的数量，判断是否需要更改参考基因组的坐标，防止重复添加坐标
            uint32_t noSynNum = 0;

            // 如果上一个节点全是共线性的，更新参考基因组的坐标
            if (preNoSynNum  == 0)
            {
                refStart = refStartTmp;
                refEnd = refEndTmp;
            }

            for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++)  // map<sample, tuple(start, end)>
            {
                string sampleTmp = iter3->first;  // 样品名
                uint32_t StartTmp = get<0>(iter3->second);  // 起始
                uint32_t EndTmp = get<1>(iter3->second);  // 终止

                // 参考基因组的直接跳过 
                if (sampleTmp == "reference")
                {
                    continue;
                }

                // 记录为非共线性的
                if (StartTmp == 0 || EndTmp == 0)
                {
                    // 非共线性数量加1
                    noSynNum++;

                    // 找非共线性的位置
                    // 先初始化字典
                    if (chrStartSampleLociVecMap.find(chrTmp) == chrStartSampleLociVecMap.end())  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
                    {
                        chrStartSampleLociVecMap[chrTmp];
                    }
                    if (chrStartSampleLociVecMap[chrTmp].find(refStart) == chrStartSampleLociVecMap[chrTmp].end())  // map<refStart, map<sample, vector<tuple<start, end> > > >
                    {
                        chrStartSampleLociVecMap[chrTmp][refStart];
                    }

                    // 如果该sample坐标被找到了，跳过
                    if (sampleNoSynFindBoolMap[sampleTmp])  // map<sample, findBool>
                    {
                        continue;
                    }

                    // 临时坐标
                    uint32_t noSynStartTmp = 1;
                    uint32_t noSynEndTmp = 1;

                    auto iter2Tmp = iter2;  // 临时节点
                    
                    if (iter2 == iter1->second.begin())  // 如果是第一个节点
                    {
                        noSynStartTmp = 1;
                    }
                    else
                    {
                        iter2Tmp = iter2;
                        iter2Tmp--;  // 向上移动一位

                        // 如果上一个节点还是0，且不为第一个节点，继续前移
                        while (iter2Tmp != iter1->second.begin() && get<1>(iter2Tmp->second[sampleTmp]) == 0)
                        {
                            iter2Tmp--;  // 向上移动一位
                        }

                        // 赋值前一个坐标
                        if (get<1>(iter2Tmp->second[sampleTmp]) != 0)
                        {
                            noSynStartTmp = get<1>(iter2Tmp->second[sampleTmp]) + 1;
                        }
                    }

                    
                    iter2Tmp = iter2;  // 初始化迭代器
                    iter2Tmp++;  // 向下移动一位

                    // 如果下一个节点还是0，且不为最后一个节点，继续后移
                    while (iter2Tmp != iter1->second.end() && get<0>(iter2Tmp->second[sampleTmp]) == 0)
                    {
                        iter2Tmp++;  // 向下移动一位
                    }

                    // 如果是最后一个节点，则为染色体长度
                    if (iter2Tmp == iter1->second.end())
                    {
                        // 染色体长度
                        noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);
                    }
                    else
                    {
                        // 下一个节点的起始坐标 -1
                        noSynEndTmp = get<0>(iter2Tmp->second[sampleTmp]) - 1;
                    }

                    // 起始必须小于终止
                    if (noSynStartTmp < noSynEndTmp)
                    {
                        // 初始化字典
                        if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                        {
                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                        }
                        
                        chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(noSynStartTmp, noSynEndTmp));  // 记录非共线性坐标

                        sampleNoSynFindBoolMap[sampleTmp] = true;  // 记录该样品被找到了
                    }
                }
                else
                {
                    /* 判断是否为最后一个节点 */
                    auto iter2Tmp = iter2;  // 临时节点
                    iter2Tmp++;  // 向下移动一位

                    // 如果是第一个节点，判断是否从 1 开始的比对。如果和参考基因组是共线性的话，不用找染色体前半部分，
                    if (iter2 == iter1->second.begin() && StartTmp == 0)
                    {
                        // 初始化字典
                        if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                        {
                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                        }

                        chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(1, StartTmp - 1));  // 赋值
                    }
                    else if (iter2Tmp == iter1->second.end())  // 如果是最后一个节点，判断是否到达染色体结尾
                    {
                        // 染色体长度
                        uint32_t noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);

                        if (EndTmp < noSynEndTmp)
                        {
                            // 初始化字典
                            
                            if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                            {
                                chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                            }

                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(EndTmp + 1, noSynEndTmp));  // 赋值
                        }
                    }
                    
                    // 更新该样品被没找到
                    sampleNoSynFindBoolMap[sampleTmp] = false;
                }
            }
            // 更新非共线性的数量
            preNoSynNum = noSynNum;
        }
    }


    // 结果转为字符串
    for (auto iter1 : chrStartSampleLociVecMap)  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
    {
        for (auto iter2 : iter1.second)  // map<refStart, map<sample, vector<tuple<start, end> > > >
        {
            // 空的节点跳过
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