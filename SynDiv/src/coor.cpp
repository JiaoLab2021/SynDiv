// g++ coor.cpp -o coor -lz
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <future>
#include <malloc.h>
#include "../include/coor.hpp"
#include "../include/ThreadPool.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

void help_coor(char* argv[]);

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

    // 调试代码
    bool debug = false;


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
                debug = true;
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


    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running." << endl;

    /* ************************************ Alter The Global Variable ************************************ */
    // 更改全局参数
    COOR::change_global(
        debug
    );

    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index." << endl;
    map<string, vector<tuple<int, int, vector<string> > > > synLocSampleVecMap;  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    synLocSampleVecMap = COOR::build_syn_idx(
        synLocFileName
    );


    /* ************************************ Find Syntenic Coordinates on query genomes ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find Syntenic Coordinates on query genomes." << endl;
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
        sleep(0.00001);  // 线程间隔

        indexTmp++;  // 索引叠加
    }

    // 多线程结果保存
    vector<COOR::synAllStructure> synLocVecOutMapTmpVec;  // COOR::synAllStructure
    for (size_t i = 0; i < chrStartSynQryLocMapVec.size(); i++)
    {
        synLocVecOutMapTmpVec.push_back(move(chrStartSynQryLocMapVec[i].get()));
    }
    chrStartSynQryLocMapVec.clear();
    vector<future<COOR::synAllStructure> >().swap(chrStartSynQryLocMapVec);
    malloc_trim(0);	// 0 is for heap memory

    // 关闭线程池
    pool.shutdown();
    
    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Save the Result." << endl;
    COOR::save_result(
        synLocSampleVecMap, 
        synLocVecOutMapTmpVec, 
        outputFileName
    );

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done." << endl;

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