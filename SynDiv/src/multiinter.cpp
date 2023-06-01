// g++ multiinter.cpp -o multiinter -lz
#include <iostream>
#include <cstring>
#include <vector>
#include <getopt.h>
#include "../include/multiinter.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

void help_multiinter(char* argv[]);

int main_multiinter(int argc, char* argv[])
{
    // 输入文件列表和名称  syri.out
    vector<string> inputFiles;
    vector<string> inputTitles;

    // 判断是否有 titlename
    bool haveTitles = false;

    // 输出文件名
    string outputFileName;

    // 调试代码
    bool debug = false;


    //Parse command line options
    if(argc <= 2)
    {
        help_multiinter(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_multiinter(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-i", 2, parameterLength)) || 
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
        else if(PARAMETER_CHECK("--debug", 7, parameterLength)) {
            if ((i) < argc) {
                debug = true;
            }
        }
    }

    if (argc <= 2 || inputFiles.size() <= 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Only " << inputFiles.size() << " syri.out file was specified. Nothing to extract, exiting." << endl;
        help_multiinter(argv);
        return 1;
    }
    if ((haveTitles == true) && (inputFiles.size() != inputTitles.size())) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of file titles (-n) does not match the number of files (-i)." << endl;
        help_multiinter(argv);
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running." << endl;


    /* ************************************ Alter The Global Variable ************************************ */
    // 是否调试代码
    MULTIINTER::debug = debug;

    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index." << endl;
    // 存储所有line的共线性信息
    map<string, map<string, vector<pair<int, int> > > > chrLineSynVecMap;  // map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
    
    int indexTmp = 0;  // 记录 lineName 的索引
    for(auto it1 : inputFiles)
    {
        string lineName;
        map<string, vector<pair<int, int> > > chrSynVecMap;
        tie(lineName, chrSynVecMap) = MULTIINTER::build_syn_index(it1);

        // 如果有 inputTitles 重新赋值
        if (haveTitles)
        {
            lineName = inputTitles[indexTmp];
        }
        
        for(auto it2 : chrSynVecMap)
        {
            string chromosome = it2.first;
            chrLineSynVecMap[chromosome][lineName] = it2.second;
        }

        indexTmp++;  // 索引叠加
    }

    /* ************************************ Find intersection ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find intersection." << endl;
    map<string, map<int, vector<tuple<int, vector<string> > > > > outChrStartEndLineVecMap;  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>
    map<int, string> idxLineMap;
    tie(outChrStartEndLineVecMap, idxLineMap) = MULTIINTER::syn_multiinter_find(
        chrLineSynVecMap
    );

    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Saving the Result." << endl;
    MULTIINTER::save_result(
        outChrStartEndLineVecMap, 
        idxLineMap, 
        outputFileName
    );

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done." << endl;
    
    return 0;
}


void help_multiinter(char* argv[])
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i FILE1 FILE2 .. FILEn [options]" << endl
       << "identifies common intervals (SYNAL) among multiple syri.out files" << endl
       << endl
       << "required arguments:" << endl
       << "    -i, --inputs        FILE      list of output files of syri (syri.out), one for multiple mate" << endl
       << endl
       << "optional arguments:" << endl
       << "    -n, --names         STRING    list of names to describe each file in -i, one for multiple mate" << endl
       << "    -o, --output        FILE      output syntenic intersection to FILE [stdout]" << endl
       << "    --debug                       debug code (false by default)" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}