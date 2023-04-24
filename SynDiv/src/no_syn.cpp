// g++ no_syn.cpp -o no_syn -lz -lpthread
#include <iostream>
#include <vector>
#include <getopt.h>
#include "../include/cal.hpp"
#include "../include/no_syn.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

void help_no_syn(char* argv[]);

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

    // 调试代码
    bool debug = false;

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
                debug = true;
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


    /* ************************************ Alter The Global Variable ************************************ */
    NOSYN::change_global(
        debug
    );

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
    SAVE::SAVE SaveClass(outputFileName);
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
