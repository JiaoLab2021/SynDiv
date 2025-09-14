// g++ SynDiv_c.cpp src/*.cpp -o SynDiv_c -lz -lpthread -O3 -std=c++17 -march=native
#include <iostream>
#include <vector>
#include "zlib.h"
#include <map>
#include <getopt.h>

#include "include/get_time.hpp"
#include "include/strip_split_join.hpp"
#include "include/sys.hpp"
#include "include/multiinter.hpp"
#include "include/coor.hpp"
#include "include/no_syn.hpp"
#include "include/cal.hpp"
#include "include/window.hpp"
#include "include/Fst.hpp"

using namespace std;

// define data
#define PROGRAM_DATA "2025/09/13"
// define version
#define PROGRAM_VERSION "1.1.5-r8"
// define author
#define PROGRAM_AUTHOR "Zezhen Du"
// define E-mail
#define PROGRAM_E_MAIL "dzz0539@gmail.com or dzz0539@163.com"


void help(char** argv);

int main(int argc, char** argv)
{
    // Record initial time
    double realtime0 = realtime();

    // Print help documentation
    if (argc == 1) {
        help(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "You are using SynDiv (v" << PROGRAM_VERSION << ")\n\n";

    // subcommand
    string subcommand = argv[1];

    if (subcommand == "-h" || subcommand == "--help") {
        help(argv);
        return 1;
    } else if (subcommand == "multiinter") {
        main_multiinter(argc, argv);
    } else if (subcommand == "coor") {
        main_coor(argc, argv);
    } else if (subcommand == "no_syn") {
        main_no_syn(argc, argv);
    } else if (subcommand == "cal") {
        main_cal(argc, argv);
    } else if (subcommand == "window") {
        main_window(argc, argv);
    } else if (subcommand == "Fst") {
        main_Fst(argc, argv);
    } else {
        cerr << "Error: ["<< argv[0] << "] command " << subcommand << " not found." << endl;
        help(argv);
        return 1;
    }

    // Print time and memory usage
    fprintf(stderr, "[SynDiv::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    
    return 0;
}

// help documentation
void help(char** argv) {
  cerr << "usage: " << argv[0] << " <command> [options]" << endl
       << endl
       << "data: " << PROGRAM_DATA << endl
       << "version: " << PROGRAM_VERSION << endl
       << "author: " << PROGRAM_AUTHOR << endl
       << endl
       << "subcommands:" << endl
       << "   multiinter     identifies common intervals (SYNAL) among multiple syri.out files" << endl
       << "   coor           retrieve syntenic coordinates on the query genome" << endl
       << "   no_syn         retrieve non-syntenic coordinates" << endl
       << "   cal            compute syntenic diversity" << endl
       << "   window         calculate the average score within window" << endl
       << "   Fst            calculate population differentiation index Fst" << endl
       << endl
       << "   -h, --help     print this help document" << endl
       << endl
       << "If you encounter any issues related to the code, please don't hesitate to contact us via email at " << PROGRAM_E_MAIL << "." << endl;
}