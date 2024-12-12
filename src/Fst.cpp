// g++ -c Fst.cpp -o Fst -lz -lpthread -O3
#include "../include/Fst.hpp"

using namespace std;

int main_Fst(int argc, char* argv[]) {
    // Syntenic diversity results files (SynDiv_c cal)
    vector<string> calFileVec;

    // Output file name
    string outputFileName;

    //Parse command line options
    if(argc <= 2) {
        help_Fst(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_Fst(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if (
            PARAMETER_CHECK("-c", 2, parameterLength) ||
            PARAMETER_CHECK("--cal", 5, parameterLength)
        ) {
            if ((i+1) < argc) {
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    calFileVec.push_back(file);
                    i++;
                    if (i < argc)
                        file = argv[i];
                }
                i--;
            }
        } else if (
            PARAMETER_CHECK("-o", 2, parameterLength) ||
            PARAMETER_CHECK("--output", 8, parameterLength)
        ) {
            if ((i+1) < argc) {
                outputFileName = argv[i + 1];
                i++;
            }
        }
    }

    if (argc <= 2) {
        help_Fst(argv);
        exit(1);
    }
    if (calFileVec.size() < 3) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The -c parameter requires at least three files, corresponding to collinearity diversity results for Subgroup 1, Subgroup 2, and the Total Group." << endl;
        help_Fst(argv);
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;

    /* ************************************ calculate population differentiation index Fst ************************************ */
    Fst::Fst FstClass(calFileVec, outputFileName);
    // Get the number of populations
    FstClass.get_number();
    // Open file
    FstClass.open_file();
    // Calculate Fst
    FstClass.cal_Fst();

    return 0;
}


void help_Fst(char* argv[]) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " -c FILE1 FILE2 .. FILEn [options]" << endl
       << "calculate population differentiation index Fst" << endl
       << endl
       << "required arguments:" << endl
       << "    -c, --cal           FILE      list of output files for syntenic diversity (output file from 'SynDiv_c cal'), one for multiple mates" << endl
       << endl
       << "optional arguments:" << endl
       << "    -o, --output        FILE      output results to FILE [stdout]" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}
