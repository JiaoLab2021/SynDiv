// g++ -c window.cpp -o window -lz -lpthread -O3
#include "../include/window.hpp"

using namespace std;

int main_window(int argc, char* argv[])
{
    // reference genome
    string referenceFileName;

    // cal calculated score
    string calFileName;

    // window and step size
    uint32_t windowSize = 5000;
    uint32_t stepSize = 1000;

    // Output file name
    string outputFileName;

    //Parse command line options
    if(argc <= 2)
    {
        help_window(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_window(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(
            PARAMETER_CHECK("-r", 2, parameterLength) || 
            PARAMETER_CHECK("--reference", 11, parameterLength)
        )
        {
            if ((i+1) < argc) {
                referenceFileName = argv[i + 1];
                i++;
            }
        }
        else if(
            PARAMETER_CHECK("-i", 2, parameterLength) || 
            PARAMETER_CHECK("--input", 7, parameterLength)
        )
        {
            if ((i+1) < argc) {
                calFileName = argv[i + 1];
                i++;
            }
        }
        else if(
            PARAMETER_CHECK("-w", 2, parameterLength) || 
            PARAMETER_CHECK("--window", 8, parameterLength)
        )
        {
            if ((i+1) < argc) {
                windowSize = stoul(argv[i + 1]);
                i++;
            }
        }
        else if(
            PARAMETER_CHECK("-s", 2, parameterLength) || 
            PARAMETER_CHECK("--step", 6, parameterLength)
        )
        {
            if ((i+1) < argc) {
                stepSize = stoul(argv[i + 1]);
                i++;
            }
        }
        else if(
            PARAMETER_CHECK("-o", 2, parameterLength) ||
            PARAMETER_CHECK("--output", 8, parameterLength)
        )
        {
            if ((i+1) < argc) {
                outputFileName = argv[i + 1];
                i++;
            }
        }
    }

    if (argc <= 2) {
        help_window(argv);
        exit(1);
    }
    if (referenceFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '-r'." << endl;
        help_window(argv);
        exit(1);
    }
    if (calFileName.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '-i'." << endl;
        help_window(argv);
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;

    /* ************************************ Buile the Index of Reference ************************************ */
    // map<chromosome, length>
    map<string, uint32_t> refLenMap = CALNAME::build_fasta_index(
        referenceFileName
    );

    /* ************************************ calculate the average score within window ************************************ */
    Window::WINDOW WindowClass(refLenMap, calFileName, outputFileName, windowSize, stepSize);
    WindowClass.make_window();
    WindowClass.cal_open();
    WindowClass.win_count();
    WindowClass.save_result();

    return 0;
}


void help_window(char* argv[])
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -r FILE -i FILE [options]" << endl
       << "calculate the average score within window" << endl
       << endl
       << "required arguments:" << endl
       << "    -r, --reference     FILE      input FASTA reference" << endl
       << "    -i, --input         FILE      syntenic diversity, output file of 'SynDiv_c cal'" << endl
       << endl
       << "optional arguments:" << endl
       << "    -w, --window        INT       window size [5000]" << endl
       << "    -s, --step          INT       step size [1000]" << endl
       << "    -o, --output        FILE      output results to FILE [stdout]" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}
