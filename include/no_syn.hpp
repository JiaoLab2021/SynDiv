#ifndef NO_SYN_HPP
#define NO_SYN_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include <tuple>
#include <regex>
#include <string.h>
#include <malloc.h>
#include <getopt.h>

#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "save.hpp"
#include "GzChunkReader.hpp"
#include "ThreadPool.hpp"
#include "cal.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

void help_no_syn(char* argv[]);
int main_no_syn(int argc, char* argv[]);

// debug
extern bool debugNoSyn;

namespace NOSYN {
    // open file
    KSEQ_INIT(gzFile, gzread)

    // Chromosome length of all query genomes
    class LENGTH
    {
    private:
        // chromosome length file
        vector<string> lengthsVec_;
        vector<string> lengthsTitles_;

        // Length information for all samples
        unordered_map<string, unordered_map<string, uint32_t> > sampleChrLenMap;  // map<sample, map<chr, chrLen> >
    public:
        LENGTH(
            vector<string> lengthsVec, 
            vector<string> lengthsTitles
        ) {
            lengthsVec_ = lengthsVec;
            lengthsTitles_ = lengthsTitles;
        }
        ~LENGTH() {}

        /**
         * @brief open file
         * 
         * @return void
        **/
        void index_lengths() {
            for (size_t i = 0; i < lengthsVec_.size(); i++) {
                // chromosome length file
                string lengthFileName = lengthsVec_[i];

                // title
                string sampleName = "";


                // If title is included, use the one in the list. If not, extract it based on the file name.
                if (lengthsTitles_.size() > 0) {
                    sampleName = lengthsTitles_[i];
                } else {
                    vector<string> lengthFileNameVec = split(lengthFileName, "/");  // path splitting
                    string baseName = lengthFileNameVec.back();  // file name
                    vector<string> baseNameVec = split(baseName, ".");  // Split by '.'
                    sampleName = baseNameVec[0];  // The sample name is the first element in the list
                }

                // initialize dictionary
                sampleChrLenMap[sampleName];

                // open file
                GzChunkReader GzChunkReaderClass(lengthFileName);
                // read line
                string line;
                while (GzChunkReaderClass.read_line(line)) {
                    // Skip empty lines
                    if (line.empty()) {
                        continue;
                    }

                    // split
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    // Check whether the number of columns meets the requirements
                    if (lineVec.size() != 2) {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Error: the number of columns is not 2 -> " << line << endl;
                        exit(1);
                    }

                    // Chromosome number and length information
                    string chromosome = lineVec[0];
                    uint32_t chrLen = stol(lineVec[1]);
                    
                    // storage
                    sampleChrLenMap[sampleName][chromosome] = chrLen;
                }
            }
        }

        /**
         * @brief Get chromosome length
         * 
         * @param sample Sample name
         * @param chr    chromosome
         * 
         * @return uint32_t
        **/
        uint32_t get_length(
            const string & sample, 
            const string & chr
        ) {
            // Find chromosome length
            // sample
            auto fintIter1 = sampleChrLenMap.find(sample);
            // If not, report an error
            if (fintIter1 == sampleChrLenMap.end()) {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: '" << sample << "' not present in chromosome length file." << endl;
                exit(1);
            }
            // find chromosome
            auto findIter2 = fintIter1->second.find(chr);
            if (findIter2 == fintIter1->second.end()) {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: the length of '" << chr << "' was not found -> " << sample << endl;
                exit(1);
            }

            // chromosome length
            return findIter2->second;
        }
    };


    // Make CALNAME::SYNCOOR a base class
    class NOSYNCOOR:public CALNAME::SYNCOOR
    {
    private:

    public:
        // Record non-collinear coordinates
        map<string, map<uint32_t, map<string, vector<tuple<uint32_t, uint32_t> > > > > chrStartSampleLociVecMap; // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >

        // Output string
        string outTxt;

        NOSYNCOOR() {}
        NOSYNCOOR(string fileName) {
            init(fileName);
        }
        ~NOSYNCOOR() {}

        void find_no_syn_coor(
            LENGTH & lengthClass
        );
    };
}

#endif