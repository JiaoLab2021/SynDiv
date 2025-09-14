#ifndef FST_HPP
#define FST_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include <numeric>
#include <tuple>
#include <string.h>
#include <getopt.h>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <iterator>

#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "save.hpp"
#include "GzChunkReader.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

void help_Fst(char* argv[]);
int main_Fst(int argc, char* argv[]);


namespace Fst {

// Calculate population differentiation index Fst
class Fst {
private:
    // Syntenic diversity results files (SynDiv_c cal)
    vector<string> calFileVec_;

    // Output file
    string outputFileName_;

    // Population number
    map<uint32_t, vector<string> > nFileVecMap_;  // map<population number, vector<syntenic diversity file> >

    vector<string> fileVec_;  // File name
    vector<uint32_t> numberVec_;  // Number of populations
    vector<std::unique_ptr<GzChunkReader> > gzReaderVec_;  // File handle


    /**
     * @brief Get the number of populations
     * 
     * @return uint32_t
     */
    uint32_t calculateN_(uint32_t states_number) {
        double discriminant = 1 + 8.0 * states_number;
        if (discriminant < 0) {
            throw std::invalid_argument("Invalid states_number: discriminant is negative.");
        }
        double n = (1 + std::sqrt(discriminant)) / 2;
        if (n != static_cast<uint32_t>(n)) {
            throw std::invalid_argument("states_number does not result in an integer n.");
        }
        return static_cast<uint32_t>(n);
    }

public:
    Fst(
        const vector<string>& calFileVec, 
        string outputFileName = ""
    ) {
        calFileVec_ = calFileVec;
        outputFileName_ = outputFileName;
    }
    ~Fst() {}


    /**
     * @brief Get the number of populations
     * 
     * @return int
     */
    int get_number() {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Retrieve the genome count for each subgroup ..." << endl << endl;

        for (const auto & calFile : calFileVec_) {
            // open file
            GzChunkReader GzChunkReaderClass(calFile);

            // read line
            string line;
            while (GzChunkReaderClass.read_line(line)) {
                // Skip empty lines
                if (line.empty() || line.find("#") != string::npos) {
                    continue;
                }

                // split
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // Check the number of columns
                if (lineVec.size() != 5) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The syntenic diversity file '"<< calFile << "' must have five columns." << endl;
                    exit(1);
                }

                // Get the number of states
                uint32_t states_number = stoul(lineVec[2]);

                // Calculate the number of populations
                uint32_t n = calculateN_(states_number);

                // Check if n is smaller than 0
                if (n <= 0) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of populations in '" << calFile << "' is " << n << ", which is less than 0." << endl;
                    exit(1);
                }

                // Save the number of populations
                nFileVecMap_[n].push_back(calFile);

                // Print the number of populations
                cerr << "           - The number of populations in '" << calFile << "' is " << n << endl;

                break;
            }
        }
        
        cerr << endl;

        return 0;
    }


    /** *
     * @brief Open file
     * 
     * @return 0
     */
    int open_file() {
        // Number of genomes and files for the total population
        auto maxElement = std::prev(nFileVecMap_.end());
        uint32_t nPopulation = maxElement->first;
        vector<string> populationFileVec = maxElement->second;
        // Check if the number of files under the current group number is greater than 1, if so, report an error
        if (populationFileVec.size() > 1) {
            cerr << "[" << __func__ << "::" << getTime() << "] Error: " << nPopulation 
                << " is the total number of population genomes, and there should be exactly 1 syntenic diversity file for this number, but currently there are " 
                << populationFileVec.size() << "." << endl;
            exit(1);
        }

        // Subpopulation number and file
        for (auto [n, fileVec] : nFileVecMap_) {
            if (n == nPopulation) {
                continue;
            }
            
            for (auto file : fileVec) {
                fileVec_.push_back(file);
                numberVec_.push_back(n);
                gzReaderVec_.push_back(std::make_unique<GzChunkReader>(file));
            }
        }

        // Population number and file
        fileVec_.push_back(populationFileVec[0]);
        numberVec_.push_back(nPopulation);
        gzReaderVec_.push_back(std::make_unique<GzChunkReader>(populationFileVec[0]));

        return 0;
    }


    /**
     * @brief Read line
     * 
     * @param syndivFileName    Syntenic diversity file name
     * @param syndivFileReader  Syntenic diversity file reader
     * 
     * @return tuple<chromosome, location, syndiv>
     */
    tuple<string, uint32_t, double> read_file(string syndivFileName, std::unique_ptr<GzChunkReader>& syndivFileReader) {
        string line;
        while (syndivFileReader->read_line(line)) {
            // Skip empty lines
            if (line.empty() || line.find("#") != string::npos) {
                continue;
            }

            // split
            std::istringstream iss(line);
            vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

            // Check the number of columns
            if (lineVec.size() != 5) {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The syntenic diversity file '"<< syndivFileName << "' must have five columns." << endl;
                exit(1);
            }

            // Location
            string chromosome = lineVec[0];
            uint32_t location = stoul(lineVec[1]);
            // Syntenic diversity
            double syndiv = stod(lineVec[4]);

            return make_tuple(chromosome, location, syndiv);
        }

        return make_tuple("", 0, 0.0);
    }


    /** *
     * @brief Calculate Fst
     * 
     * @return 0
     */
    int cal_Fst() {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating Fst ..." << endl;
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Results will be saved to '" << outputFileName_ << "' ..." << endl;

        /* ************************************ Save result ************************************ */
        SAVE SAVEClass(outputFileName_);

        stringstream outStream; // Use stringstream instead of string concatenation
        static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
        outStream.str().reserve(CACHE_SIZE);
        outStream << "#CHROM\tPOS\t";
        for (size_t i = 0; i < fileVec_.size() - 1; i++) {
            outStream << fileVec_[i] << "\t";
        }
        outStream << "Population_SynDiv\tSyn_Fst\n";

        /* ************************************ Calculate Fst ************************************ */
        uint32_t populationNum = numberVec_.size();  // The number of populations
        uint32_t nPopulation = numberVec_.back();  // Total number of populations
        vector<tuple<string, uint32_t, double> > chrLocSynDivTupVec;  // tuple<chromosome, location, syntenic diversity>
        chrLocSynDivTupVec.resize(populationNum - 1);
        string line;

        while (gzReaderVec_.back()->read_line(line)) {
            // Skip empty lines
            if (line.empty() || line.find("#") != string::npos) {
                continue;
            }

            // split
            std::istringstream iss(line);
            vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

            // Check the number of columns
            if (lineVec.size() != 5) {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The syntenic diversity file '"<< fileVec_.back() << "' must have five columns." << endl;
                exit(1);
            }

            // Location
            string populationChr = lineVec[0];
            uint32_t populationLoc = stoul(lineVec[1]);
            // Syntenic diversity
            double populationSynDiv = stod(lineVec[4]);

            // Syntenic diversity vector
            vector<double> synDivVec(populationNum, 0.0);
            
            // Subpopulation
            for (size_t i = 0; i < gzReaderVec_.size() - 1; i++) {
                // Chromosome || Location
                while (
                    populationChr != get<0>(chrLocSynDivTupVec[i]) || 
                    populationLoc > get<1>(chrLocSynDivTupVec[i])
                ) {
                    auto [subChr, subLoc, subSynDiv] = read_file(fileVec_[i], gzReaderVec_[i]);
                    // Fix issue #8: If the file has been read, exit the current while loop
                    if (subChr.empty()) {  // EOF
                        chrLocSynDivTupVec[i] = make_tuple("", 0u, 0.0);
                        synDivVec[i] = 0.0;
                        break;
                    }
                    chrLocSynDivTupVec[i] = make_tuple(subChr, subLoc, subSynDiv);
                    synDivVec[i] = subSynDiv;
                }
            }

            synDivVec.back() = populationSynDiv;

            // Calculate Fst
            // Components
            double numerator = 0.0, denominator = 0.0;
            double weightedSubPopTerm = 0.0, subPopTermSum = 0.0;

            // Total Population Term
            double totalTerm = 2 * populationSynDiv * (1 - populationSynDiv);

            // Subpopulation Contributions
            for (size_t i = 0; i < chrLocSynDivTupVec.size(); ++i) {
                auto [subChr, subLoc, subSynDiv] = chrLocSynDivTupVec[i];
                if (subChr != populationChr || subLoc != populationLoc) {
                    continue;
                }

                double subTerm = 2 * subSynDiv * (1 - subSynDiv);
                weightedSubPopTerm += (subTerm * numberVec_[i]) / nPopulation;
                subPopTermSum += subTerm;
            }

            // Numerator and Denominator
            numerator = totalTerm - weightedSubPopTerm;
            denominator = totalTerm;

            double Fst = 0.0;
            if (totalTerm != 0) {
                Fst = (numerator / denominator) + (totalTerm - subPopTermSum);
                if (Fst <= 0) {
                    Fst = 0;
                }
            }

            /* ************************************ Save result ************************************ */
            outStream << populationChr << "\t" << populationLoc << "\t";
            for (size_t i = 0; i < synDivVec.size(); i++) {
                outStream << synDivVec[i] << "\t";
            }
            outStream << Fst << "\n";

            if (outStream.tellp() >= CACHE_SIZE) {  // Cache size is 10mb
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // Clear stringstream
                outStream.str(string());
                outStream.clear();
            }
        }

        /* ************************************ Save result ************************************ */
        if (outStream.tellp() > 0) {  // Write for the last time
            string outTxt = outStream.str();
            SAVEClass.save(outTxt);
            // Clear stringstream
            outStream.str(string());
            outStream.clear();
        }

        return 0;
    }
};
}

#endif