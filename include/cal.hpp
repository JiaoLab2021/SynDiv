#ifndef CAL_HPP
#define CAL_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include <tuple>
#include <regex>
#include <memory>
#include <string.h>
#include <malloc.h>
#include <getopt.h>

#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "GzChunkReader.hpp"
#include "save.hpp"
#include "ThreadPool.hpp"


using namespace std;


// debuging code
extern bool debugCal;

// thread number
extern int threadsCal;

// Cache size for reading files
extern int32_t readBuffer;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

void help_cal(char* argv[]);
int main_cal(int argc, char* argv[]);


namespace CALNAME
{
    // kseq.h open file
    KSEQ_INIT(gzFile, gzread)

    /**
     * @brief Parse parameters
     * 
     * @param alignsTitles         aligns file title
     * @param alignsVec            aligns path Vec
     * @param syriConfigFileName   Configuration file output by syri
     * 
     * @return make_tuple(alignsMap, sampleSampleSyriMap)       map<sampleName, alignsPath>, map<sample1, map<sample2, syriOutPath> >
    */
    tuple<map<string, string>, map<string, map<string, string> > > aligns_parameter(
        const vector<string> & alignsTitles, 
        const vector<string> & alignsVec, 
        const string & syriConfigFileName
    );


    struct CALSTRUCTURE
    {
        // sampleName
        string sampleName;

        // Record the number of site collinearities
        map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >

        // Record the sample name and position corresponding to site collinearity
        map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
    };


    // Coordinate transformation
    class SYNCOOR
    {
    private:
        // Collinear coordinate file 'coor' output results
        string fileName_;

        // Determine whether 'sampleNameVec' contains all samples
        bool sampleNameVecBool = false;
    public:
        // collinear coordinates
        map<string, map<int, map<string, tuple<uint32_t, uint32_t> > > > coorChrLociSampleLociMap;  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >

        // sample name
        vector<string> sampleNameVec;

        SYNCOOR() {}
        SYNCOOR(string fileName) {
            fileName_ = fileName;
        }
        ~SYNCOOR() {}

        // Initialization interface for use by inherited classes
        void init(string fileName) {
            fileName_ = fileName;
        }

        // build index
        void open_coor() {
            // open file
            GzChunkReader GzChunkReaderClass(fileName_);

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

                // First check whether the number of columns is a multiple of three
                if (lineVec.size() % 3 != 0) {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "'Error: The number of columns is not a multiple of three -> " << fileName_ << endl;
                    exit(1);
                }

                // Record every three columns
                int idxTmp = lineVec.size() / 3;
                // temporary key 
                string chrTmp = "";
                uint32_t refStartTmp = 0;
                for (size_t i = 0; i < idxTmp; i++) {
                    string sampleNameTmp = lineVec[3*i + 0];
                    uint32_t qryStartTmp = 0;
                    uint32_t qryEndTmp = 0;

                    // Determine whether it is a number, not always 0
                    if (isdigit(lineVec[3*i + 1][0])) {
                        qryStartTmp = stoul(lineVec[3*i + 1]);
                    } if (isdigit(lineVec[3*i + 2][0])) {
                        qryEndTmp = stoul(lineVec[3*i + 2]);
                    }
                    
                    // chromosome
                    if (i == 0) {
                        chrTmp = sampleNameTmp;
                        refStartTmp = stoul(lineVec[3*i + 1]);

                        // initialize map
                        if (coorChrLociSampleLociMap.find(chrTmp) == coorChrLociSampleLociMap.end()) {
                            coorChrLociSampleLociMap[chrTmp];
                        }
                        if (coorChrLociSampleLociMap[chrTmp].find(refStartTmp) == coorChrLociSampleLociMap[chrTmp].end()) {
                            coorChrLociSampleLociMap[chrTmp][refStartTmp];
                        }

                        sampleNameTmp = "reference";
                    }

                    // Record 'sampleName', only record the first line
                    if (!sampleNameVecBool) {
                        sampleNameVec.push_back(sampleNameTmp);
                    }
                    
                    // Add to total hash table
                    coorChrLociSampleLociMap[chrTmp][refStartTmp][sampleNameTmp] = make_tuple(qryStartTmp, qryEndTmp);
                }

                // Only record the 'sampleName' of the first row
                sampleNameVecBool = true;
            }
        }
    };


    // Coordinate transformation
    class COORTRANS
    {
    private:
        // show-align output 'xx.aligns'
        string aliFileName_;

        // sampleName
        string sampleName;

        // open file
        unique_ptr<GzChunkReader> GzChunkReaderClass_;
        bool endBool = false;  // Record whether the file has been traversed

        // Log whether to loop through the next lines
        bool chrBool = true;  // Record whether the chromosome pair meets the requirements. Currently, the chromosome number is consistent. If it is inconsistent, skip it.
        bool aliBool = true;  // Record whether the start and end of the alignment are what you want
        bool findChrBool = true;  // Record whether the chromosome number you are looking for is consistent with the current chromosome number

        // Values changed while traversing the file
        // chromosome
        string chr = "";
        // The direction, start and end of alignment
        string aliRefStrand = "";
        uint32_t aliRefStart = 0;
        uint32_t aliRefEnd = 0;
        string aliQryStrand = "";
        uint32_t aliQryStart = 0;
        uint32_t aliQryEnd = 0;
        // Comparison information of the line where the file stream pointer is located
        uint32_t refStart = 0;
        string refSeq = "";
        uint32_t refEnd = 0;
        uint32_t qryStart = 0;
        string qrySeq = "";
        uint32_t qryEnd = 0;
        // The record traverses to '_ref_seq' and the corresponding coordinates
        uint32_t refLoci = 0;
        int idxTmp = -1;  // Record the index at the end of the previous find_loci() function, and the next loop can directly skip those indexes
        uint32_t qryLoci = 0;

        tuple<string, uint32_t, uint32_t> get_alignment_loc(
            string infoTmp
        );
        int next_loci();  // The coordinates of the next comparison segment
    public:
        COORTRANS() {}

        explicit COORTRANS(string aliFileName) 
            : aliFileName_(aliFileName), GzChunkReaderClass_(std::make_unique<GzChunkReader>(aliFileName_, 1024 * 1024 * readBuffer)) {}

        COORTRANS(const COORTRANS& other) = delete;  // Disable copy constructors
        COORTRANS& operator=(const COORTRANS& other) = delete;  // Disable copy assignment operator

        // move constructor
        COORTRANS(COORTRANS&& other) noexcept : 
            aliFileName_(std::move(other.aliFileName_)), 
            sampleName(std::move(other.sampleName)),
            GzChunkReaderClass_(std::move(other.GzChunkReaderClass_)),
            endBool(other.endBool),
            chrBool(other.chrBool),
            aliBool(other.aliBool),
            findChrBool(other.findChrBool),
            chr(std::move(other.chr)),
            aliRefStrand(std::move(other.aliRefStrand)),
            aliRefStart(other.aliRefStart),
            aliRefEnd(other.aliRefEnd),
            aliQryStrand(std::move(other.aliQryStrand)),
            aliQryStart(other.aliQryStart),
            aliQryEnd(other.aliQryEnd),
            refStart(other.refStart),
            refSeq(std::move(other.refSeq)),
            refEnd(other.refEnd),
            qryStart(other.qryStart),
            qrySeq(std::move(other.qrySeq)),
            qryEnd(other.qryEnd),
            refLoci(other.refLoci),
            idxTmp(other.idxTmp),
            qryLoci(other.qryLoci)
        {
            // Set the original object status to invalid
            other.GzChunkReaderClass_ = nullptr;
        }

        // move assignment operator
        COORTRANS& operator=(COORTRANS&& other) {
            if (this != &other) {
                aliFileName_ = std::move(other.aliFileName_);
                sampleName = std::move(other.sampleName);
                GzChunkReaderClass_ = std::move(other.GzChunkReaderClass_);
                endBool = other.endBool;
                chrBool = other.chrBool;
                aliBool = other.aliBool;
                findChrBool = other.findChrBool;
                chr = std::move(other.chr);
                aliRefStrand = std::move(other.aliRefStrand);
                aliRefStart = other.aliRefStart;
                aliRefEnd = other.aliRefEnd;
                aliQryStrand = std::move(other.aliQryStrand);
                aliQryStart = other.aliQryStart;
                aliQryEnd = other.aliQryEnd;
                refStart = other.refStart;
                refSeq = std::move(other.refSeq);
                refEnd = other.refEnd;
                qryStart = other.qryStart;
                qrySeq = std::move(other.qrySeq);
                qryEnd = other.qryEnd;
                refLoci = other.refLoci;
                idxTmp = other.idxTmp;
                qryLoci = other.qryLoci;
            }
            return *this;
        }
        
        // Find coordinates
        uint32_t find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };


    // Get collinear coordinates
    class SYRIOUT
    {
    private:
        // syri.out
        string fileName_ = "";

        // open file
        unique_ptr<GzChunkReader> GzChunkReaderClass_;
        bool endBool = false;  // Record whether the file has been traversed

        // Record collinearity information
        string chr = "";
        uint32_t refStart = 0;
        uint32_t refEnd = 0;
        uint32_t qryStart = 0;
        uint32_t qryEnd = 0;

        // private function
        int next_loci();  // The coordinates of the next comparison segment
    public:
        SYRIOUT() : endBool(true) {}
        explicit SYRIOUT(string fileName) : 
            fileName_(fileName), GzChunkReaderClass_(std::make_unique<GzChunkReader>(fileName_, 1024 * 1024 * readBuffer)), endBool(false) {}

        SYRIOUT(const SYRIOUT& other) = delete;  // Disable copy constructors
        SYRIOUT& operator=(const SYRIOUT& other) = delete;  // Disable copy assignment operator

        // move constructor
        SYRIOUT(SYRIOUT&& other) noexcept : 
            fileName_(std::move(other.fileName_)), 
            GzChunkReaderClass_(std::move(other.GzChunkReaderClass_)),
            endBool(other.endBool),
            chr(std::move(other.chr)),
            refStart(other.refStart),
            refEnd(other.refEnd),
            qryStart(other.qryStart),
            qryEnd(other.qryEnd) 
        {
            // Set the original object status to invalid
            other.GzChunkReaderClass_ = nullptr;
        }

        // move assignment operator
        SYRIOUT& operator=(SYRIOUT&& other) {
            if (this != &other) {
                fileName_ = std::move(other.fileName_);
                GzChunkReaderClass_ = std::move(other.GzChunkReaderClass_);
                endBool = other.endBool;
                chr = std::move(other.chr);
                refStart = other.refStart;
                refEnd = other.refEnd;
                qryStart = other.qryStart;
                qryEnd = other.qryEnd;
            }
            return *this;
        }

        int find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };


    /**
     * @brief Build a reference genome index
     * 
     * @param referenceFileName -> refgenome
     * 
     * @return refLenMap  map<chr, length>
    **/
    map<string, uint32_t> build_fasta_index(
        const string & referenceFileName
    );


    /* ************************************** calculate syntenic diversity ************************************** */
    /**
     * @brief calculate
     * 
     * @param sampleName               sample name
     * @param refLenMap                chromosome length
     * @param SynCoorTmp               collinear coordinates
     * @param sampleSampleSyriMap      syri.out output path dictionary
     * @param alignsMap                show-aligns output path dictionary
     * 
     * @return CALSTRUCTURE
    **/
    CALSTRUCTURE calculate_run(
        const string sampleName, 
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap
    );


    /* ********************************************* memory-saving mode ********************************************* */

    /**
     * @brief merge result
     * 
     * @param calOutStr                 All calculation results for a sample
     * @param refLenMap                 Chromosome length information map<string, length>
     * @param chrLociSynNumMap          Save the final result map<chr, vector<synNum> >
     * 
     * @return 0
    **/
    int merge(
        const CALSTRUCTURE& calOutStr,
        const map<string, uint32_t> & refLenMap,
        map<string, vector<uint32_t> >& chrSynNumVecMap
    );


    /**
     * @brief calculate
     * 
     * @param refLenMap                chromosome length
     * @param SynCoorTmp               collinear coordinates
     * @param sampleSampleSyriMap      syri.out output path dictionary
     * @param alignsMap                show-aligns output path dictionary
     * @param allSynNum                All combination quantities
     * @param outputFileName           Output file name
     * 
     * @return 0
    **/
    int calculate(
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap, 
        const uint32_t & allSynNum, 
        const string & outputFileName
    );


    /* ********************************************* quick mode ********************************************* */

    /**
     * @brief merge result
     * 
     * @param chromosome                         chromosome
     * @param chrLen                             chromosome length
     * @param CALSTRUCTUREVec                    Multi-threaded output results
     * 
     * @return tuple<chromosome, synOutVecTmp>   tuple<chr, vector<synNum> >
    **/
    tuple<string, vector<uint32_t> > merge_fast(
        string chromosome, 
        uint32_t chrLen, 
        const vector<CALSTRUCTURE> & CALSTRUCTUREVec
    );


    /**
     * @brief calculate
     * 
     * @param refLenMap                chromosome length
     * @param SynCoorTmp               collinear coordinates
     * @param sampleSampleSyriMap      syri.out output path dictionary
     * @param alignsMap                show-aligns output path dictionary
     * @param allSynNum                All combination quantities
     * @param outputFileName           Output file name
     * 
     * @return 0
    **/
    int calculate_fast(
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap, 
        const uint32_t & allSynNum, 
        const string & outputFileName
    );
}

#endif