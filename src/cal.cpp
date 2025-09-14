// g++ -c cal.cpp -o cal -lz -lpthread -O3 -std=c++17
#include "../include/cal.hpp"

using namespace std;

// debug code
bool debugCal = false;

// Threads
int threadsCal = 30;

// Cache size for reading files
int32_t readBuffer = 1;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

void help_cal(char* argv[]);

int main_cal(int argc, char* argv[]) {
    // reference genome
    string referenceFileName;

    // collinear coordinates
    string coorFileName;

    // no-syn coordinates
    string noSynFileName;

    // Configuration file output by syri
    string syriConfigFileName;

    // aligns file with the reference genome
    vector<string> alignsVec;
    vector<string> alignsTitles;
    // Determine whether there is a titlename
    bool haveTitles = false;

    // output filename
    string outputFileName;

    // Whether the software runs in express mode
    bool fastBool = false;

    //Parse command line options
    if(argc <= 2) {
        help_cal(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_cal(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if ((PARAMETER_CHECK("-r", 2, parameterLength)) || 
        (PARAMETER_CHECK("--reference", 11, parameterLength))) {
            if ((i+1) < argc) {
                referenceFileName = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("--coor", 6, parameterLength)) {
            if ((i+1) < argc) {
                coorFileName = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("--no_syn", 8, parameterLength)) {
            if ((i+1) < argc) {
                noSynFileName = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("--syri_outs", 11, parameterLength)) {
            if ((i+1) < argc) {
                syriConfigFileName = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("--aligns", 8, parameterLength)) {
            if ((i+1) < argc) {
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    alignsVec.push_back(file);
                    i++;
                    if (i < argc)
                        file = argv[i];
                }
                i--;
            }
        } else if ((PARAMETER_CHECK("-n", 2, parameterLength)) || 
        (PARAMETER_CHECK("--names", 7, parameterLength))) {
            if ((i+1) < argc) {
                haveTitles = true;
                i = i+1;
                string title = argv[i];
                while (title[0] != '-' && i < argc) {
                    alignsTitles.push_back(title);
                    i++;
                    if (i < argc)
                        title = argv[i];
                }
                i--;
            }
        } else if (PARAMETER_CHECK("-o", 2, parameterLength) ||
        (PARAMETER_CHECK("--output", 8, parameterLength))) {
            if ((i+1) < argc) {
                outputFileName = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("-t", 2, parameterLength) ||
        (PARAMETER_CHECK("--threads", 9, parameterLength))) {
            if ((i+1) < argc) {
                threadsCal = stoi(argv[i + 1]);
                i++;
            }
        } else if (PARAMETER_CHECK("--buffer", 8, parameterLength)) {
            if ((i+1) < argc) {
                readBuffer = stoul(argv[i + 1]);
                i++;
            }
        } else if (PARAMETER_CHECK("--fast", 6, parameterLength)) {
            if ((i) < argc) {
                fastBool = true;
            }
        } else if (PARAMETER_CHECK("--debug", 7, parameterLength)) {
            if ((i) < argc) {
                debugCal = true;
            }
        }
    }

    if (argc <= 2) {
        help_cal(argv);
        exit(1);
    }
    if (threadsCal <= 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Threads must be greater than 0 (-t)." << endl;
        help_cal(argv);
        exit(1);
    }
    if (readBuffer <= 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: ReadBuffer must be greater than 0 (--buffer)." << endl;
        help_cal(argv);
        exit(1);
    }
    if (referenceFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '-r'." << endl;
        help_cal(argv);
        exit(1);
    }
    if (coorFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '--coor'." << endl;
        help_cal(argv);
        exit(1);
    }
    if (noSynFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '--no_syn'." << endl;
        help_cal(argv);
        exit(1);
    }
    if (syriConfigFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name '--syri_outs'." << endl;
        help_cal(argv);
        exit(1);
    }
    if (alignsVec.size() <= 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Only " << alignsVec.size() << " aligns file was specified. Nothing to extract, exiting." << endl;
        help_cal(argv);
        exit(1);
    }
    if ((haveTitles == true) && (alignsVec.size() != alignsTitles.size())) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of file titles (-n) does not match the number of files (--aligns)." << endl;
        help_cal(argv);
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;


    /* ************************************ Change Global Parameters ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Change Global Parameters ..." << endl;
    // print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "debug: " << debugCal << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "threads: " << threadsCal << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "readBuffer: " << readBuffer << " MB" << endl;

    // If debugging the code, the number of threads is set to 1
    if (debugCal) { threadsCal = 1; }

    /* ************************************ Parse Parameters ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Parse Parameters ..." << endl;
    // Get the aligns/syri.out file dictionary
    map<string, string> alignsMap;  // map<sampleName, alignsPath>
    unordered_map<string, unordered_map<uint32_t, unordered_map<string, tuple<uint32_t, uint32_t> > > > insChrStartSampleLociTupMap;  // map<chr, map<start, map<sample, tuple<start, length> > > >
    map<string, map<string, string> > sampleSampleSyriMap;  // map<sample1, map<sample2, syriOutPath> >
    tie(alignsMap, insChrStartSampleLociTupMap, sampleSampleSyriMap) = CALNAME::cal_parameter(
        alignsTitles, 
        alignsVec, 
        noSynFileName, 
        syriConfigFileName
    );

    /* ************************************ Building Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building Syntenic Coordinates Index ..." << endl;
    CALNAME::SYNCOOR SynCoorTmp(coorFileName);
    SynCoorTmp.open_coor();

    // All combination quantities
    uint32_t allSynNum = SynCoorTmp.sampleNameVec.size() * (SynCoorTmp.sampleNameVec.size() - 1) / 2;  // n*(n-1)/2

    /* ************************************ Building Reference Genome Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building Reference Genome Index ..." << endl;
    map<string, uint32_t> refLenMap = CALNAME::build_fasta_index(
        referenceFileName
    );

    /* ************************************ Calculating Syntenic Diversity ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating Syntenic Diversity ..." << endl;
    if (fastBool) {
        CALNAME::calculate_fast(
            refLenMap, 
            SynCoorTmp, 
            sampleSampleSyriMap, 
            alignsMap, 
            insChrStartSampleLociTupMap, 
            allSynNum, 
            outputFileName
        );
    } else {
        CALNAME::calculate(
            refLenMap, 
            SynCoorTmp, 
            sampleSampleSyriMap, 
            alignsMap, 
            insChrStartSampleLociTupMap, 
            allSynNum, 
            outputFileName
        );
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ..." << endl;
    
    return 0;
}


void help_cal(char* argv[]) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " -r FILE --coor FILE --no_syn FILE --syri_outs FILE --aligns FILE1 FILE2 .. FILEn [options]" << endl
       << "compute syntenic diversity" << endl
       << endl
       << "required arguments:" << endl
       << "    -r, --reference     FILE      input FASTA reference" << endl
       << "    --coor              FILE      syntenic coordinates, output by 'SynDiv_c coor'" << endl
       << "    --no_syn            FILE      non-syntenic coordinates, output by 'SynDiv_c no_syn'" << endl
       << "    --syri_outs         FILE      config file for alignment, output by 'SynDiv_p' (format: sample1\tsample2\tsyri.out)" << endl
       << "    --aligns            FILE      list of output files of show-aligns (.aligns), one for multiple mate" << endl
       << endl
       << "optional arguments:" << endl
       << "    -n, --names         STRING    list of names to describe each file in --aligns, one for multiple mate" << endl
       << "    -o, --output        FILE      output results to FILE [stdout]" << endl
       << "    -t, --threads       INT       number of compute threads to use [30]" << endl
       << "    --buffer            INT       buffer size for file reading, measured in MB [1]" << endl
       << "    --fast                        enabling quick mode will increase memory consumption (false by default)" << endl
       << "    --debug                       debug code (false by default)" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}


/**
 * @brief Parse parameters (aligns + no_syn + no_syn_alignment)
 * 
 * @param alignsTitles         aligns file title
 * @param alignsVec            aligns path Vec
 * @param noSynFileName        'SynDiv_c no_syn' output file name
 * @param syriConfigFileName   'SynDiv_p' output file name
 * 
 * @return make_tuple(alignsMap, insChrStartSampleLociTupMap, sampleSampleSyriMap)       map<sampleName, alignsPath>, map<chr, map<start, map<sample, tuple<start, length> > > >, map<sample1, map<sample2, syriOutPath> >
*/
tuple<map<string, string>, unordered_map<string, unordered_map<uint32_t, unordered_map<string, tuple<uint32_t, uint32_t> > > >, map<string, map<string, string> > > CALNAME::cal_parameter(
    const vector<string> & alignsTitles, 
    const vector<string> & alignsVec, 
    const string & noSynFileName,
    const string & syriConfigFileName
) {
    /* ************************************* aligns ************************************* */
    map<string, string> alignsMap;  // map<sampleName, alignsPath>

    uint32_t indexTmp = 0;  // Index of record title
    std::regex reg1("\\.aligns(\\.gz|\\.GZ)?$", std::regex_constants::icase);  // Optimized regex

    for (const auto& alignsPathTmp : alignsVec) {
        string ailgnsTitleTmp;

        // If it contains title
        if (!alignsTitles.empty()) {
            ailgnsTitleTmp = alignsTitles[indexTmp];
        } else {
            size_t lastSlashPos = alignsPathTmp.find_last_of("/");  // Optimized string operation
            ailgnsTitleTmp = alignsPathTmp.substr(lastSlashPos + 1);
            ailgnsTitleTmp = regex_replace(ailgnsTitleTmp, reg1, "");  // Optimized regex use
        }

        // Assignment
        alignsMap[ailgnsTitleTmp] = alignsPathTmp;

        indexTmp++;  // Index iteration
    }

    /* ************************************* SynDiv_c no_syn ************************************* */
    unordered_map<string, unordered_map<uint32_t, unordered_map<string, tuple<uint32_t, uint32_t> > > > insChrStartSampleLociTupMap;  // map<chr, map<start, map<sample, tuple<start, length> > > >

    // open file
    GzChunkReader GzChunkReaderClassForSynDiv_c(noSynFileName);

    // read line
    string line;
    while (GzChunkReaderClassForSynDiv_c.read_line(line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // split
        std::istringstream iss(line);
        vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

        string chromosome = lineVec[0];
        uint32_t refStart = isdigit(lineVec[1][0]) ? stoul(lineVec[1]) : 0;

        if (refStart == 0) continue;

        auto& sampleLociTupMap = insChrStartSampleLociTupMap[chromosome][refStart];

        for (uint32_t i = 2; i < lineVec.size(); i++) {
            vector<string> sampleNameLociVec = split(lineVec[i], ":");
            if (sampleNameLociVec.size() < 2) continue;  // Skip
            string& sampleName = sampleNameLociVec[0];
            vector<string> regions = split(sampleNameLociVec[1], ";");

            uint32_t minStart = UINT32_MAX;
            uint32_t sampleLen = 0;

            for (const auto& region : regions) {
                vector<string> lociVec = split(region, "-");
                if (lociVec.size() < 2) continue;  // Skip malformed regions
                uint32_t qryStart = isdigit(lociVec[0][0]) ? stoul(lociVec[0]) : 0;
                uint32_t qryEnd = isdigit(lociVec[1][0]) ? stoul(lociVec[1]) : 0;

                if (qryStart < minStart) minStart = qryStart;
                sampleLen += qryEnd - qryStart + 1;
            }

            sampleLociTupMap[sampleName] = make_tuple(minStart, sampleLen);
        }
    }

    /* ************************************* SynDiv_p ************************************* */
    map<string, map<string, string>> sampleSampleSyriMap;  // map<sample1, map<sample2, syriOutPath> >

    // open file
    GzChunkReader GzChunkReaderClassForSynDiv_p(syriConfigFileName);

    // read line
    line.clear();
    while (GzChunkReaderClassForSynDiv_p.read_line(line)) {
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        
        // split
        std::istringstream iss(line);
        vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

        // Check if listed is 3
        if (lineVec.size() != 3) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: The dataframe in '" << syriConfigFileName << "' does not contain three columns as expected." << endl;
            exit(1);
        }

        // Sample and route information
        const string& sample1 = lineVec[0];
        const string& sample2 = lineVec[1];
        const string& syriPath = lineVec[2];

        sampleSampleSyriMap[sample1][sample2] = syriPath;
    }

    return make_tuple(alignsMap, insChrStartSampleLociTupMap, sampleSampleSyriMap);
}


/**
 * @brief Parse '-- BEGIN alignment ' field
 * 
 * @param infoTmp   Delete '-- BEGIN alignment ' and then press ' | ' to separate the fields, '+1 278 - 1703'
 * 
 * @return tuple<strand, start, end>
*/
tuple<string, uint32_t, uint32_t> CALNAME::COORTRANS::get_alignment_loc(string infoTmp) {
    static const std::regex pat(R"(([+-])1\s+(\d+)\s*-\s*(\d+))");
    std::smatch m;
    if (!std::regex_search(infoTmp, m, pat) || m.size() != 4) {
        std::cerr << "[" << __func__ << "::" << getTime() << "] Error: can't parse alignment loc: \"" << infoTmp << "\"\n";
        std::exit(1);
    }
    std::string strand = m[1].str();                // "+" or "-"
    uint32_t a = (uint32_t)std::stoul(m[2].str());  // start
    uint32_t b = (uint32_t)std::stoul(m[3].str());  // end
    if (strand == "+") return {strand, a, b};
    else               return {strand, b, a};
}

/**
 * @brief The coordinates of the next comparison segment
 * 
 * @return int  0
**/
int CALNAME::COORTRANS::next_loci() {
    if (!endBool) {
        // read line
        string line;
        if (GzChunkReaderClass_->read_line(line)) {
            // Skip empty lines
            if (line.empty()) {
                return 0;
            }

            // new chromosome pair
            if (line.find("-- Alignments between") != string::npos) {  // New aligned chromosomes    -- Alignments between chr1 and chr1
                // split string
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                string refChr = lineVec[3];  // Temporary value for ref chromosome
                string qryChr = lineVec[5];  // Temporary value for qry chromosome

                if (refChr != qryChr) {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Alignment skipped due to mismatch in chromosome names -> " << refChr << " != " << qryChr << endl;
                    chrBool = false;  // Record chromosome does not meet regulations
                } else {
                    chrBool = true;  // Record chromosome compliance
                }

                // This judgment prevents losing the first "-- BEGIN alignment"
                if (chr.length() == 0) {
                    findChrBool = true;  // the chromosome is required
                }
                
                // alignment chromosome information update
                chr = refChr;
            } else if (line.find("-- BEGIN alignment ") != string::npos && chrBool && findChrBool) {  // New aligned chromosomes   -- BEGIN alignment [ +1 1078 - 68996 | +1 1 - 67961 ], Chromosomes must conform to regulations
                line = strip(line.erase(0, 21), '\n');  // Delete '-- BEGIN alignment [ ' 21 characters in total
                line = line.erase(line.size()-2, 2);  //Delete ' ]', a total of 2 characters, located in the last two characters

                vector<string> lineVec = split(line, " | ");  // split '+1 278 - 1703 | -1 2148 - 751'

                // Get information about alignment direction, start and end, ref
                tie(aliRefStrand, aliRefStart, aliRefEnd) = get_alignment_loc(
                    lineVec[0]
                );

                // Get information about alignment direction, start and end, qry
                tie(aliQryStrand, aliQryStart, aliQryEnd) = get_alignment_loc(
                    lineVec[1]
                );
            } else if (isdigit(line[0]) != 0 && chrBool && findChrBool && aliBool) {
                // Loop only on lines starting with a number and containing syn
                // The chromosome must comply with the regulations, and the start and end of the comparison must comply with the regulations.
                /* ***************************** ref ***************************** */
                // split
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // Reset coordinates (this)
                refSeq = lineVec[1];
                // Calculate length without '.'
                uint32_t refSeqLen = std::count_if(refSeq.begin(), refSeq.end(), [](char c){ return c != '.'; });
                if (aliRefStrand == "+") {  // forward comparison
                    refStart = stoul(lineVec[0]);
                    refEnd = refStart + refSeqLen - 1;
                    refLoci = refStart - 1;  // Refresh coordinates
                } else {  // reverse comparison
                    refEnd = stoul(lineVec[0]);
                    refStart = refEnd - refSeqLen + 1;
                    refLoci = refEnd + 1;  // Refresh coordinates
                }

                /* ***************************** qry ***************************** */
                // Read next line
                GzChunkReaderClass_->read_line(line);

                // split
                std::istringstream iss1(line);
                vector<string> infoVecTmp(std::istream_iterator<std::string>{iss1}, std::istream_iterator<std::string>());
                lineVec.assign(infoVecTmp.begin(), infoVecTmp.end());

                // Refresh coordinates (this)
                qrySeq = lineVec[1];
                // Calculate length without '.'
                uint32_t qrySeqLen = std::count_if(qrySeq.begin(), qrySeq.end(), [](char c){ return c != '.'; });
                if (aliQryStrand == "+") {  // forward comparison
                    qryStart = stoul(lineVec[0]);
                    qryEnd = qryStart + qrySeqLen - 1;
                    qryLoci = qryStart - 1;  // Refresh coordinates
                } else {  // reverse comparison
                    qryEnd = stoul(lineVec[0]);
                    qryStart = qryEnd - qrySeqLen + 1;
                    qryLoci = qryEnd + 1;  // Refresh coordinates
                }

                /* ***************************** Refresh coordinates ***************************** */
                idxTmp = -1;
            }
        } else {  // File traversal completed
            // All coordinates return to zero
            // Values changed while traversing the file
            //chromosome number
            chr.clear();
            // The direction, start and end of alignment
            aliRefStrand.clear();
            aliRefStart = 0;
            aliRefEnd = 0;
            aliQryStrand.clear();
            aliQryStart = 0;
            aliQryEnd = 0;
            // Comparison information of the line where the file stream pointer is located
            refStart = 0;
            refSeq.clear();
            refEnd = 0;
            qryStart = 0;
            qrySeq.clear();
            qryEnd = 0;
            // The record traverses to '_ref_seq' and the corresponding coordinates
            refLoci = 0;
            idxTmp = -1;
            qryLoci = 0;

            endBool = true;  // Record whether the file has been traversed
        }
    }

    return 0;
}

/**
 * @brief Find loci
 * 
 * @param chr_      The chromosome you are looking for
 * @param refLoci_  The coordinates you are looking for
 * @return uint32_t  0->Traversed/not found, >0->Coordinates
**/
uint32_t CALNAME::COORTRANS::find_loci(
    string chr_, 
    uint32_t refLoci_
) {
    /* ***************************** After traversing, return 0 ***************************** */
    if (endBool) {
        return 0;
    }
    
    /* ***************************** Not finished traversing ***************************** */
    // Determine whether the chromosome is what you want
    while (chr != chr_ && !endBool) {
        findChrBool = false;  // Record the chromosome that is not required
        next_loci();  // next line
    }
    findChrBool = true;  // Records are required for chromosomes
    // If the traversal is completed, return '0'
    if (endBool) {
        return 0;
    }

    // Determine whether the alignment termination is less than  'refLoci_'
    while (aliRefEnd < refLoci_ && !endBool) {  // We have already judged whether the chromosomes match or not, so there is no need to judge here.
        aliBool = false;  // Recording does not require align
        next_loci();  // next line
    }
    aliBool = true;  // The record needs to be aligned
    // If iteration is complete, return '0'
    if (endBool) {
        return 0;
    }

    // The file pointer points to the desired line
    while (refEnd < refLoci_ && !endBool) {  // We have already judged whether the chromosomes match or not, so there is no need to judge here.
        next_loci();  // next line
    }
    // If iteration is complete, return '0'
    if (endBool) {
        return 0;
    }

    // Find the desired row
    if (refStart <= refLoci_ && refLoci_ <= refEnd && refLoci < refLoci_) {
        for (size_t idx = idxTmp + 1; idx < refSeq.length(); idx++) {
            // ref
            if (refSeq[idx] != '.') {
                if (aliRefStrand == "+") {  // forward
                    refLoci++;
                } else {  // reverse
                    refLoci--;
                }
            }

            // qry
            if (qrySeq[idx] != '.') {
                if (aliQryStrand == "+") {  // forward
                    qryLoci++;
                } else {  // reverse
                    qryLoci--;
                }
            }

            idxTmp = idx;  // refresh index

            // Determine whether it is the coordinates you are looking for
            if (refLoci == refLoci_) {
                break;  // Break out of loop
            }
        }
    } else if (refStart <= refLoci_ && refLoci_ <= refEnd && refLoci == refLoci_) {
        return qryLoci;  // Return coordinates
    } else {  // did not find
        return 0;
    }

    return qryLoci;
}


/**
 * @brief next collinear coordinate
 * 
 * @return 0
**/
int CALNAME::SYRIOUT::next_loci() {
    // Not finished traversing
    if (!endBool) {
        // read line
        string line;
        if (GzChunkReaderClass_->read_line(line)) {
            // Skip empty lines or lines without SYNAL
            if (line.empty() || line.find("SYNAL") == string::npos) {
                return 0;
            }

            // split
            std::istringstream iss(line);
            vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

            // Check if the number of columns is normal
            if (lineVec.size() != 12) {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "'" << fileName_ << "': Expected 12 columns, but got a different number. -> " << line << endl;
                exit(1);
            }

            // Assignment
            chr = lineVec[0];
            refStart = stoul(lineVec[1]);
            refEnd = stoul(lineVec[2]);
            qryStart = stoul(lineVec[6]);
            qryEnd = stoul(lineVec[7]);

            if (refStart > refEnd) std::swap(refStart, refEnd);
            if (qryStart > qryEnd) std::swap(qryStart, qryEnd);
        } else {  // File traversal completed
            // Clear all variables
            chr.clear();
            refStart = 0;
            refEnd = 0;
            qryStart = 0;
            qryEnd = 0;

            endBool = true;  // Record file traversal completed
        }
    }
    return 0;
}

/**
 * @brief Find loci
 * 
 * @param chr_      The chromosome you are looking for
 * @param refLoci_  The coordinates you are looking for
 * 
 * @return int  -1->traversal completed, 0->not found, 1->coordinates
**/
int CALNAME::SYRIOUT::find_loci(
    string chr_, 
    uint32_t refLoci_
) {
    /* ***************************** After traversing, return '-1' ***************************** */
    if (endBool) {
        return -1;
    }
    
    /* ***************************** not traversed ***************************** */
    // Determine whether the chromosome is desired
    while (chr != chr_ && !endBool) {
        next_loci();  // next line
    }
    // If the traversal is complete, return '-1'
    if (endBool) {
        return -1;
    }

    // Determine whether refEnd is less than 'refLoci_'
    while (refEnd < refLoci_ && !endBool) {  // The above has already judged whether the chromosomes match, so there is no need to judge here
        next_loci();  // next line
    }
    // If the traversal is complete, return '-1'
    if (endBool) {
        return -1;
    }

    // in the collinear interval
    if (refStart <= refLoci_ && refLoci_ <= refEnd) {
        return 1;
    } else {  // not in range
        return 0;
    }

    return 0;
}


/**
 * @brief Build a Reference Genome Index
 * 
 * @param referenceFileName -> refgenome
 * 
 * @return refLenMap  map<chr, length>
**/
map<string, uint32_t> CALNAME::build_fasta_index(
    const string & referenceFileName
)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Compute the length of each chromosome in the reference genome." << endl;

    map<string, uint32_t> refLenMap;  // Reference Genome Length Information

    // input file stream
    gzFile gzfp = gzopen(referenceFileName.c_str(), "rb");

    // open file
    if(!gzfp) {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << referenceFileName << "': File not found or maximum open file limit reached. Consider increasing the 'ulimit -n' value to proceed." << endl;
        exit(1);
    } else {
        kseq_t *ks;
        ks = kseq_init(gzfp);
    
        while( kseq_read(ks) > 0 ) {
            string chromosome = ks->name.s;
            uint32_t chrLen = ks->seq.l;
            // string sequence = ks->seq.s;

            refLenMap[chromosome] = chrLen;
        }

        // free memory, close file
        kseq_destroy(ks);
        gzclose(gzfp);
    }

    return refLenMap;
}




/* ************************************** calculate syntenic diversity ************************************** */

/**
 * @brief calculate
 * 
 * @param sampleName                   sample name
 * @param refLenMap                    chromosome length
 * @param SynCoorTmp                   collinear coordinates
 * @param sampleSampleSyriMap          syri.out output path dictionary
 * @param alignsMap                    show-aligns output path dictionary
 * @param insChrStartSampleLociTupMap  map<chr, map<start, map<sample, tuple<start, length> > > >
 * 
 * @return CALSTRUCTURE
**/
CALNAME::CALSTRUCTURE CALNAME::calculate_run(
    const string sampleName, 
    const map<string, uint32_t> & refLenMap, 
    const SYNCOOR & SynCoorTmp, 
    const map<string, map<string, string> > & sampleSampleSyriMap, 
    const map<string, string> & alignsMap, 
    const unordered_map<string, unordered_map<uint32_t, unordered_map<string, tuple<uint32_t, uint32_t> > > > & insChrStartSampleLociTupMap
) {
    // save result
    CALSTRUCTURE CALSTRUCTURETMP;

    // Sample name
    CALSTRUCTURETMP.sampleName = sampleName;

    // The number of sample
    const uint32_t sampleNum = SynCoorTmp.sampleNameVec.size();

    /* ************************************************ Find the index of sampleName ************************************************ */
    // Find the next index after sampleName
    vector<string>::const_iterator findResult = find(SynCoorTmp.sampleNameVec.begin(), SynCoorTmp.sampleNameVec.end(), sampleName);
    int idxTmp = distance(SynCoorTmp.sampleNameVec.begin(), findResult) + 1;

    // If it is the last sample, just skip it
    if (idxTmp > sampleNum - 1) {
        return CALSTRUCTURETMP;
    }
    
    /* ************************************************ Syntenic Coordinate ************************************************ */
    // The map of the collinear coordinate class
    map<string, SYRIOUT> SyriOutMap;  // map<sampleName, SYRIOUT>

    if (sampleName != "reference") {  // It is not necessary to compare the reference genome with others. If it is at the end, skip it
        // Find the sample iterator
        map<string, map<string, string> >::const_iterator findIter1 = sampleSampleSyriMap.find(sampleName);
        if (findIter1 == sampleSampleSyriMap.end()) {  // Warn if no syri.out file for this sample is submitted
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: The sample name '" << sampleName << "' is not found in the sampleSampleSyriMap." << endl;
        }

        // Judging from the post-index sample
        for (size_t i = idxTmp; i < sampleNum; i++) {
            string sampleName2 = SynCoorTmp.sampleNameVec[i];

            if (sampleName2 == "reference") {  // Reference genome without syri.out
                continue;
            }

            // Due to too few coordinates or other reasons, there is no alignment result for this sample, and an empty class is constructed
            if (findIter1 == sampleSampleSyriMap.end()) {
                // SYRIOUT
                SyriOutMap[sampleName2] = SYRIOUT();
                continue;
            }
            
            // Find the iterator of sample2
            map<string, string>::const_iterator findIter2 = findIter1->second.find(sampleName2);

            if (findIter2 == findIter1->second.end()) {  // Warn if no syri.out file for this sample is submitted
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: Neither '" << sampleName << "' nor '" << sampleName2 << "' are present in the sampleSampleSyriMap." << endl;

                // Due to too few coordinates or other reasons, there is no alignment result for this sample, and an empty class is constructed
                SyriOutMap[sampleName2] = SYRIOUT();  // SYRIOUT
                continue;
            }
            string syriOutPath = findIter2->second;

            // SYRIOUT
            SyriOutMap[sampleName2] = SYRIOUT(syriOutPath);
        }
    }

    /* ************************************************ Coordinate transformation ************************************************ */
    // Coordinate conversion class map
    map<string, COORTRANS> CoorTransMap;  // map<sampleName, COORTRANS>

    // Judging from the post-index sample
    for (size_t i = idxTmp - 1; i < sampleNum; i++) {
    
        string sampleName2 = SynCoorTmp.sampleNameVec[i];

        if (sampleName2 == "reference") {  // The reference genome has no .aligns
            continue;
        }

        // Find the aligns path
        map<string, string>::const_iterator findIter1 = alignsMap.find(sampleName2);
        if (findIter1 == alignsMap.end()) {  // If the aligns file corresponding to the sample is not submitted, an error will be reported
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: '" << sampleName2 << "' cannot be found in alignsMap." << endl;
            exit(1);
        }
        
        // Temporary structure
        CoorTransMap[sampleName2] = COORTRANS(findIter1->second);
    }

    /* ************************************************ Calculate syntenic diversity ************************************************ */
    // Record the last coordinate of each chromosome to judge whether it has reached the end of the chromosome
    map<string, uint32_t> refChrMapTmp;  // map<chr, length>

    // Traverse the total coordinate hash table
    for (const auto& iter1 : SynCoorTmp.coorChrLociSampleLociMap) {  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
        // chromosome
        string chrTmp = iter1.first;

        // Coordinates of insertion
        const auto& insChrStartSampleLociTupMapFindIter = insChrStartSampleLociTupMap.find(chrTmp);

        // initialize dictionary
        auto& sampleSynOutMap = CALSTRUCTURETMP.sampleSynOutMap[chrTmp];  // synteny + insertion
        auto& refStartAveSynMap = CALSTRUCTURETMP.sampleSynOutInsMap[chrTmp];  // insertion
        auto& sampleSynOutMapTmp =  CALSTRUCTURETMP.sampleSynOutMapTmp[chrTmp];  // debug

        // record last coordinate
        uint32_t refEndTmp = 0;

        // Record the collinearity score of the last refLoci, if it is the same, it will not be stored. Reduce memory consumption
        uint32_t preSynNum = 0;

        for(const auto& [refStart, SampleLociMap] : iter1.second) {  // map<refStart, map<sample, tuple(start, end)> >
            // start and end of the line
            uint32_t refEnd = get<1>(SampleLociMap.at("reference"));

            // Coordinates of insertion
            const auto& insStartSampleLociTupMapFindIter = 
            (insChrStartSampleLociTupMapFindIter != insChrStartSampleLociTupMap.end()) 
            ? insChrStartSampleLociTupMapFindIter->second.find(refStart) 
            : insChrStartSampleLociTupMapFindIter->second.end();

            // End coordinates for record collinearity
            refChrMapTmp[chrTmp] = refEnd;

            // If there is a distance from the previous coordinate
            if (refStart - refEndTmp > 1) {
                for (uint32_t refLoci = refEndTmp + 1; refLoci < refStart; refLoci++) {
                    uint32_t synNumTmp = 0;  // Collinearity Quantity, Temporary

                    for (uint32_t i = idxTmp; i < sampleNum; i++) {
                        string sampleName2 = SynCoorTmp.sampleNameVec[i];

                        // If the reference has no coordinates with sampleName2, then there are no coordinates
                        if (sampleName == "reference") continue;
                        
                        // if other
                        uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // The coordinates of the reference genome go to sampleName
                        uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // The coordinates of the reference genome go to sampleName2

                        if (lociA == 0 && lociB == 0) { // both del
                            synNumTmp++;  // denoted as syn

                            // debug
                            if (debugCal) {
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                            }
                        } else if (lociA != 0 && lociB == 0) {  // A->synteny, B->deletion
                            continue;
                        } else if (lociA == 0 && lociB != 0) {  // A->deletion, B->synteny
                            continue;
                        } else {  // There are coordinates, and then judge whether A and B of the coordinates are collinear
                            if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0) {  // A and B are collinear
                                synNumTmp++;  // denoted as syn

                                // debug
                                if (debugCal) {
                                    sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                    cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                }
                            }
                        }
                    }

                    // If it is different from the previous synNum, add
                    if (synNumTmp != preSynNum) {
                        sampleSynOutMap[refLoci] = synNumTmp;
                        preSynNum = synNumTmp;
                    }
                }
            } else if (sampleName != "reference" && insStartSampleLociTupMapFindIter != insChrStartSampleLociTupMapFindIter->second.end()) {  // Note: Insertions do not have corresponding coordinates on the reference genome, hence they are stored in sampleSynOutInsMap.
                const auto& insLociMapFindIter = insStartSampleLociTupMapFindIter->second.find(sampleName);

                if (insLociMapFindIter != insStartSampleLociTupMapFindIter->second.end()) {  // The sample is inserted at this site
                    uint32_t sampleStart = get<0>(insLociMapFindIter->second);  // Start of the synteny gap
                    uint32_t sampleLen = get<1>(insLociMapFindIter->second);  // Length of the synteny gap

                    uint32_t synNumInsTmp = 0;  // Amount of synteny

                    for (uint32_t i = idxTmp; i < sampleNum; i++) {
                        const string& sampleName2 = SynCoorTmp.sampleNameVec[i];

                        if (SyriOutMap[sampleName2].find_loci(chrTmp, sampleStart) <= 0) continue;  // If there is no synteny, skip;

                        auto [sampleSynChr, sampleSynStart, sampleSynEnd] = SyriOutMap[sampleName2].get_alignment_loc();

                        synNumInsTmp += sampleSynEnd - sampleSynStart + 1;  // Amount of synteny
                    }

                    // Avoid division by zero and cache the number of samples
                    uint32_t validSampleCount = sampleNum - idxTmp;
                    if (validSampleCount > 0) {
                        float value = synNumInsTmp / static_cast<float>(sampleLen * validSampleCount);
                        refStartAveSynMap[refStart] = (value > 1.0f) ? 1.0f : value;
                    }
                }
            }

            // Calculate the score of the alignmentes area
            for (uint32_t refLoci = refStart; refLoci < refEnd + 1; refLoci++) {
                uint32_t synNumTmp = 0;  // Amount of collinearity, temporary

                for (uint32_t i = idxTmp; i < sampleNum; i++) {
                    string sampleName2 = SynCoorTmp.sampleNameVec[i];
                    
                    // If there are coordinates in 'coor.txt', it means collinearity, record it directly
                    // SampleLociMap  ->  // map<sample, tuple<start, end> >
                    if (get<0>(SampleLociMap.at(sampleName)) > 0 && get<0>(SampleLociMap.at(sampleName2)) > 0) {  // both syntenic
                        synNumTmp++;  // recorded as syn

                        // debug
                        if (debugCal) {
                            sampleSynOutMapTmp[refLoci] += sampleName2 + ";";
                            cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " syn:ture" << endl;
                        }
                    } else if (get<0>(SampleLociMap.at(sampleName)) == 0 && get<0>(SampleLociMap.at(sampleName2)) > 0) {  // one syntenic
                        continue;
                    } else if (get<0>(SampleLociMap.at(sampleName)) > 0 && get<0>(SampleLociMap.at(sampleName2)) == 0) {  // one syntenic
                        continue;
                    } else {  // no syntenic
                        // If the reference has no coordinates with you, then there are no coordinates.
                        if (sampleName == "reference") continue;

                        // other
                        uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // The coordinates of the reference genome go to sampleName
                        uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // The coordinates of the reference genome go to sampleName2

                        if (lociA == 0 && lociB == 0) { // both del
                            synNumTmp++;  // recorded as syn
                            
                            // debug
                            if (debugCal) {
                                sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                            }
                        } else if (lociA != 0 && lociB == 0) {  // One has coordinates, the other does not, one del
                            continue;
                        } else if (lociA == 0 && lociB != 0) {  // One has coordinates, the other does not, one del
                            continue;
                        } else {  // There are coordinates, and then determine whether A and B corresponding to the coordinates are collinear.
                            if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0) {  // A and B are collinear
                                synNumTmp++;  // recorded as syn

                                // debug
                                if (debugCal) {
                                    sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                                    cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                                }
                            }
                        }
                    }
                }

                // If it is different from the previous synNum, add
                if (synNumTmp != preSynNum) {
                    sampleSynOutMap[refLoci] = synNumTmp;
                    preSynNum = synNumTmp;
                }
            }

            // Update coordinates
            refEndTmp = refEnd;
        }
    }

    // Determine whether the end of the chromosome has been reached
    for (const auto& [chrTmp, refPosTmp] : refChrMapTmp) {  // map<chr, length>
        // initialize dictionary
        auto& sampleSynOutMap = CALSTRUCTURETMP.sampleSynOutMap[chrTmp];
        auto& sampleSynOutMapTmp =  CALSTRUCTURETMP.sampleSynOutMapTmp[chrTmp];
        
        // Check whether there is a corresponding chromosome in the reference genome
        map<string, uint32_t>::const_iterator findIter1 = refLenMap.find(chrTmp);
        if (findIter1 == refLenMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The reference genome does not include '" << chrTmp << "'." << endl;
            exit(1);
        }

        uint32_t refLen = findIter1->second;  // chromosome length

        // If the end is reached, exit the loop
        if (refPosTmp >= refLen) break;
        
        // If the end of the chromosome is not reached
        uint32_t preSynNum = 0;
        for (uint32_t refLoci = refPosTmp + 1; refLoci < refLen + 1; refLoci++) {
            uint32_t synNumTmp = 0;  // Amount of collinearity, temporary

            for (size_t i = idxTmp; i < sampleNum; i++) {
                string sampleName2 = SynCoorTmp.sampleNameVec[i];

                // If the reference has no coordinates with you, then there are no coordinates.
                if (sampleName == "reference") continue;
                
                // other
                uint32_t lociA = CoorTransMap[sampleName].find_loci(chrTmp, refLoci);  // The coordinates of the reference genome go to sampleName
                uint32_t lociB = CoorTransMap[sampleName2].find_loci(chrTmp, refLoci);  // The coordinates of the reference genome go to sampleName2

                if (lociA == 0 && lociB == 0) {  // both del
                    synNumTmp++;  // marked as syn

                    // debug
                    if (debugCal) {
                        sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                        cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                    }
                } else if (lociA != 0 && lociB == 0) {  // One has coordinates, the other does not, one del
                    continue;
                } else if (lociA == 0 && lociB != 0) {  // One has coordinates, the other does not, one del
                    continue;
                } else {  // There are coordinates, and then determine whether A and B of the coordinates are collinear.
                    if (SyriOutMap[sampleName2].find_loci(chrTmp, lociA) > 0) {  // A and B are collinear
                        synNumTmp++;  // recorded as syn

                        // debug
                        if (debugCal) {
                            sampleSynOutMapTmp[refLoci] += sampleName2 + ":" + to_string(lociA) + "-" + to_string(lociB) + ";";
                            cerr << chrTmp << " " << sampleName << " " << sampleName2 << " " << refLoci << " lociA:" << lociA << " lociB:" << lociB << " syn:ture" << endl;
                        }
                    }
                }
            }

            // If it is different from the previous synNum, add
            if (synNumTmp != preSynNum) {
                sampleSynOutMap[refLoci] = synNumTmp;
                preSynNum = synNumTmp;
            }
        }
    }

    return CALSTRUCTURETMP;
}




/* ********************************************* memory-saving mode ********************************************* */

/**
 * @brief merge result
 * 
 * @param calOutStr                 All calculation results for a sample
 * @param refLenMap                 Chromosome length information, map<string, length>
 * @param chrSynNumVecMap           Save the final result map<chr, vector<synNum> >
 * @param chrLociAveSynNumInsMap    Save the final result map<chr, map<refStart, synNumAve> > (insertion)
 * 
 * @return 0
**/
int CALNAME::merge(
    const CALSTRUCTURE& calOutStr,
    const map<string, uint32_t> & refLenMap,
    map<string, vector<uint32_t> >& chrSynNumVecMap,
    map<string, map<uint32_t, float> > & chrLociAveSynNumInsMap
) {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Combine the computation results of " << calOutStr.sampleName << "." << endl;

    // synteny + deletion
    for(const auto& [chromosome, chrLen] : refLenMap) {  // map<string, length>
        // Save the final result map<chr, vector<synNum> >
        auto& SynNumVec = chrSynNumVecMap[chromosome];

        // If it is merged for the first time, first initialize the length of SynNumVec to the length of the chromosome.
        if (SynNumVec.empty()) {
            SynNumVec = vector<uint32_t>(chrLen, 0);
        }
        
        // Iterator corresponding to chromosome
        map<string, map<uint32_t, uint32_t> >::const_iterator iter0 = calOutStr.sampleSynOutMap.find(chromosome);
        map<string, map<uint32_t, string> >::const_iterator iter3 = calOutStr.sampleSynOutMapTmp.find(chromosome);

        // Fix issue #8: If the chromosome does not exist or is empty, skip it directly
        if (iter0 == calOutStr.sampleSynOutMap.end() || iter0->second.empty()) {
            continue;
        }

        // Temporary iterator that records the scores of the current and next positions
        map<uint32_t, uint32_t>::const_iterator iter1 = iter0->second.begin();
        map<uint32_t, uint32_t>::const_iterator iter2 = std::next(iter1);
        
        // Record the synNum of the previous node
        uint32_t preSynNumTmp = 0;

        // Circular chromosome length, each position added separately
        for (uint32_t refLoci = 1; refLoci < chrLen + 1; refLoci++) {
            // Determine whether the current position points to the first iterator
            if (iter0 != calOutStr.sampleSynOutMap.end() && refLoci == iter1->first) {  // If so, assign the number corresponding to the iterator
                preSynNumTmp = iter1->second;
            } else if (iter0 != calOutStr.sampleSynOutMap.end() && iter2 != iter0->second.end() && refLoci == iter2->first) {  // Go to the next node and update the iterator
                iter1++;  // update iterator
                iter2++;  // update iterator

                preSynNumTmp = iter1->second;  // Update coordinates
            }

            // site frequency overlay
            SynNumVec[refLoci - 1] += preSynNumTmp;

            // debug
            if (debugCal) {
                if (preSynNumTmp > 0) {
                    cerr << chromosome << " " << refLoci << " " << calOutStr.sampleName << " " << preSynNumTmp << " " << iter3->second.at(refLoci) << endl;
                }
            }
        }
    }

    // deletion
    for (auto& refStartAveSynMapIter : calOutStr.sampleSynOutInsMap) {  // map<chr, map<refLoci, synNumAve> >
        for (auto & StartAveSynIter : refStartAveSynMapIter.second) {  // map<refLoci, synNumAve>
            auto& LociAveSynNumInsMap = chrLociAveSynNumInsMap[refStartAveSynMapIter.first];
            auto findIter = LociAveSynNumInsMap.find(StartAveSynIter.first);
            if (findIter == LociAveSynNumInsMap.end()) {
                LociAveSynNumInsMap[StartAveSynIter.first] = StartAveSynIter.second;
            } else {
                findIter->second += StartAveSynIter.second;
            }
        }
    }

    return 0;
}


/**
 * @brief calculate
 * 
 * @param refLenMap                    Reference genome length
 * @param SynCoorTmp                   collinear coordinates
 * @param sampleSampleSyriMap          syri.out output path dictionary
 * @param alignsMap                    show-aligns output path dictionary
 * @param insChrStartSampleLociTupMap  map<chr, map<start, map<sample, tuple<start, length> > > >
 * @param allSynNum                    All combination quantities
 * @param outputFileName               output file name
 * 
 * @return 0
**/
int CALNAME::calculate(
    const map<string, uint32_t> & refLenMap, 
    const SYNCOOR & SynCoorTmp, 
    const map<string, map<string, string> > & sampleSampleSyriMap, 
    const map<string, string> & alignsMap, 
    const unordered_map<string, unordered_map<uint32_t, unordered_map<string, tuple<uint32_t, uint32_t> > > > & insChrStartSampleLociTupMap, 
    const uint32_t & allSynNum, 
    const string & outputFileName
) {
    // final result
    map<string, vector<uint32_t> > chrSynNumVecMap;  // map<chr, vector<synNum> >
    map<string, map<uint32_t, float> > chrLociAveSynNumInsMap;  // map<refChr, map<refStart, aveSynSum> >

    /* ********************************************** calculate syntenic diversity ********************************************** */
    // process pool
    ThreadPool pool(threadsCal);
    const int MAX_THREADS_NUM = threadsCal*2;  // The most stored element in a multi-threaded vector, write out when it is greater than

    // init
    pool.init();

    // Save the results of multiple threads
    /*
        map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >
        map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
    */
    vector<future<CALSTRUCTURE> > calOutStrFutureVec;
    
    // sample name
    for (const auto& iter1 : SynCoorTmp.sampleNameVec) {  // vector<sampleName>
        string sampleNameTmp = iter1;

        // last sample skipped
        if (sampleNameTmp == SynCoorTmp.sampleNameVec.back()) continue;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate: " << sampleNameTmp << endl;

        // Submit and save results in multiple threads
        calOutStrFutureVec.push_back(
            pool.submit(
                calculate_run, 
                sampleNameTmp, 
                ref(refLenMap), 
                ref(SynCoorTmp), 
                ref(sampleSampleSyriMap), 
                ref(alignsMap), 
                ref(insChrStartSampleLociTupMap)
            )
        );

        // If the length of calOutStrVecTmp is greater than the threshold, write it first and then submit the task
        if (calOutStrFutureVec.size() >= MAX_THREADS_NUM) {
            // Traversing the results returned by multiple threads
            for (auto&& calOutStrFuture : calOutStrFutureVec) {  // vector<future<CALSTRUCTURE> >
                CALSTRUCTURE calOutStr = move(calOutStrFuture.get());  // future<CALSTRUCTURE>
                // merged result
                merge(calOutStr, refLenMap, chrSynNumVecMap, chrLociAveSynNumInsMap);
            }

            calOutStrFutureVec.clear();
            vector<future<CALSTRUCTURE> >().swap(calOutStrFutureVec);

            malloc_trim(0); // 0 is for heap memory
        }
    }

    // The last combined results
    if (calOutStrFutureVec.size() > 0) {
        // Traversing the results returned by multiple threads
        for (auto&& calOutStrFuture : calOutStrFutureVec) {  // vector<future<CALSTRUCTURE> >
            CALSTRUCTURE calOutStr = move(calOutStrFuture.get());  // future<CALSTRUCTURE>
            // merged result
            merge(calOutStr, refLenMap, chrSynNumVecMap, chrLociAveSynNumInsMap);
        }
        calOutStrFutureVec.clear();
        vector<future<CALSTRUCTURE> >().swap(calOutStrFutureVec);

        malloc_trim(0); // 0 is for heap memory
    }

    /* ********************************************** save the result (synteny + deletion) ********************************************** */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;

    SAVE SAVEClass(outputFileName);

    stringstream outStream;  // Use stringstream instead of string concatenation
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);
    outStream << "#CHROM\tPOS\tAll_States\tSyntenic_States\tSyntenic_Diversity\n";

    for (const auto& [chromosome, SynNumVec] : chrSynNumVecMap) {  // map<chr, vector<synNum> >
        uint32_t loci = 1;
        for (const auto& synNum : SynNumVec) {  // vector<synNum>
            outStream << chromosome << "\t" << loci << "\t" << allSynNum << "\t" << synNum << "\t" << 1.0 - synNum/(double)allSynNum << "\n";

            if (outStream.tellp() >= CACHE_SIZE) {  // Cache size is 10mb
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // empty stringstream
                outStream.str(string());
                outStream.clear();
            }
            ++loci;
        }

        if (outStream.tellp() > 0) {  // write one last time
            string outTxt = outStream.str();
            SAVEClass.save(outTxt);
            // empty stringstream
            outStream.str(string());
            outStream.clear();
        }
    }

    /* ********************************************** save the result (insertion) ********************************************** */
    string outputFileNameIns = outputFileName;
    if (!outputFileNameIns.empty()) outputFileNameIns += ".ins";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results (insertion) are being saved to '" << outputFileNameIns << "'" << endl;
    SAVE SAVEInsClass(outputFileNameIns);
    // total number of samples
    const uint32_t sampleNum = SynCoorTmp.sampleNameVec.size();
    // Score
    outStream << "#Insertion\n#CHROM\tPOS\tSyntenic_Diversity\n";
    for (auto LociAveSynNumInsMapIter : chrLociAveSynNumInsMap) {  // map<refChr, map<refStart, aveSynSum> >
        for (auto AveSynNumInsMapIter : LociAveSynNumInsMapIter.second) {  // map<refStart, aveSynSum>
            outStream << LociAveSynNumInsMapIter.first << "\t" << AveSynNumInsMapIter.first << "\t" << 1.0 - AveSynNumInsMapIter.second/(float)sampleNum << "\n";

            if (outStream.tellp() >= CACHE_SIZE) {  // Cache size is 10mb
                string outTxt = outStream.str();
                SAVEInsClass.save(outTxt);
                // empty stringstream
                outStream.str(string());
                outStream.clear();
            }
        }
    }

    if (outStream.tellp() > 0) {  // write one last time
        string outTxt = outStream.str();
        SAVEInsClass.save(outTxt);
        // empty stringstream
        outStream.str(string());
        outStream.clear();
    }
    

    // Close the thread pool
    pool.shutdown();

    return 0;
}



/* ********************************************* quick mode ********************************************* */

/**
 * @brief merge result
 * 
 * @param chromosome                         chromosome
 * @param chrLen                             chromosome length
 * @param CALSTRUCTUREVec                    output of threads
 * 
 * @return tuple<chromosome, synOutVecTmp, chrLociAveSynNumInsMap>   tuple<chr, vector<synNum>, map<refStart, synNumAve> >
**/
tuple<string, vector<uint32_t>, map<uint32_t, float> > CALNAME::merge_fast(
    string chromosome, 
    uint32_t chrLen, 
    const vector<CALSTRUCTURE> & CALSTRUCTUREVec
) {
    // save result
    vector<uint32_t> synOutVecTmp(chrLen, 0);  //  vector<synNum>

    map<uint32_t, float> LociAveSynNumInsMap;  // Save the final result map<refStart, synNumAve> (insertion)
    
    // synteny and deletion
    for (const auto& CALSTRUCTURETmp : CALSTRUCTUREVec) {  // vector<CALSTRUCTURE>
        // Iterator corresponding to chromosome
        map<string, map<uint32_t, uint32_t> >::const_iterator iter0 = CALSTRUCTURETmp.sampleSynOutMap.find(chromosome);
        map<string, map<uint32_t, string> >::const_iterator iter3 = CALSTRUCTURETmp.sampleSynOutMapTmp.find(chromosome);

        // Fix issue #8: If the chromosome does not exist or is empty, skip it directly
        if (iter0 == CALSTRUCTURETmp.sampleSynOutMap.end() || iter0->second.empty()) {
            continue;
        }

        // Temporary iterator that records the scores of the current and next positions
        map<uint32_t, uint32_t>::const_iterator iter1 = iter0->second.begin();
        map<uint32_t, uint32_t>::const_iterator iter2 = std::next(iter1);
        
        // Record the synNum of the previous node
        uint32_t preSynNumTmp = 0;

        // Circular chromosome length, each position added separately
        for (uint32_t refLoci = 1; refLoci < chrLen + 1; refLoci++) {
            // Determine whether the current position points to the first iterator
            if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end() && refLoci == iter1->first) {  // If so, assign the number corresponding to the iterator
                preSynNumTmp = iter1->second;
            } else if (iter0 != CALSTRUCTURETmp.sampleSynOutMap.end() && iter2 != iter0->second.end() && refLoci == iter2->first) {  // Go to the next node and update the iterator
                iter1++;  // update iterator
                iter2++;  // update iterator

                preSynNumTmp = iter1->second;  // Update coordinates
            }

            // site frequency overlay
            synOutVecTmp[refLoci - 1] += preSynNumTmp;

            // debug
            if (debugCal) {
                if (preSynNumTmp > 0)
                {
                    cerr << chromosome << " " << refLoci << " " << CALSTRUCTURETmp.sampleName << " " << preSynNumTmp << " " << iter3->second.at(refLoci) << endl;
                }
            }
        }
    }

    // deletion
    for (const auto& calOutStr : CALSTRUCTUREVec) {  // vector<calOutStr>
        for (auto& refStartAveSynMapIter : calOutStr.sampleSynOutInsMap) {  // map<chr, map<refLoci, synNumAve> >
            if (refStartAveSynMapIter.first != chromosome) continue;  // Skip if chromosome doesn't match
            for (auto & StartAveSynIter : refStartAveSynMapIter.second) {  // map<refLoci, synNumAve>
                auto findIter = LociAveSynNumInsMap.find(StartAveSynIter.first);
                if (findIter == LociAveSynNumInsMap.end()) {
                    LociAveSynNumInsMap[StartAveSynIter.first] = StartAveSynIter.second;
                } else {
                    findIter->second += StartAveSynIter.second;
                }
            }
        }
    }

    return make_tuple(chromosome, synOutVecTmp, LociAveSynNumInsMap);
}


/**
 * @brief calculate
 * 
 * @param refLenMap                    chromosome length
 * @param SynCoorTmp                   collinear coordinates
 * @param sampleSampleSyriMap          syri.out output path dictionary
 * @param alignsMap                    show-aligns output path dictionary
 * @param insChrStartSampleLociTupMap  map<chr, map<start, map<sample, tuple<start, length> > > >
 * @param allSynNum                    All combination quantities
 * @param outputFileName               Output file name
 * 
 * @return 0
**/
int CALNAME::calculate_fast(
    const map<string, uint32_t> & refLenMap, 
    const SYNCOOR & SynCoorTmp, 
    const map<string, map<string, string> > & sampleSampleSyriMap, 
    const map<string, string> & alignsMap, 
    const unordered_map<string, unordered_map<uint32_t, unordered_map<string, tuple<uint32_t, uint32_t> > > > & insChrStartSampleLociTupMap, 
    const uint32_t & allSynNum, 
    const string & outputFileName
) {
    /* ********************************************** calculate syntenic diversity ********************************************** */
    // thread pool
    ThreadPool pool(threadsCal);

    // init
    pool.init();

    // Save multi-threaded results
    /*
        map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >
        map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
    */
    vector<future<CALSTRUCTURE> > calOutStrVecTmp;
    
    // sample name
    for (const auto& iter1 : SynCoorTmp.sampleNameVec) {  // vector<sampleName>
        string sampleNameTmp = iter1;

        // last sample skipped
        if (sampleNameTmp == SynCoorTmp.sampleNameVec[SynCoorTmp.sampleNameVec.size() - 1]) {
            continue;
        }

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate: " << sampleNameTmp << endl;

        // Submit and save results through multiple threads
        calOutStrVecTmp.push_back(
            pool.submit(
                calculate_run, 
                sampleNameTmp, 
                ref(refLenMap), 
                ref(SynCoorTmp), 
                ref(sampleSampleSyriMap), 
                ref(alignsMap), 
                ref(insChrStartSampleLociTupMap)
            )
        );

        sleep(0.1);  // Thread interval
    }

    // Multi-threaded result saving
    vector<CALSTRUCTURE> CALSTRUCTUREVecTmp;
    for (size_t i = 0; i < calOutStrVecTmp.size(); i++) {
        CALSTRUCTUREVecTmp.push_back(move(calOutStrVecTmp[i].get()));
    }
    calOutStrVecTmp.clear();
    vector<future<CALSTRUCTURE> >().swap(calOutStrVecTmp);
    malloc_trim(0);	// 0 is for heap memory


    /* ********************************************** merge the result ********************************************** */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Merge the result." << endl;
    vector<future<tuple<string, vector<uint32_t>, map<uint32_t, float> > > > mergeOutVec;
    for (const auto& iter1 : refLenMap) {  // map<chromosome, length>
        string chromosome = iter1.first;
        uint32_t chrLen = iter1.second;

        // Submit and save results through multiple threads
        mergeOutVec.push_back(
            pool.submit(
                merge_fast, 
                chromosome, 
                chrLen, 
                ref(CALSTRUCTUREVecTmp)
            )
        );

        sleep(0.1);  // Thread interval
    }
    

    /* ********************************************** save the result ********************************************** */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;
    string outputFileNameIns = outputFileName;
    if (!outputFileNameIns.empty()) outputFileNameIns += ".ins";
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results (insertion) are being saved to '" << outputFileNameIns << "'" << endl;

    SAVE SAVEClass(outputFileName);
    SAVE SAVEInsClass(outputFileNameIns);

    stringstream outStream; // Use stringstream instead of string concatenation
    stringstream outStreamIns; // Use stringstream instead of string concatenation
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);
    outStreamIns.str().reserve(CACHE_SIZE);
    // total number of samples
    const uint32_t sampleNum = SynCoorTmp.sampleNameVec.size();
    outStream << "#CHROM\tPOS\tAll_States\tSyntenic_States\tSyntenic_Diversity\n";
    outStreamIns << "#Insertion\n#CHROM\tPOS\tSyntenic_Diversity\n";

    for (size_t i = 0; i < mergeOutVec.size(); ++i) {  // vector<future<tuple<string, vector<uint32_t> > > > 
        auto [chromosome, calOutVec, LociAveSynNumInsMap] = move(mergeOutVec[i].get());  // Get multi-threaded results

        // synteny and deletion
        uint32_t loci = 1;
        for (const auto& calOut : calOutVec) {  // // vector<uint32_t>
            outStream << chromosome << "\t" << loci << "\t" << allSynNum << "\t" << calOut << "\t" << 1.0 - calOut/(double)allSynNum << "\n";

            if (outStream.tellp() >= CACHE_SIZE) {  // Cache size is 10mb
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // Clear stringstream
                outStream.str(string());
                outStream.clear();
            }
            ++loci;
        }

        if (outStream.tellp() > 0) {  // Write for the last time
            string outTxt = outStream.str();
            SAVEClass.save(outTxt);
            // Clear stringstream
            outStream.str(string());
            outStream.clear();
        }

        // insertion
        for (auto AveSynNumInsMapIter : LociAveSynNumInsMap) {  // map<refStart, aveSynSum> (insertion)
            outStreamIns << chromosome << "\t" << AveSynNumInsMapIter.first << "\t" << 1.0 - AveSynNumInsMapIter.second/(float)sampleNum << "\n";

            if (outStreamIns.tellp() >= CACHE_SIZE) {  // Cache size is 10mb
                string outTxt = outStreamIns.str();
                SAVEInsClass.save(outTxt);
                // empty stringstream
                outStreamIns.str(string());
                outStreamIns.clear();
            }
        }

        if (outStreamIns.tellp() > 0) {  // write one last time
            string outTxt = outStreamIns.str();
            SAVEInsClass.save(outTxt);
            // empty stringstream
            outStreamIns.str(string());
            outStreamIns.clear();
        }
    }

    // Close thread pool
    pool.shutdown();

    return 0;
}