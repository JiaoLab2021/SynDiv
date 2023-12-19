// g++ -c src/coor.cpp -O3 -std=c++17
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <future>
#include <malloc.h>
#include "../include/coor.hpp"
#include "../include/ThreadPool.hpp"

using namespace std;

// Global variable
int thresholdLength = 1000;  // Collinear coordinates ref and qry length ratio threshold

// Whether to debug code
bool debugCoor = false;


int main_coor(int argc, char* argv[])
{
    // loc Output file
    string synLocFileName;

    // Enter a list of files and a name  show-aligns
    vector<string> inputFiles;
    vector<string> inputTitles;
    // Determine whether there is titlename
    bool haveTitles = false;

    // Output file name
    string outputFileName;

    // threads number
    int threads = 30;

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
                debugCoor = true;
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


    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;

    
    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index ..." << endl;
    map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > > synLocSampleVecMap;  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    synLocSampleVecMap = COOR::build_syn_idx(
        synLocFileName
    );


    /* ************************************ Find Syntenic Coordinates on query genomes ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find Syntenic Coordinates on query genomes ..." << endl;
    ThreadPool pool(threads);  // Thread pool

    // Initializes the thread pool
    pool.init();

    // Save the result of multiple threads
    vector<future<COOR::synAllStructure> > chrStartSynQryLocMapVec;

    // Obtain the coordinates of collinearity on qry
    uint32_t indexTmp = 0;  // Record the index of the sample name
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
        // Multithreading submits and saves results
        chrStartSynQryLocMapVec.push_back(
            pool.submit(
                COOR::get_syn_coor, 
                sampleName, 
                it, 
                ref(synLocSampleVecMap), 
                false
            )
        );

        indexTmp++;  // Index stack
    }

    // Multithreaded result saving
    map<string, unordered_map<string, unordered_map<uint32_t, tuple<uint32_t, uint32_t> > > > sampleChrStartSynQryLocMap;
    for (size_t i = 0; i < chrStartSynQryLocMapVec.size(); i++)
    {
        COOR::synAllStructure synAllStructureTmp = move(chrStartSynQryLocMapVec[i].get());
        sampleChrStartSynQryLocMap[move(synAllStructureTmp.sampleName)] = move(synAllStructureTmp.chrStartSynQryLocMap);
    }
    chrStartSynQryLocMapVec.clear();
    vector<future<COOR::synAllStructure> >().swap(chrStartSynQryLocMapVec);
    malloc_trim(0);	// 0 is for heap memory

    // Close thread pool
    pool.shutdown();
    
    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;
    COOR::save_result(
        synLocSampleVecMap, 
        sampleChrStartSynQryLocMap, 
        outputFileName
    );

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ..." << endl;

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



/**
 * @brief Parse the '-- BEGIN alignment 'field
 * 
 * @param informationTmp    '+1 278 - 1703'
 * 
 * @return tuple<strand, start, end>
*/
tuple<string, uint32_t, uint32_t> COOR::get_alignment_loc(string informationTmp)
{
    // Get alignment direction, start and end information
    istringstream iss(informationTmp);  // +1 278 - 1703
    string strandTmp;
    uint32_t startTmp = 0;
    uint32_t endTmp = 0;
    char sep;

    iss >> strandTmp >> startTmp >> sep >> endTmp;

    // For reverse alignment, switch the start and end coordinates
    if (strandTmp == "-1") {
        swap(startTmp, endTmp);
    }

    return make_tuple(strandTmp, startTmp, endTmp);
}


/**
 * @brief Find the coordinates of collinearity on qry
 * 
 * @param synLoc        Collinear coordinates
 * @param refStart      Start of the ref
 * @param refEnd        End of ref
 * @param refSeq        ref sequence
 * @param qryStrand     qry Compare directions
 * @param qryStart      Indicates the start of qry
 * @param qryEnd        End of qry
 * @param qrySeq        qry
 * 
 * @return qrySynLoc coordinates of syn on qry, 0- Not found
*/
uint32_t COOR::find_qry_syn(
    const uint32_t & synLoc, 
    uint32_t refStart, 
    uint32_t refEnd, 
    const string & refSeq, 
    const string & AliQryStrand, 
    uint32_t qryStart, 
    uint32_t qryEnd, 
    const string & qrySeq
)
{
    uint32_t qrySynLoc = 0;

    if (refStart <= synLoc && synLoc <= refEnd)  // If you include collinear coordinates
    {
        uint32_t iIdx = 0;  // Record the location corresponding to synStart

        // Subtract 1 from the coordinates
        if (refStart > 0) {
            --refStart;
        }

        for (size_t i = 0; i < refSeq.size(); ++i)  // Cyclic ref sequence
        {
            if (refSeq[i] != '.')
            {
                ++refStart;  // ref coordinate increment
            }

            if (refStart == synLoc)  // Find the index of syn coordinates
            {
                iIdx = i;
                break;  // Exit the loop when the coordinates are found
            }
        }

        if (AliQryStrand == "+1")  // qry is for forward comparison
        {
            // Subtract 1 from the coordinates
            if (qryStart > 0) {
                --qryStart;
            }

            for (size_t i = 0; i < qrySeq.size(); ++i)
            {
                if (qrySeq[i] != '.')
                {
                    ++qryStart;  // ref coordinate increment
                }
                
                if (i == iIdx)  // Loop to this position and exit
                {
                    break;  // Exit the loop when the coordinates are found
                }
            }
            qrySynLoc = qryStart;
        }
        else  // qry is the reverse comparison
        {
            ++qryEnd;  // Add 1 to the coordinates
            for (size_t i = 0; i < qrySeq.size(); ++i)
            {
                if (qrySeq[i] != '.')
                {
                    // qry coordinate decrement
                    if (qryEnd > 0) {
                        --qryEnd;
                    }
                }
                
                if (i == iIdx)  // Loop to this position and exit
                {
                    break;
                }
            }
            qrySynLoc = qryEnd;
        }
    }

    return qrySynLoc;
}


/**
 * @brief Find the coordinates of collinearity on qry
 *
 * @param chrStartSynQryLocMap            synAllStructure
 * @param refChr                          Chromosome number
 * @param qrySynStartTmp                  find_qry_syn found qrySynStart
 * @param qrySynEndTmp                    find_qry_syn found qrySynEnd
 * @param refStart                        Start of the line ref
 * @param refEnd                          Terminates the ref of the alignment line
 * @param qryStart                        Start of the row qry
 * @param qryEnd                          Terminates the qry of the tape comparison row
 * @param synStart                        The ref start of the syn you are looking for
 * @param synEnd                          ref of syn currently sought terminates
 * @param qrySynStart                     Indicates the qry start of syn
 * @param qrySynEnd                       The qry of syn is terminated
 * @param whileBool                       determines if a while loop is still needed
 * @param aliRowStartNum                  The approximate number of lines in which the syn starts
 * @param aliRowEndNum                    Approximate number of lines in the alignment where syn terminates
 *
 * @return qrySynLoc coordinates of syn on qry, 0- not found
**/
uint32_t COOR::syn_all_loc_push(
    synAllStructure & chrStartSynQryLocMap, 
    const string & refChr, 
    const uint32_t & qrySynStartTmp, 
    const uint32_t & qrySynEndTmp, 
    const uint32_t & refStart, 
    const uint32_t & refEnd, 
    const uint32_t & qryStart, 
    const uint32_t & qryEnd, 
    const uint32_t & synStart, 
    const uint32_t & synEnd, 
    uint32_t & qrySynStart, 
    uint32_t & qrySynEnd, 
    bool & whileBool, 
    int32_t & aliRowStartNum, 
    int32_t & aliRowEndNum
)
{
    whileBool = false;

    if (qrySynStartTmp != 0) {
        if (qrySynStart > 0) {  // Having previously found qrySynStart, select the sequence closest to preQrySynStart
            if (abs(static_cast<int64_t>(qrySynStart) - static_cast<int64_t>(synStart)) > abs(static_cast<int64_t>(qrySynStartTmp) - static_cast<int64_t>(synStart)))
            {
                qrySynStart = qrySynStartTmp;
            } else if (debugCoor) {
                cerr << "skip_start:" << qrySynStartTmp << endl;
            }
        }
        else  // If you don't find it, assign it
        {
            qrySynStart = qrySynStartTmp;
        }
        
        aliRowStartNum = aliRowEndNum;  // Reset aliRowStartNum and jump to aliRowEndNum to find the end coordinates
    }
    if (qrySynEndTmp != 0 && qrySynStart != 0)  // The beginning must be found first
    {
        if (qrySynEnd > 0)  // Having previously found qrySynEnd, select the sequence closest to preQrySynEnd
        {
            if (abs(static_cast<int64_t>(qrySynEnd) - static_cast<int64_t>(synEnd)) > abs(static_cast<int64_t>(qrySynEndTmp) - static_cast<int64_t>(synEnd)))
            {
                qrySynEnd = qrySynEndTmp;
            }
            else if (debugCoor)
            {
                cerr << "skip_end:" << qrySynEndTmp << endl;
            }
        }
        else  // If you don't find it, assign it
        {
            qrySynEnd = qrySynEndTmp;
        }

        // If it is greater than 0, then assign a value, and if the difference multiple of the comparison length is smaller than the thresholdLength, then assign a value
        uint32_t refSynLen = (synEnd > synStart) ? synEnd - synStart : synStart - synEnd;
        uint32_t qrySynLen = (qrySynEnd > qrySynStart) ? qrySynEnd - qrySynStart : qrySynStart - qrySynEnd;
        if (qrySynEnd > 0 && max(refSynLen, qrySynLen)/(float)min(refSynLen, qrySynLen) < thresholdLength) {
            chrStartSynQryLocMap.chrStartSynQryLocMap[refChr][synStart] = make_tuple(
                min(qrySynStart, qrySynEnd), 
                max(qrySynStart, qrySynEnd)
            );
        }
    }

    // debug
    if (qrySynStartTmp != 0 || qrySynEndTmp != 0) {
        if (debugCoor) {
            cerr << "ref:" << refChr << "-" << refStart << "-" << refEnd << "\t" << 
                "qry:" << qryStart << "-" << qryEnd << "\t" << 
                "refSyn:" << synStart << "-" << synEnd << "\t" <<
                "qrySyn:" << qrySynStart << "-" << qrySynEnd << "\t" <<
                "synOut:" << qrySynStartTmp << "-" << qrySynEndTmp << endl;
        }
    }

    return 0;
}


/**
 * @brief Update temporary syn coordinates
 * 
 * @param sampleName            Sample name, used to determine whether the line is to be judged
 * @param synLocSampleVecMap    map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param synChr                chromosome
 * @param synIdx                syn is indexed in the vector of the chromosome
 * @param synStart              Start position of syn
 * @param synEnd                End position of syn
 * @param qrySynStart           Start position of syn on qry
 * @param qrySynEnd             End position of syn on qry
 * @param whileBool             Determine whether to jump out of the while loop
 * 
 * @return 0
*/
int COOR::renew_syn_loc(
    const string& sampleName, 
    const map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > >& synLocSampleVecMap, 
    const string& synChr, 
    uint32_t& synIdx, 
    uint32_t& synStart, 
    uint32_t& synEnd, 
    uint32_t& qrySynStart, 
    uint32_t& qrySynEnd, 
    bool& whileBool
)
{
    map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > >::const_iterator findIter = synLocSampleVecMap.find(synChr);  //  Look for chromosomes in a collinear map
    if (findIter != synLocSampleVecMap.end()) {  // Got it.
        auto& synLocSampleVec = findIter->second;  // vector<tuple<uint32_t, uint32_t, vector<string> > >

        if (synIdx < synLocSampleVec.size()) {  // Prevention of transgression
            // If the collinearity does not have a corresponding sample name, proceed to the next coordinate
            while (find(get<2>(synLocSampleVec[synIdx]).begin(), get<2>(synLocSampleVec[synIdx]).end(), sampleName) == get<2>(synLocSampleVec[synIdx]).end()) {
                synIdx++;  // Update coordinate
                // If you cross the line, zero the coordinates and exit the function
                if (synIdx >= synLocSampleVec.size()) {
                    synIdx = 0;
                    synStart = 0;
                    synEnd = 0;
                    qrySynStart = 0;
                    qrySynEnd = 0;
                    whileBool = false;

                    return 0;
                }
            }
            
            synStart = get<0>(synLocSampleVec[synIdx]);
            synEnd = get<1>(synLocSampleVec[synIdx]);
            synIdx++;  // Update coordinate
            qrySynStart = 0;
            qrySynEnd = 0;
            whileBool = true;
        } else {
            synIdx = 0;
            synStart = 0;
            synEnd = 0;
            qrySynStart = 0;
            qrySynEnd = 0;
            whileBool = false;
        }
    }
    
    return 0;
}


/**
 * @brief The syn coordinate index common to all samples in the component
 * 
 * @param inputFileName         syn_loc output file
 * 
 * @return synLocSampleVecMap   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
*/
map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > > COOR::build_syn_idx(const string & inputFileName)
{
    // Save collinear coordinates
    map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > > synLocSampleVecMap;  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >

    // open syntenic coordinations file
    GzChunkReader GzChunkReaderClass(inputFileName);

    // read line
    string line;
    while (GzChunkReaderClass.read_line(line)) {
        // Skip empty line
        if (line.empty()) {
            continue;
        }

        // split
        std::istringstream iss(line);
        vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

        synLocSampleVecMap[lineVec[0]].push_back(make_tuple(stoul(lineVec[1]), stoul(lineVec[2]), split(lineVec[4], ",")));  // Store the coordinates of collinearity on the ref and the sample name for that site
    }

    return synLocSampleVecMap;
}


/**
 * @brief Get the coordinates of syn on qry
 * 
 * @param sampleName                 sample name
 * @param inputFileName              show-aligns output file
 * @param synLocSampleVecMap         map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param findRevBool                Whether to traverse the reverse alignment sequence
 * 
 * @return chrStartSynQryLocMap      synAllStructure
*/
COOR::synAllStructure COOR::get_syn_coor(
    string sampleName, 
    const string& inputFileName, 
    const map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > >& synLocSampleVecMap, 
    const bool& findRevBool
)
{
    // Obtain the coordinates of collinearity on qry
    // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
    synAllStructure chrStartSynQryLocMap;

    // If no sample name is submitted, extract according to the file
    if (sampleName.empty()) {
        vector<string> inputFileNameVecTmp = split(inputFileName, "/");  // Path splitting
        sampleName = inputFileNameVecTmp.back();

        std::regex reg1;
        if (sampleName.substr(sampleName.length() - 3) == ".gz") {
            reg1 = std::regex(".aligns.gz");
        } else if (sampleName.substr(sampleName.length() - 3) == ".GZ") {
            reg1 = std::regex(".aligns.GZ");
        } else {
            reg1 = std::regex(".aligns");
        }

        sampleName = regex_replace(sampleName, reg1, "");  // 'An-1.aligns' Delete '.aligns'
    }

    // Record sample name
    chrStartSynQryLocMap.sampleName = sampleName;

    // First initialize the sylLocVecOutMap
    for (const auto& [chromosome, locationSampleVec] : synLocSampleVecMap) {  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
        for (const auto& [synStart, synEnd, sampleVec] : locationSampleVec) {  // vector<tuple<refStart, refEnd, vector<sample> > >
            chrStartSynQryLocMap.chrStartSynQryLocMap[chromosome].emplace(synStart, make_tuple(0, 0));
        }
    }
    
    // Temporary syn coordinates
    string synChr = "";
    uint32_t synIdx = 0;  // Used to extract collinear coordinates
    uint32_t synStart = 0;
    uint32_t synEnd = 0;
    // Stores the coordinates of syn on qry
    uint32_t qrySynStart;
    uint32_t qrySynEnd;

    // Temporary chromosome numbers and coordinates
    string refChr = "";
    string qryChr = "";
    // Information on alignmenmt
    string AliRefStrand = "";
    string AliQryStrand = "";
    uint32_t aliRefStart = 0;
    uint32_t aliRefEnd = 0;
    uint32_t aliQryStart = 0;
    uint32_t aliQryEnd = 0;
    // Boolean value (false) used to determine whether to loop the alignment
    bool aliBool = false;
    // The coordinates used to record the syn are in the line of this alignment
    int32_t aliRowStartNum = 0;
    int32_t aliRowEndNum = INT32_MAX;

    // Line-specific information
    uint32_t refStart = 0;
    uint32_t refEnd = 0;
    string refSeq = "";
    uint32_t qryStart = 0;
    uint32_t qryEnd = 0;
    string qrySeq = "";

    // Constructed temporary bool value, not used, only submitted as a parameter
    bool whileBoolTmp = true;

    // open aligns file
    GzChunkReader GzChunkReaderClass(inputFileName);

    // Record loop index. even ->ref odd ->qry
    int64_t forIdx = -1;

    // read line
    string line;
    
    while (GzChunkReaderClass.read_line(line)) {
        // Skip empty line
        if (line.empty()) {
            continue;
        }
        
        if (
            line == "\n" || 
            line[0] == ' ' || 
            line[0] == '=' || 
            line.find("--   END alignment ") != string::npos
        ) {  // ' /=/\n/END' alignment is skipped
            continue;
        } else {  // Get the soon-Aligns coordinates
            if (line.find("-- Alignments between") != string::npos) {  // New comparison chromosome    -- Alignments between chr1 and chr1
                // split
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                refChr = lineVec[3];  // Temporary value of the ref chromosome
                qryChr = lineVec[5];  // Temporary value of the qry chromosome
                
                // Collinear information update
                synChr = refChr;
                synIdx = 0;
                synStart = 0;
                synEnd = 0;
                renew_syn_loc(
                    sampleName, 
                    synLocSampleVecMap, 
                    synChr, 
                    synIdx, 
                    synStart, 
                    synEnd, 
                    qrySynStart, 
                    qrySynEnd, 
                    whileBoolTmp
                );

                // Reset coordinate
                AliRefStrand = "";
                AliQryStrand = "";

                // Empty string
                line.clear();
                string().swap(line);
                continue;
            } else if (line.find("-- BEGIN alignment ") != string::npos) {  // New alignment coordinates    -- BEGIN alignment [ +1 1078 - 68996 | +1 1 - 67961 ]
                line = strip(line.erase(0, 21), '\n');  // Delete '-- BEGIN alignment [' contains 21 characters
                line = line.erase(line.size()-2, 2);  // Remove ']' a total of 2 characters, position in the last two

                vector<string> lineVec = split(line, " | ");  // Partition '+1 278-1703-1 2148-751 '

                // Get alignment direction, start and end information  ref
                string AliRefStrandTmp;
                uint32_t aliRefStartTmp;
                uint32_t aliRefEndTmp;
                tie(AliRefStrandTmp, aliRefStartTmp, aliRefEndTmp) = get_alignment_loc(
                    lineVec[0]
                );

                // Get alignment direction, start and end information  qry
                string AliQryStrandTmp;
                uint32_t aliQryStartTmp;
                uint32_t aliQryEndTmp;
                tie(AliQryStrandTmp, aliQryStartTmp, aliQryEndTmp) = get_alignment_loc(
                    lineVec[1]
                );

                // Reset the index of the sequence
                forIdx = -1;

                // Skip this read before the coordinates are collinear.
                if (aliRefEndTmp < synStart) {
                    aliBool = false;
                }
                // Determines whether to skip the reverse complementary sequence and skip this read.
                else if (!findRevBool && (AliRefStrandTmp == "-1" || AliQryStrandTmp == "-1")) {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Reverse alignment skipped -> " << synChr << " " << line << endl;
                    aliBool = false;
                } else if (refChr != qryChr) {
                    // ref and qry have different chromosomes, skip this read.
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Skipped alignment due to chromosomal mismatch -> " << refChr << " != " << qryChr << endl;
                    aliBool = false;
                } else {
                    // Reset coordinate
                    AliRefStrand = AliRefStrandTmp;
                    aliRefStart = aliRefStartTmp;
                    aliRefEnd = aliRefEndTmp;
                    AliQryStrand = AliQryStrandTmp;
                    aliQryStart = aliQryStartTmp;
                    aliQryEnd = aliQryEndTmp;
                    aliBool = true;

                    // Calculate the approximate line where syn is located,
                    aliRowStartNum = ((static_cast<int32_t>(synStart) - static_cast<int32_t>(aliRefStart))/49)*2-10;
                    aliRowEndNum = ((static_cast<int32_t>(synEnd) - static_cast<int32_t>(aliRefStart))/49)*2-10;
                }
            } else if (isdigit(line[0]) != 0 && aliBool && synEnd > 0) {
                // Only loops that start with a number and contain syn
                // If the map is done, again skip it, and when it's done it means 'synEnd==0'.
                forIdx++;  // The number of rows is even ref odd qry

                // If the line is smaller than the threshold, skip
                if (forIdx < aliRowStartNum) {
                    continue;
                }

                // Split string
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                if (( forIdx & 1 ) == 0) {  // even
                    // Reset coordinates (this)
                    refSeq = lineVec[1];
                    // Temporary ref sequence
                    string refSeqTmp = refSeq;
                    //Remove '.' to calculate length
                    refSeqTmp.erase(remove(refSeqTmp.begin(), refSeqTmp.end(), '.'), refSeqTmp.end());
                    if (AliRefStrand == "+1") {  // Forward alignment
                        refStart = stoul(lineVec[0]);
                        refEnd = refStart + refSeqTmp.size() - 1;
                    } else {  // Reverse alignment
                        refEnd = stoul(lineVec[0]);
                        refStart = refEnd - refSeqTmp.size() + 1;
                    }
                } else {  // odd
                    // Reset coordinates (this)
                    qrySeq = lineVec[1];
                    // Temporary ref sequence
                    string qrySeqTmp = qrySeq;
                    // Remove '.' to calculate length
                    qrySeqTmp.erase(remove(qrySeqTmp.begin(), qrySeqTmp.end(), '.'), qrySeqTmp.end());
                    if (AliQryStrand == "+1") {  // Forward alignment
                        qryStart = stoul(lineVec[0]);
                        qryEnd = qryStart + qrySeqTmp.size() - 1;
                    } else {  // Reverse alignment
                        qryEnd = stoul(lineVec[0]);
                        qryStart = qryEnd - qrySeqTmp.size() + 1;
                    }
                    
                    // If the colinear interval is before that interval, update the colinear coordinates, but to prevent an endless loop, so when 'synEnd==0' means that the dictionary has been traversed, skip
                    while (synEnd > 0 && synEnd < refStart) {
                        renew_syn_loc(
                            sampleName, 
                            synLocSampleVecMap, 
                            synChr, 
                            synIdx, 
                            synStart, 
                            synEnd, 
                            qrySynStart, 
                            qrySynEnd, 
                            whileBoolTmp
                        );

                        // Calculate the approximate line where syn is located,
                        aliRowStartNum = ((static_cast<int32_t>(synStart) - static_cast<int32_t>(aliRefStart))/49)*2-10;
                        aliRowEndNum = ((static_cast<int32_t>(synEnd) - static_cast<int32_t>(aliRefStart))/49)*2-10;

                        // Check whether syn is included in the comparison. Skip the lines that do not contain SYN
                        if (aliRefEnd < synStart) {  // Skip this read before the coordinates are collinear
                            aliBool = false;
                        }
                    }

                    // Determines whether to enter the while loop
                    bool whileBool = true;

                    // Find the coordinates and push them into the diagram
                    while (
                        ((refStart <= synStart && synStart <= refEnd) || 
                        (refStart <= synEnd && synEnd <= refEnd)) &&
                        whileBool
                    ) {
                        uint32_t qrySynStartTmp = find_qry_syn(
                            synStart, 
                            refStart, 
                            refEnd, 
                            refSeq, 
                            AliQryStrand, 
                            qryStart, 
                            qryEnd, 
                            qrySeq
                        );

                        uint32_t qrySynEndTmp = find_qry_syn(
                            synEnd, 
                            refStart, 
                            refEnd, 
                            refSeq, 
                            AliQryStrand, 
                            qryStart, 
                            qryEnd, 
                            qrySeq
                        );

                        // See if you find it. If you find it, add it to the total hash table
                        syn_all_loc_push(
                            chrStartSynQryLocMap, 
                            refChr, 
                            qrySynStartTmp, 
                            qrySynEndTmp, 
                            refStart, 
                            refEnd, 
                            qryStart, 
                            qryEnd, 
                            synStart, 
                            synEnd, 
                            qrySynStart, 
                            qrySynEnd, 
                            whileBool, 
                            aliRowStartNum, 
                            aliRowEndNum
                        );

                        // If both coordinates are found, or if synEnd has been reached, update the syn coordinates
                        if (
                            (qrySynStart > 0 && qrySynEnd > 0) || 
                            (refStart <= synEnd && synEnd <= refEnd) || 
                            synEnd < refStart
                        ) {
                            renew_syn_loc(
                                sampleName, 
                                synLocSampleVecMap, 
                                synChr, 
                                synIdx, 
                                synStart, 
                                synEnd, 
                                qrySynStart, 
                                qrySynEnd, 
                                whileBool
                            );
                            // Check whether syn is included in the comparison. Skip the lines that do not contain SYN
                            if (aliRefEnd < synStart) {  // Before the coordinates are collinear, skip this read.
                                aliBool = false;
                            }
                            // Calculate the approximate line where syn is located,
                            aliRowStartNum = ((static_cast<int32_t>(synStart) - static_cast<int32_t>(aliRefStart))/49)*2-10;
                            aliRowEndNum = ((static_cast<int32_t>(synEnd) - static_cast<int32_t>(aliRefStart))/49)*2-10;
                        }
                    }
                }
            }
        }
    }

    return chrStartSynQryLocMap;
}


/**
 * @brief Save the result
 * 
 * @param synLocSampleVecMap            build_syn_idx index to build   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
 * @param sampleChrStartSynQryLocMap    get_syn_coor Return value   map<sampleName, chrStartSynQryLocMap>   map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
 * @param outputFileName                Output file name
 * 
 * @return chrStartSynQryLocMap      synAllStructure
*/
int COOR::save_result(
    const map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > >& synLocSampleVecMap,
    const map<string, unordered_map<string, unordered_map<uint32_t, tuple<uint32_t, uint32_t> > > >& sampleChrStartSynQryLocMap,
    const string& outputFileName
) {
    SAVE SAVEClass(outputFileName);

    stringstream outStream; // Use stringstream instead of string concatenation
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // The cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);

    // save result
    for (const auto& [chromosome, locationSampleVec] : synLocSampleVecMap) {  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
        // Record the last collinear coordinate to prevent collinearity between intervals
        uint64_t preSynEnd = 0;

        for (const auto& [synStart, synEnd, sampleVec] : locationSampleVec) {  // vector<tuple<refStart, refEnd, vector<sample> > >
            // Otherwise, check whether there is an interval
            if(preSynEnd != 0 && synStart > preSynEnd + 1) {
                outStream << chromosome << '\t' << preSynEnd + 1 << '\t' << synStart - 1;
                
                for (const auto& it3 : sampleChrStartSynQryLocMap) {  // map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
                    outStream << '\t' << it3.first << '\t' << 0 << '\t' << 0;
                }
                outStream << '\n';
            }

            preSynEnd = synEnd;  // Record the coordinates of the previous collinearity

            // Information about the current node
            outStream << chromosome << '\t' << synStart << '\t' << synEnd;
            for (const auto& it3 : sampleChrStartSynQryLocMap) {  // map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
                const auto& chrStartSynQryLoc = it3.second.at(chromosome).at(synStart);
                outStream << '\t' << it3.first << '\t' << get<0>(chrStartSynQryLoc) << '\t' << get<1>(chrStartSynQryLoc);
            }
            outStream << '\n';

            if (outStream.tellp() >= CACHE_SIZE) {  // The cache size is 10mb
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // Clearing a stringstream
                outStream.str(string());
                outStream.clear();
            }
        }
    }

    if (outStream.tellp() > 0) {  // Write for the last time
        string outTxt = outStream.str();
        SAVEClass.save(outTxt);
        // Clearing a stringstream
        outStream.str(string());
        outStream.clear();
    }

    return 0;
}