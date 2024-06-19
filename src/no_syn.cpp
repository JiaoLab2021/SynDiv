// g++ no_syn.cpp -o no_syn -lz -lpthread
#include "../include/no_syn.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

// debug
bool debugNoSyn = false;

int main_no_syn(int argc, char* argv[])
{
    // syntenic coordinates
    string coorFileName;

    // chromosome length file
    vector<string> lengthsVec;
    vector<string> lengthsTitles;
    // Determine whether there is titlename
    bool haveTitles = false;

    // output file name
    string outputFileName;

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
                debugNoSyn = true;
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

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;


   /* ************************************ Build Chromosome Length Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Chromosome Length Index ..." << endl;
    NOSYN::LENGTH lengthClass(lengthsVec, lengthsTitles);
    lengthClass.index_lengths();

    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index ..." << endl;
    NOSYN::NOSYNCOOR NoSynCoor(coorFileName);
    NoSynCoor.open_coor();

    /* ************************************ Find no-syntenic Coordinates on query genomes ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find No-syntenic Coordinates on query genomes ..." << endl;
    NoSynCoor.find_no_syn_coor(
        lengthClass
    );

    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;
    SAVE SaveClass(outputFileName);
    SaveClass.save(NoSynCoor.outTxt);  // save result

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ..." << endl;
    
    return 0;
}


void help_no_syn(char* argv[])
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " --coor FILE --lengths FILE1 FILE2 .. FILEn [options]" << endl
       << "retrieve non-syntenic coordinates" << endl
       << endl
       << "required arguments:" << endl
       << "    --coor              FILE      syntenic coordinates, output file of 'SynDiv_c coor'" << endl
       << "    --lengths           FILE      chromosome length files for query genomes (format: chr\tlength), one for multiple mate" << endl
       << endl
       << "optional arguments:" << endl
       << "    -n, --names         STRING    list of names to describe each file in --lengths, one for multiple mate" << endl
       << "    -o, --output        FILE      output coordinates to FILE [stdout]" << endl
       << "    --debug                       debug code (false by default)" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}


/**
 * @brief Find non-syntenic coordinates
 * 
 * @param lengthClass  Sample chromosome length class
 * 
 * @return void
**/
void NOSYN::NOSYNCOOR::find_no_syn_coor(
    LENGTH & lengthClass
)
{   
    // Traverse syntenic coordinates
    for (auto iter1 = coorChrLociSampleLociMap.begin(); iter1 != coorChrLociSampleLociMap.end(); iter1++) {  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
        string chrTmp = iter1->first;  //chromosome

        // The number of non-synteny determines whether the coordinates of the reference genome need to be changed to prevent repeated addition of coordinates.
        uint32_t preNoSynNum = 0;

        uint32_t refStart;  // Start of reference genome
        uint32_t refEnd;  // End of reference genome

        // Used to determine whether the non-synteny coordinates of the sample in a small area have been found
        unordered_map<string, bool> sampleNoSynFindBoolMap;  // map<sample, bool>   true->found it   false->did not find

        // Initialize whether the non-synteniy coordinates of each sample are found
        for (auto iter1 : sampleNameVec) {  // map<sample, bool>
            if (iter1 == "reference") {
                continue;
            }

            sampleNoSynFindBoolMap[iter1] = false;
        }


        for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {  // map<refStart, map<sample, tuple(start, end)> >
            uint32_t refStartTmp;  // Start of reference genome
            uint32_t refEndTmp;  // End of reference genome

            // The starting and ending positions of the node
            auto findIter1 = iter2->second.find("reference");
            if (findIter1 != iter2->second.end()) {
                refStartTmp = get<0>(findIter1->second);
                refEndTmp = get<1>(findIter1->second);
            } else {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: 'reference' not present in 'coorChrLociSampleLociMap' -> " << iter2->first << endl;
                exit(1);
            }
            
            // The number of non-synteny determines whether the coordinates of the reference genome need to be changed to prevent repeated addition of coordinates.
            uint32_t noSynNum = 0;

            // If the previous nodes are all synteny, update the coordinates of the reference genome
            if (preNoSynNum  == 0) {
                refStart = refStartTmp;
                refEnd = refEndTmp;
            }

            for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {  // map<sample, tuple(start, end)>
                string sampleTmp = iter3->first;  // sample name
                uint32_t StartTmp = get<0>(iter3->second);  // start
                uint32_t EndTmp = get<1>(iter3->second);  // end

                // Direct skipping of reference genome
                if (sampleTmp == "reference") {
                    continue;
                }

                // Recorded as non-synteny
                if (StartTmp == 0 || EndTmp == 0) {
                    // Add 1 to the number of non-synteny
                    noSynNum++;

                    // Find non-synteny locations
                    // Initialize the dictionary first
                    auto& startSampleLociVecMap = chrStartSampleLociVecMap.emplace(chrTmp, map<uint32_t, map<string, vector<tuple<uint32_t, uint32_t> > > >()).first->second;  // map<string, map<uint32_t, map<string, vector<tuple<uint32_t, uint32_t> > > > >
                    auto& sampleLociVecMap = startSampleLociVecMap.emplace(refStart, map<string, vector<tuple<uint32_t, uint32_t> > >()).first->second;  // map<uint32_t, map<string, vector<tuple<uint32_t, uint32_t> > > >

                    // If the sample coordinates are found, skip
                    if (sampleNoSynFindBoolMap[sampleTmp]) {  // map<sample, findBool>
                        continue;
                    }

                    // Temporary coordinates
                    uint32_t noSynStartTmp = 1;
                    uint32_t noSynEndTmp = 1;

                    auto startIter = iter2;  // start iter

                    if (iter2 == iter1->second.begin()) {  // If it is the first node
                        noSynStartTmp = 1;
                    } else {
                        startIter = iter2;
                        startIter--;  // Move up one position

                        // If the previous node is still 0 and is not the first node, continue to move forward.
                        while (startIter != iter1->second.begin() && get<1>(startIter->second[sampleTmp]) == 0) {
                            startIter--;  // Move up one position
                        }

                        // Assign the previous coordinate
                        if (get<1>(startIter->second[sampleTmp]) != 0) {
                            noSynStartTmp = get<1>(startIter->second[sampleTmp]) + 1;
                        }
                    }

                    auto endIter = iter2;  // end iter
                    endIter++;  // Move down one position

                    // If the next node is still 0 and is not the last node, continue to move backward.
                    while (endIter != iter1->second.end() && get<0>(endIter->second[sampleTmp]) == 0) {
                        endIter++;  // Move down one position
                    }

                    // If it is the last node, it is the length of the chromosome
                    if (endIter == iter1->second.end()) {
                        // chromosome length
                        noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);
                    } else {
                        // Starting coordinate of next node -1
                        noSynEndTmp = get<0>(endIter->second[sampleTmp]) - 1;
                    }

                    // start must be less than end
                    if (noSynStartTmp < noSynEndTmp) {
                        // initialize dictionary
                        auto& tupVec = sampleLociVecMap.emplace(sampleTmp, vector<tuple<uint32_t, uint32_t> >()).first->second;  // map<string, vector<tuple<uint32_t, uint32_t> > >

                        // If the length becomes shorter instead of longer, it means the coordinates are wrong and need to be corrected.
                        int64_t noSynLenTmp = static_cast<int64_t>(noSynEndTmp) - static_cast<int64_t>(noSynStartTmp);
                        // start
                        if (noSynLenTmp > 10000 && startIter != iter1->second.begin()) {
                            for (size_t i = 0; i < 20; i++) {
                                startIter--;
                                
                                if (
                                    get<1>(startIter->second[sampleTmp]) != 0 && 
                                    abs(static_cast<int64_t>(noSynEndTmp) - static_cast<int64_t>(get<1>(startIter->second[sampleTmp]))) - noSynLenTmp < -1000
                                ) {
                                    noSynStartTmp = get<1>(startIter->second[sampleTmp]) + 1;
                                    noSynLenTmp = abs(static_cast<int64_t>(noSynEndTmp) - static_cast<int64_t>(noSynStartTmp));
                                    break;
                                }

                                if (startIter == iter1->second.begin()) break;
                            }
                        }

                        // end
                        if (noSynLenTmp > 10000 && endIter != iter1->second.end()) {
                            for (size_t i = 0; i < 20; i++) {
                                endIter++;
                                if (endIter == iter1->second.end()) break;

                                if (
                                    get<0>(endIter->second[sampleTmp]) != 0 && 
                                    abs(static_cast<int64_t>(get<0>(endIter->second[sampleTmp])) - static_cast<int64_t>(noSynStartTmp)) - noSynLenTmp < -1000
                                ) {
                                    noSynEndTmp = get<0>(endIter->second[sampleTmp]) - 1;
                                    break;
                                }
                            }
                        }
                        
                        tupVec.push_back(make_tuple(min(noSynStartTmp, noSynEndTmp), max(noSynStartTmp, noSynEndTmp)));  // Record non-collinear coordinates

                        sampleNoSynFindBoolMap[sampleTmp] = true;  // Record that the sample was found
                    }
                } else {
                    // Determine whether it is the last node
                    auto iter2Tmp = iter2;  // Temporary node
                    iter2Tmp++;  // Move down one position

                    // If it is the first node, determine whether the comparison starts from 1. If it is collinear with the reference genome, there is no need to look for the first half of the chromosome.
                    if (iter2 == iter1->second.begin() && StartTmp == 0) {
                        // initialize dictionary
                        auto& sampleLociVecMap = chrStartSampleLociVecMap[chrTmp][refStart];
                        if (sampleLociVecMap.find(sampleTmp) == sampleLociVecMap.end()) {
                            sampleLociVecMap[sampleTmp];
                        }
                        sampleLociVecMap[sampleTmp].push_back(make_tuple(1, StartTmp - 1));
                    } else if (iter2Tmp == iter1->second.end()) {  // If it is the last node, determine whether it reaches the end of the chromosome
                        // chromosome length
                        uint32_t noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);
                        if (EndTmp < noSynEndTmp) {
                            // initialize dictionary
                            auto& sampleLociVecMap = chrStartSampleLociVecMap[chrTmp][refStart];
                            if (sampleLociVecMap.find(sampleTmp) == sampleLociVecMap.end()) {
                                sampleLociVecMap[sampleTmp];
                            }
                            sampleLociVecMap[sampleTmp].push_back(make_tuple(EndTmp + 1, noSynEndTmp));
                        }
                    }
                    
                    // Update The sample was not found
                    sampleNoSynFindBoolMap[sampleTmp] = false;
                }
            }
            // Update the number of non-synteny
            preNoSynNum = noSynNum;
        }
    }


    // saving the result
    stringstream outStream;  // Use stringstream instead of string concatenation
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10;  // Cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);

    // Convert result to string
    for (const auto& [chr, startSampleLociVecMap] : chrStartSampleLociVecMap) {  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
        for (const auto& [refStart, sampleLociVecMap] : startSampleLociVecMap) {  // map<refStart, map<sample, vector<tuple<start, end> > > >
            if (sampleLociVecMap.empty()) continue;

            outStream << chr << "\t" << to_string(refStart);

            for (const auto& [sample, lociVec] : sampleLociVecMap) {  // map<sample, vector<tuple<start, end> > >
                outStream << "\t" << sample << ":";
            
                for (size_t i = 0; i < lociVec.size(); ++i) {
                    if (i > 0) outStream << ";";

                    outStream << to_string(get<0>(lociVec[i])) << "-" << to_string(get<1>(lociVec[i]));
                }
                outStream << "\t";
            }
            outStream << "\n";
        }

        if (outStream.tellp() >= CACHE_SIZE)  // Cache size is 10mb
        {
            outTxt += outStream.str();
            // Clear stringstream
            outStream.str(string());
            outStream.clear();
        }
    }

    if (outStream.tellp() > 0) {  // Write for the last time
        outTxt += outStream.str();
        // Clear stringstream
        outStream.str(string());
        outStream.clear();
    }
}