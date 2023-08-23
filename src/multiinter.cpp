// g++ -c multiinter.cpp -o multiinter -lz -O3
#include "../include/multiinter.hpp"

using namespace std;

// Debug the code
bool debugMultiinter = false;

int main_multiinter(int argc, char* argv[])
{
    // Enter the file list and name syri.out
    vector<string> inputFiles;
    vector<string> inputTitles;

    // Check whether there is a titlename
    bool haveTitles = false;

    // Output file name
    string outputFileName;

    //Parse command line options
    if(argc <= 2)
    {
        help_multiinter(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_multiinter(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-i", 2, parameterLength)) || 
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
        else if(PARAMETER_CHECK("--debug", 7, parameterLength)) {
            if ((i) < argc) {
                debugMultiinter = true;
            }
        }
    }

    if (argc <= 2 || inputFiles.size() <= 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Only " << inputFiles.size() << " syri.out file was specified. Nothing to extract, exiting." << endl;
        help_multiinter(argv);
        return 1;
    }
    if ((haveTitles == true) && (inputFiles.size() != inputTitles.size())) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of file titles (-n) does not match the number of files (-i)." << endl;
        help_multiinter(argv);
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ..." << endl;

    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index ..." << endl;

    // Stores collinear information for all lines
    map<string, map<string, vector<pair<uint32_t, uint32_t> > > > chrLineSynVecMap;  // map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
    
    uint32_t indexTmp = 0;  // Record the lineName index

    for(auto it1 : inputFiles) {
        string lineName;
        map<string, vector<pair<uint32_t, uint32_t> > > chrSynVecMap;
        tie(lineName, chrSynVecMap) = MULTIINTER::build_syn_index(it1);

        // Reassign inputTitles if they exist
        if (haveTitles) {
            lineName = inputTitles[indexTmp];
        }
        
        for(auto it2 : chrSynVecMap) {
            string chromosome = it2.first;
            chrLineSynVecMap[chromosome][lineName] = it2.second;
        }

        indexTmp++;  // Index overlay
    }

    /* ************************************ Find intersection ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find intersection ..." << endl;
    map<string, map<uint32_t, vector<tuple<uint32_t, vector<string> > > > > outChrStartEndLineVecMap;  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>
    map<uint32_t, string> idxLineMap;
    tie(outChrStartEndLineVecMap, idxLineMap) = MULTIINTER::syn_multiinter_find(
        chrLineSynVecMap
    );

    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Results are being saved to '" << outputFileName << "'" << endl;
    MULTIINTER::save_result(
        outChrStartEndLineVecMap, 
        idxLineMap, 
        outputFileName
    );

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ..." << endl;
    
    return 0;
}


void help_multiinter(char* argv[])
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i FILE1 FILE2 .. FILEn [options]" << endl
       << "identifies common intervals (SYNAL) among multiple syri.out files" << endl
       << endl
       << "required arguments:" << endl
       << "    -i, --inputs        FILE      list of output files of syri (syri.out), one for multiple mate" << endl
       << endl
       << "optional arguments:" << endl
       << "    -n, --names         STRING    list of names to describe each file in -i, one for multiple mate" << endl
       << "    -o, --output        FILE      output syntenic intersection to FILE [stdout]" << endl
       << "    --debug                       debug code (false by default)" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}


/**
 * @brief build syn index for syri.out
 * 
 * @param inputFileName   the output of syri
 * 
 * @return pair<lineName, chrSynVecMap>        pair<lineName, map<chromosome, vector<pair<refStart, refEnd>>>>
**/
pair<string, map<string, vector<pair<uint32_t, uint32_t> > > > MULTIINTER::build_syn_index(
    const string & inputFileName
)
{
    // The name of the breed
    string lineName = inputFileName;

    vector<string> inputFileNameVecTmp = split(inputFileName, "/");  // Path splitting
    lineName = inputFileNameVecTmp.back();  // The last one

    string prefix = ".syri.out";
    auto iter1 = lineName.find(prefix);
    if(iter1 != string::npos)
    {
        lineName.replace(lineName.begin() + iter1, lineName.begin() + iter1 + prefix.size(), "");
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building index: " << inputFileName << " ...\n";  // print log

    map<string, vector<pair<uint32_t, uint32_t> > > chrSynVecMap;  // map<chromosome, vector<pair<refStart, refEnd>>>

    // open file
    GzChunkReader GzChunkReaderClass(inputFileName);

    // Coordinates of the previous collinearity
    string preChromosome;
    uint32_t PreRefStart = 0;
    uint32_t PreRefEnd = 0;

    // read line
    string line;

    while (GzChunkReaderClass.read_line(line))
    {
        // Skip the blank line
        if (line.empty() || line.find("SYNAL") == string::npos)
        {
            continue;
        }

        // Split
        std::istringstream iss(line);
        vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

        string chromosome = lineVec[0];
        uint32_t refStart = stoul(lineVec[1]);
        uint32_t refEnd = stoul(lineVec[2]);

        if(preChromosome != chromosome && preChromosome.size() > 0)  // If you change chromosomes/Skip the first chromosome. Zero clearing
        {
            chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));
            PreRefStart = refStart;
            PreRefEnd = refEnd;
        }
        else
        {
            if(PreRefStart == 0)  // First syn on chromosome 1
            {
                PreRefStart = refStart;
                PreRefEnd = refEnd;
            }
            else if(refStart <= PreRefEnd)  // The previous one contains the interval
            {
                if(refEnd > PreRefEnd)  // Part contains
                {
                    PreRefStart = PreRefStart;
                    PreRefEnd = refEnd;
                }
                else  // Fully included
                {
                    PreRefStart = PreRefStart;
                    PreRefEnd = PreRefEnd;
                }
            }
            else  // No overlap
            {
                chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));
                PreRefStart = refStart;
                PreRefEnd = refEnd;
            }
        }

        preChromosome = chromosome;
    }

    // Last coordinate added
    chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));

    return make_pair(lineName, chrSynVecMap);
}


/**
 * @brief Find set
 * 
 * @param chrLineSynVecMap          map<string, map<string, vector<pair<uint32_t, uint32_t> > > >, map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
 * 
 * @return pair(outChrStartEndLineVecMap, idxLineMap)   pair(map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>, store the index and name of the line)
**/
pair<map<string, map<uint32_t, vector<tuple<uint32_t, vector<string> > > > >, map<uint32_t, string> > MULTIINTER::syn_multiinter_find(
    const map<string, map<string, vector<pair<uint32_t, uint32_t> > > > & chrLineSynVecMap
)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Searching ..." << endl;

    // save result
    map<string, map<uint32_t, vector<tuple<uint32_t, vector<string> > > > > outChrStartEndLineVecMap;  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>

    map<uint32_t, string> idxLineMap;  // Store the index and name of the line

    // Ergodic chromosome
    for(auto it1 : chrLineSynVecMap)  // map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
    {
        // it1.first -> chromosome
        string chromosome = it1.first;

        // it1.second -> map<lineName, vector<pair<refStart, refEnd>>>

        uint32_t lineIdx = 0;  // lineName Index
        
        // All collinear lists and index information corresponding to chromosomes are extracted, and the first coordinate is assigned as the first of the syn list
        map<string, tuple<vector<pair<uint32_t, uint32_t> >, uint32_t, pair<uint32_t, uint32_t> > > LineSynVecIdxSynMap;  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
        for(auto it2 : it1.second)  // map<lineName, vector<pair<refStart, refEnd>>>
        {
            // it2.first -> lineName
            // it2.second -> vector<pair<refStart, refEnd>>

            LineSynVecIdxSynMap[it2.first] = make_tuple(it2.second, 0, it2.second[0]);

            idxLineMap[lineIdx] = it2.first;  // Store line information
            lineIdx++;  // Index plus one
        }


        // To find, first construct a temporary collinear list, storing the first collinear interval of each sample
        vector<pair<uint32_t, uint32_t> > synVec;  // Store collinear intervals, find a set
        
        uint32_t lineNum = idxLineMap.size();  // sample Indicates the number of samples to determine whether to terminate the loop
        uint32_t endNum = 0;  // Used to determine whether to terminate the while loop

        for(auto it1 : LineSynVecIdxSynMap)  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
        {
            // it1.first -> lineName
            string lineName = it1.first;

            // it1.second -> tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> >
            synVec.push_back(get<2>(it1.second));

            // If debugging code, print log
            if (debugMultiinter)
            {
                cerr << chromosome << " " << lineName << ":" << get<2>(it1.second).first << "-" << get<2>(it1.second).second << endl;
            }
            
            if(get<1>(it1.second) == get<0>(it1.second).size())  // Check whether the loop is complete
            {
                endNum++;
            }
        }
        
        while(endNum < lineNum)  // End the loop when endNum == lineNum
        {
            endNum = 0;  // Reset count

            // Find intersection
            uint32_t outRefStart = 0;
            uint32_t outRefEnd = 0;
            vector<string> outLineVec;
            tie(outRefStart, outRefEnd, outLineVec) = loc_find(
                synVec, 
                idxLineMap
            );

            // Print log if debugging code
            if (debugMultiinter)
            {
                cerr << chromosome << " loc:" << outRefStart << "-" << outRefEnd << endl << chromosome << " outLine:";

                for (auto outLine : outLineVec)
                {
                    cerr << " " << outLine;
                }
                cerr << endl << endl;
            }

            // Store intersection results
            outChrStartEndLineVecMap[chromosome][outRefStart].push_back(make_tuple(outRefEnd, outLineVec));

            // Empty coordinates
            synVec.clear();
            vector<pair<uint32_t, uint32_t > >().swap(synVec);

            // Update collinear indexes and coordinates
            for(auto it1 : LineSynVecIdxSynMap)  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
            {
                // it1.first -> lineName
                string lineName = it1.first;

                // Collinear coordinates corresponding to varieties, used to determine whether to update the coordinates
                vector<pair<uint32_t, uint32_t> > lineSynVec = get<0>(it1.second);
                uint32_t lineSynIdx = get<1>(it1.second);
                pair<uint32_t, uint32_t> lineSycLocTmp = get<2>(it1.second);

                // Update coordinate
                if(lineSycLocTmp.first >= outRefEnd)  // If this collinear interval is not used, go straight to the next variety
                {
                    synVec.push_back(lineSycLocTmp);  // Update the coordinates and use the coordinates of the other one
                }
                else if(outRefEnd >= lineSycLocTmp.second)  // The collinear coordinates are used up
                {
                    lineSynIdx++;
                    if(lineSynIdx >= lineSynVec.size())  // Check whether the loop is complete
                    {
                        lineSycLocTmp = make_pair(0, 0);  // Total return to zero
                        endNum++;
                    }
                    else
                    {
                        lineSycLocTmp = lineSynVec[lineSynIdx];  // Collinearity points to the next coordinate
                    }

                    synVec.push_back(lineSycLocTmp);
                    LineSynVecIdxSynMap[lineName] = make_tuple(lineSynVec, lineSynIdx, lineSycLocTmp);  // Update the coordinates of the hash table
                }
                else  // Collinear coordinates are used half of the time
                {
                    lineSycLocTmp.first = outRefEnd;

                    synVec.push_back(lineSycLocTmp);  // Update coordinate
                    LineSynVecIdxSynMap[lineName] = make_tuple(lineSynVec, lineSynIdx, lineSycLocTmp);  // Update the coordinates of the hash table
                }

                // If debugging code, print log
                if (debugMultiinter)
                {
                    cerr << chromosome << " " << lineName << ":" << lineSycLocTmp.first << "-" << lineSycLocTmp.second << " endNum:" << endNum << endl;
                }
            }
        }
    }

    return make_pair(outChrStartEndLineVecMap, idxLineMap);
}


/**
 * @brief Find set _find
 * 
 * @param synVec      vector<pair<uint32_t, uint32_t > >, vector<pair<refStart, refEnd> >
 * 
 * @return tuple<uint32_t, uint32_t, vector<string> >  make_tuple(outStart, outEnd, lineName)
**/
tuple<uint32_t, uint32_t, vector<string> > MULTIINTER::loc_find(
    const vector<pair<uint32_t, uint32_t > > synVec, 
    const map<uint32_t, string> & idxLineMap
)
{
    // Collinear interval
    uint32_t outStart = UINT32_MAX;
    uint32_t outEnd = UINT32_MAX;

    vector<string> lineNameVec;  // Index of varieties corresponding to collinear intervals
    

    for(auto it1 : synVec)  // Find the lowest collinear coordinates, which must be in the starting coordinates
    {
        if(it1.first == 0 && it1.second == 0) continue;  // If both are 0, the line loop is complete

        outStart = min(outStart, it1.first);  // The lowest value in all collinearity must be in the starting coordinates
    }

    for(auto it1 : synVec)  // Find the collinearity second smallest coordinates
    {
        if(it1.first == 0 && it1.second == 0) continue;  // If both are 0, the line loop is complete

        if (it1.first != outStart)  // It cannot be the minimum value
        {
            outEnd = min(outEnd, it1.first);
        }
        outEnd = min(outEnd, it1.second);
    }

    if (outEnd == UINT32_MAX)  // Prevent situations where the end is the same as the start
    {
        outEnd = outStart;
    }
    

    uint32_t lineIdx = 0;
    for(auto it1 : synVec)  // Find the minimum interval corresponding to the variety information
    {
        if(it1.first == 0 && it1.second == 0)  // If both are 0, the line loop is complete
        {
            lineIdx++;
            continue;
        }

        if(it1.first >= outStart && it1.first < outEnd)
        {
            string lineName;
            auto iter = idxLineMap.find(lineIdx);
            if(iter != idxLineMap.end())
            {
                lineName = iter->second;
            }
            else
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "keyError: " << lineIdx << " not in idxLineMap." << endl;
                exit(1);
            }
            lineNameVec.push_back(lineName);
        }
        lineIdx++;
    }
    
    return make_tuple(outStart, outEnd, lineNameVec);
}


/**
 * @brief save result
 * 
 * @param outChrStartEndLineVecMap    const map<chr, map<refStart, vector<tuple<refEnd, vector<lineName> > > > >
 * @param idxLineMap                  map<index, lineName>
 * @param outputName                  Output file name
 * 
 * @return tuple<uint32_t, uint32_t, vector<string> >  make_tuple(outStart, outEnd, lineName)
**/
int MULTIINTER::save_result(
    const map<string, map<uint32_t, vector<tuple<uint32_t, vector<string> > > > > & outChrStartEndLineVecMap, 
    map<uint32_t, string> idxLineMap, 
    const string & outputName
)
{
    SAVE SAVEClass(outputName);

    stringstream outStream; // Use stringstream instead of string concatenation
    static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // The cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);

    // Add lineName to the list
    vector<string> lineNameVec;
    for(auto it : idxLineMap)
    {
        lineNameVec.push_back(it.second); 
    }

    // save result
    for(auto it1 : outChrStartEndLineVecMap)  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName> > > > >
    {
        for(auto it2 : it1.second)  // map<refStart, vector<tuple<refEnd, vector<lineName> > > >
        {
            for(auto it3 : it2.second)  // vector<tuple<refEnd, vector<lineName> > >
            {
                outStream << it1.first << "\t" << to_string(it2.first) << "\t" + to_string(get<0>(it3)) << "\t" << to_string(get<1>(it3).size()) << "\t";

                // Collinearity contains index of breed names
                uint32_t idxTmp=0;

                // Output (Column 5 and beyond)
                vector<string> lineNameVecTmp;
                vector<int> boolVecTmp;
                for (size_t i = 0; i < lineNameVec.size(); i++)  // vector<lineName>
                {
                    if(lineNameVec[i] == get<1>(it3)[idxTmp])  // Contains the corresponding variety
                    {
                        lineNameVecTmp.push_back(lineNameVec[i]);
                        boolVecTmp.push_back(1);
                        idxTmp++;
                    }
                    else
                    {
                        boolVecTmp.push_back(0);
                    }

                    // If you go to the end of the vector, jump out of the for loop
                    if (idxTmp == get<1>(it3).size())
                    {
                        break;
                    }  
                }

                // If the result is empty, skip the site
                if (lineNameVecTmp.size() == 0)
                {
                    continue;
                }
                
                outStream << join(lineNameVecTmp, ",") << "\t" + join(boolVecTmp, "\t") << "\n";

                if (outStream.tellp() >= CACHE_SIZE)  // The cache size is 10mb
                {
                    string outTxt = outStream.str();
                    SAVEClass.save(outTxt);
                    // Clearing a stringstream
                    outStream.str(string());
                    outStream.clear();
                }
            }
        }
    }

    if (outStream.tellp() > 0)  // Write for the last time
    {
        string outTxt = outStream.str();
        SAVEClass.save(outTxt);
        // Clearing a stringstream
        outStream.str(string());
        outStream.clear();
    }
    
    return 0;
}