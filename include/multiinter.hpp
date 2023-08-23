#ifndef MULTIINTER_HPP
#define MULTIINTER_HPP
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include <climits>
#include <cstring>
#include <getopt.h>

#include "save.hpp"
#include "GzChunkReader.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

extern bool debugMultiinter;

void help_multiinter(char* argv[]);
int main_multiinter(int argc, char* argv[]);

namespace MULTIINTER
{
    /**
     * @brief build syn index for syri.out
     * 
     * @param inputFileName   the output of syri
     * 
     * @return pair<lineName, chrSynVecMap>        pair<lineName, map<chromosome, vector<pair<refStart, refEnd>>>>
    **/
    pair<string, map<string, vector<pair<uint32_t, uint32_t> > > > build_syn_index(
        const string & inputFileName
    );


    /**
     * @brief Find set
     * 
     * @param chrLineSynVecMap          map<string, map<string, vector<pair<uint32_t, uint32_t> > > >, map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
     * 
     * @return pair(outChrStartEndLineVecMap, idxLineMap)   pair(map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>, store the index and name of the line)
    **/
    pair<map<string, map<uint32_t, vector<tuple<uint32_t, vector<string> > > > >, map<uint32_t, string> > syn_multiinter_find(
        const map<string, map<string, vector<pair<uint32_t, uint32_t> > > > & chrLineSynVecMap
    );
    

    /**
     * @brief Find set _find
     * 
     * @param synVec      vector<pair<uint32_t, uint32_t > >, vector<pair<refStart, refEnd> >
     * 
     * @return tuple<uint32_t, uint32_t, vector<string> >  make_tuple(outStart, outEnd, lineName)
    **/
    tuple<uint32_t, uint32_t, vector<string> > loc_find(
        const vector<pair<uint32_t, uint32_t > > synVec, 
        const map<uint32_t, string> & idxLineMap
    );


    /**
     * @brief save result
     * 
     * @param outChrStartEndLineVecMap    const map<chr, map<refStart, vector<tuple<refEnd, vector<lineName> > > > >
     * @param idxLineMap                  map<index, lineName>
     * @param outputName                  Output file name
     * 
     * @return tuple<uint32_t, uint32_t, vector<string> >  make_tuple(outStart, outEnd, lineName)
    **/
    int save_result(
        const map<string, map<uint32_t, vector<tuple<uint32_t, vector<string> > > > > & outChrStartEndLineVecMap, 
        map<uint32_t, string> idxLineMap, 
        const string & outputName
    );

}

#endif