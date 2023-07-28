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
    pair<string, map<string, vector<pair<int, int> > > > build_syn_index(
        const string & inputFileName
    );


    /**
        * @brief 找合集
        * 
        * @param chrLineSynVecMap          map<string, map<string, vector<pair<int, int> > > >, map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
        * 
        * @return pair(outChrStartEndLineVecMap, idxLineMap)   pair(map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>, 存储line的索引和名字)
    **/
    pair<map<string, map<int, vector<tuple<int, vector<string> > > > >, map<int, string> > syn_multiinter_find(
        const map<string, map<string, vector<pair<int, int> > > > & chrLineSynVecMap
    );
    

    /**
        * @brief 找集合_find
        * 
        * @param synVec      vector<pair<int, int > >, vector<pair<refStart, refEnd> >
        * 
        * @return tuple<int, int, vector<string> >  make_tuple(outStart, outEnd, lineName)
    **/
    tuple<int, int, vector<string> > loc_find(
        const vector<pair<int, int > > synVec, 
        const map<int, string> & idxLineMap
    );


    /**
        * @brief 存储结果
        * 
        * @param outChrStartEndLineVecMap      const map<chr, map<refStart, vector<tuple<refEnd, vector<lineName> > > > >
        * @param idxLineMap                  map<index, lineName>
        * @param outputName                  输出文件名
        * 
        * @return tuple<int, int, vector<string> >  make_tuple(outStart, outEnd, lineName)
    **/
    int save_result(
        const map<string, map<int, vector<tuple<int, vector<string> > > > > & outChrStartEndLineVecMap, 
        map<int, string> idxLineMap, 
        const string & outputName
    );

}

#endif