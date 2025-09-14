#ifndef COOR_HPP
#define COOR_HPP
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <regex>
#include <getopt.h>
#include <future>
#include <malloc.h>
#include <iterator>

#include "ThreadPool.hpp"
#include "zlib.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "save.hpp"
#include "GzChunkReader.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

// Global variable
extern int thresholdLength;  // Collinear coordinates ref and qry length ratio threshold

// Whether to debug code
extern bool debugCoor;

void help_coor(char* argv[]);
int main_coor(int argc, char* argv[]);

namespace COOR
{
    struct synAllStructure
    {
        string sampleName;
        unordered_map<string, unordered_map<uint32_t, tuple<uint32_t, uint32_t> > > chrStartSynQryLocMap;  // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
    };


    /**
     * @brief Parse the '-- BEGIN alignment 'field
     * 
     * @param informationTmp    '+1 278 - 1703'
     * 
     * @return tuple<strand, start, end>
    */
    tuple<string, uint32_t, uint32_t> get_alignment_loc(string informationTmp);


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
    uint32_t find_qry_syn(
        const uint32_t & synLoc, 
        uint32_t refStart, 
        uint32_t refEnd, 
        const string & refSeq, 
        const string & AliQryStrand, 
        uint32_t qryStart, 
        uint32_t qryEnd, 
        const string & qrySeq
    );


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
    uint32_t syn_all_loc_push(
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
    );


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
    int renew_syn_loc(
        const string& sampleName, 
        const map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > >& synLocSampleVecMap, 
        const string& synChr, 
        uint32_t& synIdx, 
        uint32_t& synStart, 
        uint32_t& synEnd, 
        uint32_t& qrySynStart, 
        uint32_t& qrySynEnd, 
        bool& whileBool
    );


    /**
     * @brief The syn coordinate index common to all samples in the component
     * 
     * @param inputFileName         syn_loc output file
     * 
     * @return synLocSampleVecMap   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    */
    map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > > build_syn_idx(const string & inputFileName);


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
    synAllStructure get_syn_coor(
        string sampleName, 
        const string& inputFileName, 
        const map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > >& synLocSampleVecMap, 
        const bool& findRevBool
    );


    /**
     * @brief Save the result
     * 
     * @param synLocSampleVecMap            build_syn_idx index to build   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
     * @param sampleChrStartSynQryLocMap    get_syn_coor Return value   map<sampleName, chrStartSynQryLocMap>   map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
     * @param outputFileName                Output file name
     * 
     * @return chrStartSynQryLocMap      synAllStructure
    */
    int save_result(
        const map<string, vector<tuple<uint32_t, uint32_t, vector<string> > > >& synLocSampleVecMap,
        const map<string, unordered_map<string, unordered_map<uint32_t, tuple<uint32_t, uint32_t> > > >& sampleChrStartSynQryLocMap,
        const string& outputFileName
    );
}

#endif