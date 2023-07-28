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

#include "ThreadPool.hpp"
#include "zlib.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "save.hpp"
#include "GzChunkReader.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

// 全局变量
extern int thresholdLength;  // 共线性坐标ref和qry长度比阈值

// 是否调试代码
extern bool debugCoor;

void help_coor(char* argv[]);
int main_coor(int argc, char* argv[]);

namespace COOR
{
    struct synAllStructure
    {
        string sampleName;
        unordered_map<string, unordered_map<int, tuple<int, int> > > chrStartSynQryLocMap;  // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
    };


    /**
     * @brief 解析 '-- BEGIN alignment ' 字段
     * 
     * @param informationTmp    '+1 278 - 1703'
     * 
     * @return tuple<strand, start, end>
    */
    tuple<string, int, int> get_alignment_loc(string informationTmp);


    /**
     * @brief 找qry上共线性的坐标
     * 
     * @param synLoc        共线性的坐标
     * @param refStart      ref的起始
     * @param refEnd        ref的终止
     * @param refSeq        ref的序列
     * @param qryStrand     qry比对方向
     * @param qryStart      qry的起始
     * @param qryEnd        qry的终止
     * @param qrySeq        qry的序列
     * 
     * @return qrySynLoc    qry上syn的坐标， 0-没找到
    */
    int find_qry_syn(
        const int & synLoc, 
        int refStart, 
        int refEnd, 
        const string & refSeq, 
        const string & AliQryStrand, 
        int qryStart, 
        int qryEnd, 
        const string & qrySeq
    );


    /**
     * @brief 找qry上共线性的坐标
     * 
     * @param chrStartSynQryLocMap  synAllStructure
     * @param refChr                染色体号
     * @param qrySynStartTmp        find_qry_syn找的qrySynStart
     * @param qrySynEndTmp          find_qry_syn找的qrySynEnd
     * @param refStart              该比对行ref的起始
     * @param refEnd                该比对行ref的终止
     * @param qryStart              该比对行qry的起始
     * @param qryEnd                带比对行qry的终止
     * @param synStart              目前要找的syn的ref起始
     * @param synEnd                目前要找的syn的ref终止
     * @param qrySynStart           syn的qry起始
     * @param qrySynEnd             syn的qry终止
     * @param whileBool             判断是否还需要while循环
     * @param aliRowStartNum        syn起始位置在该alignment的大致行数
     * @param aliRowEndNum          syn终止位置在该alignment的大致行数
     * 
     * @return qrySynLoc    qry上syn的坐标， 0-没找到
    */
    int syn_all_loc_push(
        synAllStructure & chrStartSynQryLocMap, 
        const string & refChr, 
        const int & qrySynStartTmp, 
        const int & qrySynEndTmp, 
        const int & refStart, 
        const int & refEnd, 
        const int & qryStart, 
        const int & qryEnd, 
        const int & synStart, 
        const int & synEnd, 
        int & qrySynStart, 
        int & qrySynEnd, 
        bool & whileBool, 
        int & aliRowStartNum, 
        int & aliRowEndNum
    );


    /**
     * @brief 更新临时的syn坐标
     * 
     * @param sampleName            样品名，用于判断该行是否要进行判断
     * @param synLocSampleVecMap    map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
     * @param synChr                染色体
     * @param synIdx                syn在该染色体的vector中索引
     * @param synStart              syn的起始
     * @param synEnd                syn的终止
     * @param qrySynStart           qry上syn的起始
     * @param qrySynEnd             qry上syn的终止
     * @param whileBool             判断是否跳出while循环
     * 
     * @return 0
    */
    int renew_syn_loc(
        const string& sampleName, 
        const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap, 
        const string& synChr, 
        int& synIdx, 
        int& synStart, 
        int& synEnd, 
        int& qrySynStart, 
        int& qrySynEnd, 
        bool& whileBool
    );


    /**
     * @brief 构件所有sample共有的syn坐标索引
     * 
     * @param inputFileName         syn_loc输出文件
     * 
     * @return synLocSampleVecMap   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    */
    map<string, vector<tuple<int, int, vector<string> > > > build_syn_idx(const string & inputFileName);


    /**
     * @brief 得到syn在qry上的坐标
     * 
     * @param sampleName                 样品名
     * @param inputFileName              show-aligns 输出文件
     * @param synLocSampleVecMap         map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
     * @param findRevBool                是否遍历反向比对序列
     * 
     * @return chrStartSynQryLocMap      synAllStructure
    */
    synAllStructure get_syn_coor(
        string sampleName, 
        const string& inputFileName, 
        const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap, 
        const bool& findRevBool
    );


    /**
     * @brief 保存结果
     * 
     * @param synLocSampleVecMap            build_syn_idx构建的索引   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
     * @param sampleChrStartSynQryLocMap    get_syn_coor返回值   map<sampleName, chrStartSynQryLocMap>   map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
     * @param outputFileName                输出文件名
     * 
     * @return chrStartSynQryLocMap      synAllStructure
    */
    int save_result(
        const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap,
        const map<string, unordered_map<string, unordered_map<int, tuple<int, int> > > >& sampleChrStartSynQryLocMap,
        const string& outputFileName
    );
}

#endif