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

// ȫ�ֱ���
extern int thresholdLength;  // ����������ref��qry���ȱ���ֵ

// �Ƿ���Դ���
extern bool debugCoor;

void help_coor(char* argv[]);
int main_coor(int argc, char* argv[]);

namespace COOR
{
    struct synAllStructure
    {
        string sampleName;
        unordered_map<string, unordered_map<int64_t, tuple<int64_t, int64_t> > > chrStartSynQryLocMap;  // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
    };


    /**
     * @brief ���� '-- BEGIN alignment ' �ֶ�
     * 
     * @param informationTmp    '+1 278 - 1703'
     * 
     * @return tuple<strand, start, end>
    */
    tuple<string, int64_t, int64_t> get_alignment_loc(string informationTmp);


    /**
     * @brief ��qry�Ϲ����Ե�����
     * 
     * @param synLoc        �����Ե�����
     * @param refStart      ref����ʼ
     * @param refEnd        ref����ֹ
     * @param refSeq        ref������
     * @param qryStrand     qry�ȶԷ���
     * @param qryStart      qry����ʼ
     * @param qryEnd        qry����ֹ
     * @param qrySeq        qry������
     * 
     * @return qrySynLoc    qry��syn�����꣬ 0-û�ҵ�
    */
    int64_t find_qry_syn(
        const int64_t & synLoc, 
        int64_t refStart, 
        int64_t refEnd, 
        const string & refSeq, 
        const string & AliQryStrand, 
        int64_t qryStart, 
        int64_t qryEnd, 
        const string & qrySeq
    );


    /**
     * @brief ��qry�Ϲ����Ե�����
     * 
     * @param chrStartSynQryLocMap  synAllStructure
     * @param refChr                Ⱦɫ���
     * @param qrySynStartTmp        find_qry_syn�ҵ�qrySynStart
     * @param qrySynEndTmp          find_qry_syn�ҵ�qrySynEnd
     * @param refStart              �ñȶ���ref����ʼ
     * @param refEnd                �ñȶ���ref����ֹ
     * @param qryStart              �ñȶ���qry����ʼ
     * @param qryEnd                ���ȶ���qry����ֹ
     * @param synStart              ĿǰҪ�ҵ�syn��ref��ʼ
     * @param synEnd                ĿǰҪ�ҵ�syn��ref��ֹ
     * @param qrySynStart           syn��qry��ʼ
     * @param qrySynEnd             syn��qry��ֹ
     * @param whileBool             �ж��Ƿ���Ҫwhileѭ��
     * @param aliRowStartNum        syn��ʼλ���ڸ�alignment�Ĵ�������
     * @param aliRowEndNum          syn��ֹλ���ڸ�alignment�Ĵ�������
     * 
     * @return qrySynLoc    qry��syn�����꣬ 0-û�ҵ�
    */
    int64_t syn_all_loc_push(
        synAllStructure & chrStartSynQryLocMap, 
        const string & refChr, 
        const int64_t & qrySynStartTmp, 
        const int64_t & qrySynEndTmp, 
        const int64_t & refStart, 
        const int64_t & refEnd, 
        const int64_t & qryStart, 
        const int64_t & qryEnd, 
        const int64_t & synStart, 
        const int64_t & synEnd, 
        int64_t & qrySynStart, 
        int64_t & qrySynEnd, 
        bool & whileBool, 
        int64_t & aliRowStartNum, 
        int64_t & aliRowEndNum
    );


    /**
     * @brief ������ʱ��syn����
     * 
     * @param sampleName            ��Ʒ���������жϸ����Ƿ�Ҫ�����ж�
     * @param synLocSampleVecMap    map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
     * @param synChr                Ⱦɫ��
     * @param synIdx                syn�ڸ�Ⱦɫ���vector������
     * @param synStart              syn����ʼ
     * @param synEnd                syn����ֹ
     * @param qrySynStart           qry��syn����ʼ
     * @param qrySynEnd             qry��syn����ֹ
     * @param whileBool             �ж��Ƿ�����whileѭ��
     * 
     * @return 0
    */
    int renew_syn_loc(
        const string& sampleName, 
        const map<string, vector<tuple<int64_t, int64_t, vector<string> > > >& synLocSampleVecMap, 
        const string& synChr, 
        int64_t& synIdx, 
        int64_t& synStart, 
        int64_t& synEnd, 
        int64_t& qrySynStart, 
        int64_t& qrySynEnd, 
        bool& whileBool
    );


    /**
     * @brief ��������sample���е�syn��������
     * 
     * @param inputFileName         syn_loc����ļ�
     * 
     * @return synLocSampleVecMap   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    */
    map<string, vector<tuple<int64_t, int64_t, vector<string> > > > build_syn_idx(const string & inputFileName);


    /**
     * @brief �õ�syn��qry�ϵ�����
     * 
     * @param sampleName                 ��Ʒ��
     * @param inputFileName              show-aligns ����ļ�
     * @param synLocSampleVecMap         map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
     * @param findRevBool                �Ƿ��������ȶ�����
     * 
     * @return chrStartSynQryLocMap      synAllStructure
    */
    synAllStructure get_syn_coor(
        string sampleName, 
        const string& inputFileName, 
        const map<string, vector<tuple<int64_t, int64_t, vector<string> > > >& synLocSampleVecMap, 
        const bool& findRevBool
    );


    /**
     * @brief ������
     * 
     * @param synLocSampleVecMap            build_syn_idx����������   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
     * @param sampleChrStartSynQryLocMap    get_syn_coor����ֵ   map<sampleName, chrStartSynQryLocMap>   map<sampleName, map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > > >
     * @param outputFileName                ����ļ���
     * 
     * @return chrStartSynQryLocMap      synAllStructure
    */
    int save_result(
        const map<string, vector<tuple<int64_t, int64_t, vector<string> > > >& synLocSampleVecMap,
        const map<string, unordered_map<string, unordered_map<int64_t, tuple<int64_t, int64_t> > > >& sampleChrStartSynQryLocMap,
        const string& outputFileName
    );
}

#endif