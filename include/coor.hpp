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
        unordered_map<string, unordered_map<int, tuple<int, int> > > chrStartSynQryLocMap;  // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
    };


    /**
     * @brief ���� '-- BEGIN alignment ' �ֶ�
     * 
     * @param informationTmp    '+1 278 - 1703'
     * 
     * @return tuple<strand, start, end>
    */
    tuple<string, int, int> get_alignment_loc(string informationTmp);


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
     * @brief ��������sample���е�syn��������
     * 
     * @param inputFileName         syn_loc����ļ�
     * 
     * @return synLocSampleVecMap   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    */
    map<string, vector<tuple<int, int, vector<string> > > > build_syn_idx(const string & inputFileName);


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
        const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap, 
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
        const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap,
        const map<string, unordered_map<string, unordered_map<int, tuple<int, int> > > >& sampleChrStartSynQryLocMap,
        const string& outputFileName
    );
}

#endif