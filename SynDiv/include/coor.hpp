#ifndef COOR_HPP
#define COOR_HPP
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <regex>
#include "zlib.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "save.hpp"
#include "GzChunkReader.hpp"

using namespace std;


namespace COOR
{
    // ȫ�ֱ���
    int ratio = 1000;  // ����������ref��qry���ȱ���ֵ


    struct synAllStructure
    {
        string sampleName;
        unordered_map<string, unordered_map<int, tuple<int, int> > > chrStartSynQryLocMap;  // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
    };

    // �Ƿ���Դ���
    bool debug = false;

    /**
     * @brief ����ȫ�ֲ���
     * 
     * @param debugTmp     �Ƿ���Դ���
     * 
     * @return 0
    */
    int change_global(
        bool debugTmp
    )
    {
        debug = debugTmp;

        return 0;
    }


    /**
     * @brief ���� '-- BEGIN alignment ' �ֶ�
     * 
     * @param informationTmp    '+1 278 - 1703'
     * 
     * @return tuple<strand, start, end>
    */
    tuple<string, int, int> get_alignment_loc(string informationTmp)
    {
        // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ
        istringstream iss(informationTmp);  // +1 278 - 1703
        string strandTmp;
        int startTmp, endTmp;
        char sep;

        iss >> strandTmp >> startTmp >> sep >> endTmp;

        // ���Ϊ����ȶԣ�������ʼ����ֹ����
        if (strandTmp == "-1") {
            swap(startTmp, endTmp);
        }

        return make_tuple(strandTmp, startTmp, endTmp);
    }


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
    )
    {
        int qrySynLoc = 0;

        if (refStart <= synLoc && synLoc <= refEnd)  // ��������˹���������
        {
            int iIdx = 0;  // ��¼synStart��Ӧ��λ��

            --refStart;  // �����ȼ�1
            for (size_t i = 0; i < refSeq.size(); ++i)  // ѭ��ref����
            {
                if (refSeq[i] != '.')
                {
                    ++refStart;  // ref�������
                }

                if (refStart == synLoc)  // ��syn���������
                {
                    iIdx = i;
                    break;  // �ҵ�������˳���ѭ��
                }
            }

            if (AliQryStrand == "+1")  // qry������ȶ�ʱ
            {
                --qryStart;  // �����ȼ�1
                for (size_t i = 0; i < qrySeq.size(); ++i)
                {
                    if (qrySeq[i] != '.')
                    {
                        ++qryStart;  // ref�������
                    }
                    
                    if (i == iIdx)  // ѭ������λ�ú��˳�
                    {
                        break;  // �ҵ�������˳���ѭ��
                    }
                }
                qrySynLoc = qryStart;
            }
            else  // qry�Ƿ���ȶ�ʱ
            {
                ++qryEnd;  // �����ȼ�1
                for (size_t i = 0; i < qrySeq.size(); ++i)
                {
                    if (qrySeq[i] != '.')
                    {
                        --qryEnd;  // ref����ݼ�
                    }
                    
                    if (i == iIdx)  // ѭ������λ�ú��˳�
                    {
                        break;  // �ҵ�������˳���ѭ��
                    }
                }
                qrySynLoc = qryEnd;
            }
        }

        return qrySynLoc;
    }


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
    )
    {
        whileBool = false;

        if (qrySynStartTmp != 0) {
            if (qrySynStart > 0) {  // ֮ǰ�ҵ���qrySynStart��ѡ����preQrySynStart���������
                if (abs(qrySynStart - synStart) > abs(qrySynStartTmp - synStart))
                {
                    qrySynStart = qrySynStartTmp;
                } else if (debug) {
                    cerr << "skip_start:" << qrySynStartTmp << endl;
                }
            }
            else  // û���ҵ��Ļ�ֱ�Ӹ�ֵ
            {
                qrySynStart = qrySynStartTmp;
            }
            
            aliRowStartNum = aliRowEndNum;  // ����aliRowStartNum������aliRowEndNum����ֹ����
        }
        if (qrySynEndTmp != 0 && qrySynStart != 0)  // ��ʼ�������ҵ�
        {
            if (qrySynEnd > 0)  // ֮ǰ�ҵ���qrySynEnd��ѡ����preQrySynEnd���������
            {
                if (abs(qrySynEnd - synEnd) > abs(qrySynEndTmp - synEnd))
                {
                    qrySynEnd = qrySynEndTmp;
                }
                else if (debug)
                {
                    cerr << "skip_end:" << qrySynEndTmp << endl;
                }
            }
            else  // û���ҵ��Ļ�ֱ�Ӹ�ֵ
            {
                qrySynEnd = qrySynEndTmp;
            }

            // ����0���ٸ�ֵ���ȶԳ�������С��ratio���ٸ�ֵ
            int refSynLen = abs(synEnd - synStart);
            int qrySynLen = abs(qrySynEnd - qrySynStart);
            if (qrySynEnd > 0 && max(refSynLen, qrySynLen)/(float)min(refSynLen, qrySynLen) < ratio) {
                chrStartSynQryLocMap.chrStartSynQryLocMap[refChr][synStart] = make_tuple(
                    min(qrySynStart, qrySynEnd), 
                    max(qrySynStart, qrySynEnd)
                );
            }
        }

        // debug
        if (qrySynStartTmp != 0 || qrySynEndTmp != 0)
        {
            if (debug)
            {
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
    )
    {
        map<string, vector<tuple<int, int, vector<string> > > >::const_iterator findIter = synLocSampleVecMap.find(synChr);  //  �ڹ�����map����Ⱦɫ��
        if (findIter != synLocSampleVecMap.end())  // �ҵ���
        {
            auto& synLocSampleVec = findIter->second;  // vector<tuple<int, int, vector<string> > >

            if (synIdx < synLocSampleVec.size())  // ��ֹԽ��
            {
                // ����ù�����û�ж�Ӧ����Ʒ����������һ������
                while (find(get<2>(synLocSampleVec[synIdx]).begin(), get<2>(synLocSampleVec[synIdx]).end(), sampleName) == get<2>(synLocSampleVec[synIdx]).end())
                {
                    synIdx++;  // ��������
                    // ���Խ�磬�������겢��������
                    if (synIdx >= synLocSampleVec.size())
                    {
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
                synIdx++;  // ��������
                qrySynStart = 0;
                qrySynEnd = 0;
                whileBool = true;
            }
            else
            {
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
     * @brief ��������sample���е�syn��������
     * 
     * @param inputFileName         syn_loc����ļ�
     * 
     * @return synLocSampleVecMap   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
    */
    map<string, vector<tuple<int, int, vector<string> > > > build_syn_idx(const string & inputFileName)
    {
        // ���湲��������
        map<string, vector<tuple<int, int, vector<string> > > > synLocSampleVecMap;  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >

        // open syntenic coordinations file
        GzChunkReader GzChunkReaderClass(inputFileName);

        // read line
        string line;
        while (GzChunkReaderClass.read_line(line))
        {
            // ��������
            if (line.empty())
            {
                continue;
            }

            // ���
            std::istringstream iss(line);
            vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

            synLocSampleVecMap[lineVec[0]].push_back(make_tuple(stoi(lineVec[1]), stoi(lineVec[2]), split(lineVec[4], ",")));  // �洢��������ref�ϵ������Լ���λ�����Ʒ����
        }

        return synLocSampleVecMap;
    }


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
    )
    {
        // ��ȡ��������qry�ϵ�����
        // map<chr, map<synStart, tuple<qrySynStart, qrySynEnd> > >
        synAllStructure chrStartSynQryLocMap;

        // ���û���ύ��Ʒ���������ļ���ȡ
        if (sampleName.empty())
        {
            vector<string> inputFileNameVecTmp = split(inputFileName, "/");  // ·�����
            sampleName = inputFileNameVecTmp[inputFileNameVecTmp.size() - 1];
            std::regex reg1(".aligns");  // �滻
            sampleName = regex_replace(sampleName, reg1, "");  // 'An-1.aligns' ɾ�� '.aligns'
        }

        // ��¼��Ʒ����
        chrStartSynQryLocMap.sampleName = sampleName;

        // �ȳ�ʼ��sylLocVecOutMap
        for (const auto& [chromosome, locationSampleVec] : synLocSampleVecMap)  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
        {
            for (const auto& [synStart, synEnd, sampleVec] : locationSampleVec)  // vector<tuple<refStart, refEnd, vector<sample> > >
            {
                chrStartSynQryLocMap.chrStartSynQryLocMap[chromosome].emplace(synStart, make_tuple(0, 0));
            }
        }
        
        // ��ʱ��syn����
        string synChr = "";
        int synIdx = 0;  // ������ȡ�����Ե�����
        int synStart = 0;
        int synEnd = 0;
        // �洢qry��syn������
        int qrySynStart;
        int qrySynEnd;

        // ��ʱ��Ⱦɫ��ź�����
        string refChr = "";
        string qryChr = "";
        // alignmenmt����Ϣ
        string AliRefStrand = "";
        string AliQryStrand = "";
        int aliRefStart = 0;
        int aliRefEnd = 0;
        int aliQryStart = 0;
        int aliQryEnd = 0;
        // �����ж��Ƿ�ѭ������alignment�Ĳ���ֵ  (false)
        bool aliBool = false;
        // ���ڼ�¼syn�������ڸ�alignment�ĵڼ���
        int aliRowStartNum = 0;
        int aliRowEndNum = INT32_MAX;

        // �����е���Ϣ
        int refStart = 0;
        int refEnd = 0;
        string refSeq = "";
        int qryStart = 0;
        int qryEnd = 0;
        string qrySeq = "";

        // �������ʱboolֵ����ʹ�ã�ֻ����Ϊ�����ύ
        bool whileBoolTmp = true;

        // open aligns file
        GzChunkReader GzChunkReaderClass(inputFileName);

        // ��¼ѭ��������  ż��->ref  ����->qry
        int forIdx = -1;

        // read line
        string line;
        
        while (GzChunkReaderClass.read_line(line))
        {
            // ��������
            if (line.empty())
            {
                continue;
            }
            
            if (line == "\n" || 
                line[0] == ' ' || 
                line[0] == '=' || 
                line.find("--   END alignment ") != string::npos)  // �ո�/=/\n/END alignment��ֱ������
            {
                continue;
            }
            else  // ��ȡsown-aligns����
            {
                if (line.find("-- Alignments between") != string::npos)  // �µıȶ�Ⱦɫ��    -- Alignments between chr1 and chr1
                {
                    // ����ַ���
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    refChr = lineVec[3];  // refȾɫ�����ʱֵ
                    qryChr = lineVec[5];  // qryȾɫ�����ʱֵ
                    
                    // ��������Ϣ����
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

                    // ��������
                    AliRefStrand = "";
                    AliQryStrand = "";

                    // ����ַ���
                    line.clear();
                    string().swap(line);
                    continue;
                }
                else if (line.find("-- BEGIN alignment ") != string::npos)  // �µıȶ�����    -- BEGIN alignment [ +1 1078 - 68996 | +1 1 - 67961 ]
                {
                    line = strip(line.erase(0, 21), '\n');  // ɾ�� '-- BEGIN alignment [ ' ��21���ַ�
                    line = line.erase(line.size()-2, 2);  // ɾ�� ' ]' ��2���ַ���λ�����������

                    vector<string> lineVec = split(line, " | ");  // �ָ� '+1 278 - 1703 -1 2148 - 751'

                    // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ  ref
                    string AliRefStrandTmp;
                    int aliRefStartTmp;
                    int aliRefEndTmp;
                    tie(AliRefStrandTmp, aliRefStartTmp, aliRefEndTmp) = get_alignment_loc(
                        lineVec[0]
                    );

                    // ��ȡ�ȶԷ�����ʼ����ֹ����Ϣ  qry
                    string AliQryStrandTmp;
                    int aliQryStartTmp;
                    int aliQryEndTmp;
                    tie(AliQryStrandTmp, aliQryStartTmp, aliQryEndTmp) = get_alignment_loc(
                        lineVec[1]
                    );

                    // �������е�����
                    forIdx = -1;

                    // �����ڹ�����֮ǰ����������read��
                    if (aliRefEndTmp < synStart)
                    {
                        aliBool = false;
                    }
                    // �ж��Ƿ��������򻥲����У���������read��
                    else if (!findRevBool && (AliRefStrandTmp == "-1" || AliQryStrandTmp == "-1"))
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Warning: Reverse alignment skipped -> " << synChr << " " << line << endl;
                        aliBool = false;
                    }
                    // ref��qry��Ⱦɫ�岻һ������������read��
                    else if (refChr != qryChr)
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Warning: Skipped alignment due to chromosomal mismatch -> " << refChr << " != " << qryChr << endl;
                        aliBool = false;
                    }
                    else
                    {
                        // ��������
                        AliRefStrand = AliRefStrandTmp;
                        aliRefStart = aliRefStartTmp;
                        aliRefEnd = aliRefEndTmp;
                        AliQryStrand = AliQryStrandTmp;
                        aliQryStart = aliQryStartTmp;
                        aliQryEnd = aliQryEndTmp;
                        aliBool = true;

                        // ����syn���ڵĴ����У�
                        aliRowStartNum = ((synStart - aliRefStart)/49)*2-10;
                        aliRowEndNum = ((synEnd - aliRefStart)/49)*2-10;
                    }
                }
                // ֻѭ�������ֿ�ͷ�Ͱ���syn����
                // ���map�������ˣ�ͬ�������������������  'synEnd==0'
                else if (isdigit(line[0]) != 0 && aliBool && synEnd > 0)   
                {
                    forIdx++;  // ��¼����  ż��ref ����qry

                    // ���С����ֵ��������
                    if (forIdx < aliRowStartNum)
                    {
                        continue;
                    }

                    // ����ַ���
                    std::istringstream iss(line);
                    vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                    if (( forIdx & 1 ) == 0)  // ż��
                    {
                        // �������� (this)
                        refSeq = lineVec[1];
                        // ��ʱ��ref����
                        string refSeqTmp = refSeq;
                        // ɾ��  '.'  ���㳤��
                        refSeqTmp.erase(remove(refSeqTmp.begin(), refSeqTmp.end(), '.'), refSeqTmp.end());
                        if (AliRefStrand == "+1")  // ����ȶ�
                        {
                            refStart = stoi(lineVec[0]);
                            refEnd = refStart + refSeqTmp.size() - 1;
                        }
                        else  // ����ȶ�
                        {
                            refEnd = stoi(lineVec[0]);
                            refStart = refEnd - refSeqTmp.size() + 1;
                        }
                    }
                    else
                    {
                        // �������� (this)
                        qrySeq = lineVec[1];
                        // ��ʱ��ref����
                        string qrySeqTmp = qrySeq;
                        // ɾ��  '.'  ���㳤��
                        qrySeqTmp.erase(remove(qrySeqTmp.begin(), qrySeqTmp.end(), '.'), qrySeqTmp.end());
                        if (AliQryStrand == "+1")  // ����ȶ�
                        {
                            qryStart = stoi(lineVec[0]);
                            qryEnd = qryStart + qrySeqTmp.size() - 1;
                        }
                        else  // ����ȶ�
                        {
                            qryEnd = stoi(lineVec[0]);
                            qryStart = qryEnd - qrySeqTmp.size() + 1;
                        }
                        
                        // ��������������ڸ�����֮ǰ�ˣ����¹����Ե����꣬����Ҫ��ֹ��ѭ�������Ե�'synEnd==0'ʱ�����������ֵ��ˣ�����
                        while (synEnd > 0 && synEnd < refStart)
                        {
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

                            // ����syn���ڵĴ����У�
                            aliRowStartNum = ((synStart - aliRefStart)/49)*2-10;
                            aliRowEndNum = ((synEnd - aliRefStart)/49)*2-10;

                            // �жϸ����ȶ��Ƿ������syn���������±ߵ��ж�����
                            if (aliRefEnd < synStart)  // �����ڹ�����֮ǰ����������read
                            {
                                aliBool = false;
                            }
                        }

                        // �ж��Ƿ����whileѭ��
                        bool whileBool = true;

                        // �����꣬push����ͼ��
                        while (((refStart <= synStart && synStart <= refEnd) || 
                                (refStart <= synEnd && synEnd <= refEnd)) &&
                                whileBool)
                        {
                            int qrySynStartTmp = find_qry_syn(
                                synStart, 
                                refStart, 
                                refEnd, 
                                refSeq, 
                                AliQryStrand, 
                                qryStart, 
                                qryEnd, 
                                qrySeq
                            );

                            int qrySynEndTmp = find_qry_syn(
                                synEnd, 
                                refStart, 
                                refEnd, 
                                refSeq, 
                                AliQryStrand, 
                                qryStart, 
                                qryEnd, 
                                qrySeq
                            );

                            // ���Ƿ��ҵ����ҵ������ܹ�ϣ��������
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

                            // ����������궼�ҵ��������Ѿ�����synEnd�����У�����syn������
                            if ((qrySynStart > 0 && qrySynEnd > 0) || 
                                (refStart <= synEnd && synEnd <= refEnd) || 
                                synEnd < refStart)
                            {
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
                                // �жϸ����ȶ��Ƿ������syn���������±ߵ��ж�����
                                if (aliRefEnd < synStart)  // �����ڹ�����֮ǰ����������read
                                {
                                    aliBool = false;
                                }
                                // ����syn���ڵĴ����У�
                                aliRowStartNum = ((synStart - aliRefStart)/49)*2-10;
                                aliRowEndNum = ((synEnd - aliRefStart)/49)*2-10;
                            }
                        }
                    }
                }
            }
        }

        return chrStartSynQryLocMap;
    }


    /**
     * @brief ������
     * 
     * @param synLocSampleVecMap       build_syn_idx����������   map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
     * @param synLocVecOutMapTmpVec    get_syn_coor����ֵ   vector<synAllStructure>
     * @param outputFileName             ����ļ���
     * 
     * @return chrStartSynQryLocMap      synAllStructure
    */
    int save_result(
        const map<string, vector<tuple<int, int, vector<string> > > >& synLocSampleVecMap,
        const vector<COOR::synAllStructure>& synLocVecOutMapTmpVec,
        const string& outputFileName
    )
    {
        // ������
        stringstream outTxtStream;
        for (const auto& [chromosome, locationSampleVec] : synLocSampleVecMap)  // map<chr, vector<tuple<refStart, refEnd, vector<sample> > > >
        {
            // ��¼��һ�����������꣬��ֹ������֮���м��
            uint64_t preSynEnd = 0;

            for (const auto& [synStart, synEnd, sampleVec] : locationSampleVec)  // vector<tuple<refStart, refEnd, vector<sample> > >
            {
                // �����ж��Ƿ��м��
                if(preSynEnd != 0 && synStart - preSynEnd > 1)
                {
                    outTxtStream << chromosome << '\t' << preSynEnd + 1 << '\t' << synStart - 1;
                    for (const auto& it3 : synLocVecOutMapTmpVec)  // vector<synAllStructure>
                    {
                        outTxtStream << '\t' << it3.sampleName << '\t' << 0 << '\t' << 0;
                    }
                    outTxtStream << '\n';
                }

                preSynEnd = synEnd;  // ��¼��һ�������Ե�����

                // ��ǰ�ڵ����Ϣ
                outTxtStream << chromosome << '\t' << synStart << '\t' << synEnd;
                for (const auto& it3 : synLocVecOutMapTmpVec)  // vector<synAllStructure>
                {
                    const auto& chrStartSynQryLoc = it3.chrStartSynQryLocMap.at(chromosome).at(synStart);
                    outTxtStream << '\t' << it3.sampleName << '\t' << get<0>(chrStartSynQryLoc) << '\t' << get<1>(chrStartSynQryLoc);
                }
                outTxtStream << '\n';
            }
        }

        // �������Ļ�򱣴浽�ļ�
        string outTxt = strip(outTxtStream.str(), '\n');
        SAVE::SAVE SaveClass(outputFileName);
        SaveClass.save(outTxt);

        return 0;
    }
}

#endif