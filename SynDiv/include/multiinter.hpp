#ifndef MULTIINTER_hpp
#define MULTIINTER_hpp
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include "save.hpp"
#include "GzChunkReader.hpp"

using namespace std;

namespace MULTIINTER
{
    bool debug;
    // bool debug = false;
    // bool debug = true;

    pair<string, map<string, vector<pair<int, int> > > > build_syn_index(
        const string & inputFileName
    );

    tuple<int, int, vector<string> > loc_find(
        const vector<pair<int, int > > synVec, 
        const map<int, string> & idxLineMap
    );

    int save_result(
        const map<string, map<int, vector<tuple<int, vector<string> > > > > & outChrStartEndLineVecMap, 
        map<int, string> idxLineMap, 
        const string & outputName
    );



    /**
     * @brief build syn index for syri.out
     * 
     * @param inputFileName   the output of syri
     * 
     * @return pair<lineName, chrSynVecMap>        pair<lineName, map<chromosome, vector<pair<refStart, refEnd>>>>
    **/
    pair<string, map<string, vector<pair<int, int> > > > build_syn_index(
        const string & inputFileName
    )
    {
        // Ʒ�ֵ�����
        string lineName = inputFileName;

        vector<string> inputFileNameVecTmp = split(inputFileName, "/");  // ·�����
        lineName = inputFileNameVecTmp[inputFileNameVecTmp.size() - 1];  // ���һ��

        string prefix = ".syri.out";
        auto iter1 = lineName.find(prefix);
        if(iter1 != string::npos)
        {
            lineName.replace(lineName.begin() + iter1, lineName.begin() + iter1 + prefix.size(), "");
        }
        
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Building index: " << inputFileName << ".\n";  // print log

        map<string, vector<pair<int, int> > > chrSynVecMap;  // map<chromosome, vector<pair<refStart, refEnd>>>

        // open file
        GzChunkReader GzChunkReaderClass(inputFileName);

        // ��һ�������Ե�����
        string preChromosome;
        int PreRefStart = 0;
        int PreRefEnd = 0;

        // read line
        string line;

        while (GzChunkReaderClass.read_line(line))
        {
            // ��������
            if (line.empty() || line.find("SYNAL") == string::npos)
            {
                continue;
            }

            // ���
            std::istringstream iss(line);
            vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

            string chromosome = lineVec[0];
            int refStart = stoi(lineVec[1]);
            int refEnd = stoi(lineVec[2]);

            if(preChromosome != chromosome && preChromosome.size() > 0)  // �����Ⱦɫ����/��һ��Ⱦɫ���Թ�������
            {
                chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));
                PreRefStart = refStart;
                PreRefEnd = refEnd;
            }
            else
            {
                if(PreRefStart == 0)  // ��һ��Ⱦɫ���һ��syn
                {
                    PreRefStart = refStart;
                    PreRefEnd = refEnd;
                }
                else if(refStart <= PreRefEnd)  // ��һ������������
                {
                    if(refEnd > PreRefEnd)  // ���ְ���
                    {
                        PreRefStart = PreRefStart;
                        PreRefEnd = refEnd;
                    }
                    else  // ��ȫ����
                    {
                        PreRefStart = PreRefStart;
                        PreRefEnd = PreRefEnd;
                    }
                }
                else  // ���ص�
                {
                    chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));
                    PreRefStart = refStart;
                    PreRefEnd = refEnd;
                }
            }

            preChromosome = chromosome;
        }

        // ���һ���������
        chrSynVecMap[preChromosome].push_back(make_pair(PreRefStart, PreRefEnd));

        return make_pair(lineName, chrSynVecMap);
    }


    /**
        * @brief �Һϼ�
        * 
        * @param chrLineSynVecMap          map<string, map<string, vector<pair<int, int> > > >, map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
        * 
        * @return pair(outChrStartEndLineVecMap, idxLineMap)   pair(map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>, �洢line������������)
    **/
    pair<map<string, map<int, vector<tuple<int, vector<string> > > > >, map<int, string> > syn_multiinter_find(
        const map<string, map<string, vector<pair<int, int> > > > & chrLineSynVecMap
    )
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Searching." << endl;
        // ������
        map<string, map<int, vector<tuple<int, vector<string> > > > > outChrStartEndLineVecMap;  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>

        map<int, string> idxLineMap;  // �洢line������������

        // ����Ⱦɫ��
        for(auto it1 : chrLineSynVecMap)  // map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
        {
            // it1.first -> chromosome
            string chromosome = it1.first;

            // it1.second -> map<lineName, vector<pair<refStart, refEnd>>>

            int lineIdx = 0;  // lineName������
            
            // ��Ⱦɫ���Ӧ�����й������б��������Ϣ��ȡ����������ֵ��һ������Ϊsyn�б�ĵ�һ��
            map<string, tuple<vector<pair<int, int> >, int, pair<int, int> > > LineSynVecIdxSynMap;  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
            for(auto it2 : it1.second)  // map<lineName, vector<pair<refStart, refEnd>>>
            {
                // it2.first -> lineName
                // it2.second -> vector<pair<refStart, refEnd>>

                LineSynVecIdxSynMap[it2.first] = make_tuple(it2.second, 0, it2.second[0]);

                idxLineMap[lineIdx] = it2.first;  // �洢line��Ϣ
                lineIdx++;  // ������1
            }


            // Ѱ�ң��ȹ���һ����ʱ�Ĺ������б��洢ÿ��������һ������������
            vector<pair<int, int> > synVec;  // �洢���������䣬�Һϼ�
            
            int lineNum = idxLineMap.size();  // sample�����������ж��Ƿ���ֹѭ��
            int endNum = 0;  // �����ж��Ƿ���ֹwhileѭ��

            for(auto it1 : LineSynVecIdxSynMap)  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
            {
                // it1.first -> lineName
                string lineName = it1.first;

                // it1.second -> tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> >
                synVec.push_back(get<2>(it1.second));

                // ������Դ���Ļ�����ӡlog
                if (debug)
                {
                    cerr << chromosome << " " << lineName << ":" << get<2>(it1.second).first << "-" << get<2>(it1.second).second << endl;
                }
                
                if(get<1>(it1.second) == get<0>(it1.second).size())  // �ж��Ƿ�ѭ�����
                {
                    endNum++;
                }
            }
            
            while(endNum < lineNum)  // �� endNum == lineNum ʱ����ѭ��
            {
                endNum = 0;  // ���ü���

                // �ҽ���
                int outRefStart = 0;
                int outRefEnd = 0;
                vector<string> outLineVec;
                tie(outRefStart, outRefEnd, outLineVec) = loc_find(
                    synVec, 
                    idxLineMap
                );

                // ������Դ���Ļ�����ӡlog
                if (debug)
                {
                    cerr << chromosome << " loc:" << outRefStart << "-" << outRefEnd << endl << chromosome << " outLine:";

                    for (auto outLine : outLineVec)
                    {
                        cerr << " " << outLine;
                    }
                    cerr << endl << endl;
                }

                // �洢�������
                outChrStartEndLineVecMap[chromosome][outRefStart].push_back(make_tuple(outRefEnd, outLineVec));

                // �������
                synVec.clear();
                vector<pair<int, int > >().swap(synVec);

                // ���¹���������������
                for(auto it1 : LineSynVecIdxSynMap)  // map<lineName, tuple<vector<pair<refStart, refEnd> >, index, pair<refStart, refEnd> > >
                {
                    // it1.first -> lineName
                    string lineName = it1.first;

                    // Ʒ�ֶ�Ӧ�Ĺ��������꣬�����ж��Ƿ��������
                    vector<pair<int, int> > lineSynVec = get<0>(it1.second);
                    int lineSynIdx = get<1>(it1.second);
                    pair<int, int> lineSycLocTmp = get<2>(it1.second);

                    // ��������
                    if(lineSycLocTmp.first >= outRefEnd)  // ���û�õ��ù��������䣬ֱ����һ��Ʒ��
                    {
                        synVec.push_back(lineSycLocTmp);  // �������꣬������һ��������
                    }
                    else if(outRefEnd >= lineSycLocTmp.second)  // ���������걻������
                    {
                        lineSynIdx++;
                        if(lineSynIdx >= lineSynVec.size())  // �ж��Ƿ�ѭ�����
                        {
                            lineSycLocTmp = make_pair(0, 0);  // ȫ����
                            endNum++;
                        }
                        else
                        {
                            lineSycLocTmp = lineSynVec[lineSynIdx];  // ������ָ����һ������
                        }

                        synVec.push_back(lineSycLocTmp);
                        LineSynVecIdxSynMap[lineName] = make_tuple(lineSynVec, lineSynIdx, lineSycLocTmp);  // ���¹�ϣ�������
                    }
                    else  // ���������걻����һ��
                    {
                        lineSycLocTmp.first = outRefEnd;

                        synVec.push_back(lineSycLocTmp);  // ��������
                        LineSynVecIdxSynMap[lineName] = make_tuple(lineSynVec, lineSynIdx, lineSycLocTmp);  // ���¹�ϣ�������
                    }

                    // ������Դ���Ļ�����ӡlog
                    if (debug)
                    {
                        cerr << chromosome << " " << lineName << ":" << lineSycLocTmp.first << "-" << lineSycLocTmp.second << " endNum:" << endNum << endl;
                    }
                }
            }
        }

        return make_pair(outChrStartEndLineVecMap, idxLineMap);
    }


    /**
        * @brief �Ҽ���_find
        * 
        * @param synVec      vector<pair<int, int > >, vector<pair<refStart, refEnd> >
        * 
        * @return tuple<int, int, vector<string> >  make_tuple(outStart, outEnd, lineName)
    **/
    tuple<int, int, vector<string> > loc_find(
        const vector<pair<int, int > > synVec, 
        const map<int, string> & idxLineMap
    )
    {
        // ����������
        int outStart = INT_MAX;
        int outEnd = INT_MAX;

        vector<string> lineNameVec;  // �����������Ӧ��Ʒ������
        

        for(auto it1 : synVec)  // �ҹ�������С�����꣬�϶�����ʼ������
        {
            if(it1.first == 0 && it1.second == 0) continue;  // �����Ϊ0�������lineѭ�����

            outStart = min(outStart, it1.first);  // ���й������е���Сֵ���϶�����ʼ������
        }

        for(auto it1 : synVec)  // �ҹ����Եڶ�С������
        {
            if(it1.first == 0 && it1.second == 0) continue;  // �����Ϊ0�������lineѭ�����

            if (it1.first != outStart)  // ��������Сֵ
            {
                outEnd = min(outEnd, it1.first);
            }
            outEnd = min(outEnd, it1.second);
        }

        if (outEnd == INT_MAX)  // ��ֹ������ֹ����ʼһ�������
        {
            outEnd = outStart;
        }
        

        int lineIdx = 0;
        for(auto it1 : synVec)  // ����С�����Ӧ��Ʒ����Ϣ
        {
            if(it1.first == 0 && it1.second == 0)  // �����Ϊ0�������lineѭ�����
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
        * @brief �洢���
        * 
        * @param outChrStartEndLineVecMap      const map<chr, map<refStart, vector<tuple<refEnd, vector<lineName> > > > >
        * @param idxLineMap                  map<index, lineName>
        * @param outputName                  ����ļ���
        * 
        * @return tuple<int, int, vector<string> >  make_tuple(outStart, outEnd, lineName)
    **/
    int save_result(
        const map<string, map<int, vector<tuple<int, vector<string> > > > > & outChrStartEndLineVecMap, 
        map<int, string> idxLineMap, 
        const string & outputName
    )
    {
        SAVE::SAVE SAVEClass(outputName);

        stringstream outStream; // ʹ�� stringstream �����ַ���ƴ��
        static const uint64_t CACHE_SIZE = 1024 * 1024 * 10; // �����СΪ 10mb
        outStream.str().reserve(CACHE_SIZE);

        // ��lineName��ӵ��б���
        vector<string> lineNameVec;
        for(auto it : idxLineMap)
        {
            lineNameVec.push_back(it.second); 
        }

        // ������
        for(auto it1 : outChrStartEndLineVecMap)  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName> > > > >
        {
            for(auto it2 : it1.second)  // map<refStart, vector<tuple<refEnd, vector<lineName> > > >
            {
                for(auto it3 : it2.second)  // vector<tuple<refEnd, vector<lineName> > >
                {
                    outStream << it1.first << "\t" << to_string(it2.first) << "\t" + to_string(get<0>(it3)) << "\t" << to_string(get<1>(it3).size()) << "\t";

                    // �����Ժ��е�Ʒ����������
                    int idxTmp=0;

                    // ���(�����м�֮�������)
                    vector<string> lineNameVecTmp;
                    vector<int> boolVecTmp;
                    for (size_t i = 0; i < lineNameVec.size(); i++)  // vector<lineName>
                    {
                        if(lineNameVec[i] == get<1>(it3)[idxTmp])  // ���ж�Ӧ��Ʒ��
                        {
                            lineNameVecTmp.push_back(lineNameVec[i]);
                            boolVecTmp.push_back(1);
                            idxTmp++;
                        }
                        else
                        {
                            boolVecTmp.push_back(0);
                        }

                        // �����vector�Ľ�β������forѭ��
                        if (idxTmp == get<1>(it3).size())
                        {
                            break;
                        }  
                    }

                    // ������Ϊ�գ�������λ��
                    if (lineNameVecTmp.size() == 0)
                    {
                        continue;
                    }
                    
                    outStream << join(lineNameVecTmp, ",") << "\t" + join(boolVecTmp, "\t") << "\n";

                    if (outStream.tellp() >= CACHE_SIZE)  // �����СΪ 10mb
                    {
                        string outTxt = outStream.str();
                        SAVEClass.save(outTxt);
                        // ��� stringstream
                        outStream.str(string());
                        outStream.clear();
                    }
                }
            }
        }

        if (outStream.tellp() >= 0)  // ���дһ��
        {
            string outTxt = outStream.str();
            SAVEClass.save(outTxt);
            // ��� stringstream
            outStream.str(string());
            outStream.clear();
        }
        
        return 0;
    }
}

#endif