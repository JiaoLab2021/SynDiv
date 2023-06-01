// g++ -c multiinter.cpp -o multiinter -lz -O3
#include "../include/multiinter.hpp"

using namespace std;

 // ���Դ���
bool debugMultiinter = false;

int main_multiinter(int argc, char* argv[])
{
    // �����ļ��б������  syri.out
    vector<string> inputFiles;
    vector<string> inputTitles;

    // �ж��Ƿ��� titlename
    bool haveTitles = false;

    // ����ļ���
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

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running." << endl;

    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index." << endl;
    // �洢����line�Ĺ�������Ϣ
    map<string, map<string, vector<pair<int, int> > > > chrLineSynVecMap;  // map<chromosome, map<lineName, vector<pair<refStart, refEnd>>>>
    
    int indexTmp = 0;  // ��¼ lineName ������
    for(auto it1 : inputFiles)
    {
        string lineName;
        map<string, vector<pair<int, int> > > chrSynVecMap;
        tie(lineName, chrSynVecMap) = MULTIINTER::build_syn_index(it1);

        // ����� inputTitles ���¸�ֵ
        if (haveTitles)
        {
            lineName = inputTitles[indexTmp];
        }
        
        for(auto it2 : chrSynVecMap)
        {
            string chromosome = it2.first;
            chrLineSynVecMap[chromosome][lineName] = it2.second;
        }

        indexTmp++;  // ��������
    }

    /* ************************************ Find intersection ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find intersection." << endl;
    map<string, map<int, vector<tuple<int, vector<string> > > > > outChrStartEndLineVecMap;  // map<chr, map<refStart, vector<tuple<refEnd, vector<lineName>>>>>
    map<int, string> idxLineMap;
    tie(outChrStartEndLineVecMap, idxLineMap) = MULTIINTER::syn_multiinter_find(
        chrLineSynVecMap
    );

    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Saving the Result." << endl;
    MULTIINTER::save_result(
        outChrStartEndLineVecMap, 
        idxLineMap, 
        outputFileName
    );

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done." << endl;
    
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
pair<string, map<string, vector<pair<int, int> > > > MULTIINTER::build_syn_index(
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
pair<map<string, map<int, vector<tuple<int, vector<string> > > > >, map<int, string> > MULTIINTER::syn_multiinter_find(
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
            if (debugMultiinter)
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
            if (debugMultiinter)
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
    * @brief �Ҽ���_find
    * 
    * @param synVec      vector<pair<int, int > >, vector<pair<refStart, refEnd> >
    * 
    * @return tuple<int, int, vector<string> >  make_tuple(outStart, outEnd, lineName)
**/
tuple<int, int, vector<string> > MULTIINTER::loc_find(
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
int MULTIINTER::save_result(
    const map<string, map<int, vector<tuple<int, vector<string> > > > > & outChrStartEndLineVecMap, 
    map<int, string> idxLineMap, 
    const string & outputName
)
{
    SAVE SAVEClass(outputName);

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