// g++ no_syn.cpp -o no_syn -lz -lpthread
#include "../include/no_syn.hpp"

using namespace std;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// ���Դ���
bool debugNoSyn = false;

int main_no_syn(int argc, char* argv[])
{
    // ����������
    string coorFileName;

    // Ⱦɫ�峤���ļ�
    vector<string> lengthsVec;
    vector<string> lengthsTitles;
    // �ж��Ƿ��� titlename
    bool haveTitles = false;

    // ����ļ���
    string outputFileName;

    //Parse command line options
    if(argc <= 2)
    {
        help_no_syn(argv);
        exit(1);
    }

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            help_no_syn(argv);
            exit(1);
        }
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("--coor", 6, parameterLength)) {
            if ((i+1) < argc) {
                coorFileName = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("--lengths", 9, parameterLength)) {
            if ((i+1) < argc) {
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    lengthsVec.push_back(file);
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
                    lengthsTitles.push_back(title);
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
                debugNoSyn = true;
            }
        }
    }

    if (argc <= 2) {
        help_no_syn(argv);
        exit(1);
    }
    if (coorFileName.size() == 0) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Missing file name (--coor)." << endl;
        help_no_syn(argv);
        exit(1);
    }
    if (lengthsVec.size() <= 1) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: Only " << lengthsVec.size() << " length file was specified. Nothing to extract, exiting." << endl;
        help_no_syn(argv);
        exit(1);
    }
    if ((haveTitles == true) && (lengthsVec.size() != lengthsTitles.size())) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: The number of file titles (-n) does not match the number of files (--lengths)." << endl;
        help_no_syn(argv);
        exit(1);
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running." << endl;


   /* ************************************ Build Chromosome Length Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Chromosome Length Index." << endl;
    NOSYN::LENGTH lengthClass(lengthsVec, lengthsTitles);
    lengthClass.index_lengths();

    /* ************************************ Build Syntenic Coordinates Index ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Build Syntenic Coordinates Index." << endl;
    NOSYN::NOSYNCOOR NoSynCoor(coorFileName);
    NoSynCoor.open_coor();

    /* ************************************ Find no-syntenic Coordinates on query genomes ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Find No-syntenic Coordinates on query genomes." << endl;
    NoSynCoor.find_no_syn_coor(
        lengthClass
    );

    /* ************************************ Save The Result ************************************ */
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Save the Result." << endl;
    SAVE SaveClass(outputFileName);
    SaveClass.save(NoSynCoor.outTxt);  // ����

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done." << endl;
    
    return 0;
}


void help_no_syn(char* argv[])
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " --coor FILE --lengths FILE1 FILE2 .. FILEn [options]" << endl
       << "retrieve non-syntenic coordinates" << endl
       << endl
       << "required arguments:" << endl
       << "    --coor              FILE      syntenic coordinates, output file of coor" << endl
       << "    --lengths           FILE      chromosome length files for query genomes (format: chr\tlength), one for multiple mate" << endl
       << endl
       << "optional arguments:" << endl
       << "    -n, --names         STRING    list of names to describe each file in --lengths, one for multiple mate" << endl
       << "    -o, --output        FILE      output coordinates to FILE [stdout]" << endl
       << "    --debug                       debug code (false by default)" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}


/**
 * @brief Ѱ�ҷǹ����Ե�����
 * 
 * @param lengthClass  ��ƷȾɫ�峤��class
 * 
 * @return void
**/
void NOSYN::NOSYNCOOR::find_no_syn_coor(
    LENGTH & lengthClass
)
{   
    // ��������������
    for (auto iter1 = coorChrLociSampleLociMap.begin(); iter1 != coorChrLociSampleLociMap.end(); iter1++)  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >
    {
        string chrTmp = iter1->first;  // Ⱦɫ���

        // �ǹ����Ե��������ж��Ƿ���Ҫ���Ĳο�����������꣬��ֹ�ظ��������
        uint32_t preNoSynNum = 0;

        uint32_t refStart;  // �ο����������ʼ
        uint32_t refEnd;  // �ο����������ֹ


        // �����ж���Ʒ��С�����еķǹ����������Ƿ��ҵ�
        unordered_map<string, bool> sampleNoSynFindBoolMap;  // map<sample, bool>   true->�ҵ���   false->û�ҵ�
        // ��ʼ��ÿ����Ʒ�ķǹ����������Ƿ��ҵ�
        for (auto iter1 : sampleNameVec)  // map<sample, bool>
        {
            if (iter1 == "reference")
            {
                continue;
            }

            sampleNoSynFindBoolMap[iter1] = false;
        }


        for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)  // map<refStart, map<sample, tuple(start, end)> >
        {
            uint32_t refStartTmp;  // �ο����������ʼ
            uint32_t refEndTmp;  // �ο����������ֹ

            // �ڵ����ʼ����ֹλ��
            auto findIter1 = iter2->second.find("reference");
            if (findIter1 != iter2->second.end())
            {
                refStartTmp = get<0>(findIter1->second);
                refEndTmp = get<1>(findIter1->second);
            }
            else
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Error: 'reference' not present in 'coorChrLociSampleLociMap' -> " << iter2->first << endl;
                exit(1);
            }
            
            // �ǹ����Ե��������ж��Ƿ���Ҫ���Ĳο�����������꣬��ֹ�ظ��������
            uint32_t noSynNum = 0;

            // �����һ���ڵ�ȫ�ǹ����Եģ����²ο������������
            if (preNoSynNum  == 0)
            {
                refStart = refStartTmp;
                refEnd = refEndTmp;
            }

            for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++)  // map<sample, tuple(start, end)>
            {
                string sampleTmp = iter3->first;  // ��Ʒ��
                uint32_t StartTmp = get<0>(iter3->second);  // ��ʼ
                uint32_t EndTmp = get<1>(iter3->second);  // ��ֹ

                // �ο��������ֱ������ 
                if (sampleTmp == "reference")
                {
                    continue;
                }

                // ��¼Ϊ�ǹ����Ե�
                if (StartTmp == 0 || EndTmp == 0)
                {
                    // �ǹ�����������1
                    noSynNum++;

                    // �ҷǹ����Ե�λ��
                    // �ȳ�ʼ���ֵ�
                    if (chrStartSampleLociVecMap.find(chrTmp) == chrStartSampleLociVecMap.end())  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
                    {
                        chrStartSampleLociVecMap[chrTmp];
                    }
                    if (chrStartSampleLociVecMap[chrTmp].find(refStart) == chrStartSampleLociVecMap[chrTmp].end())  // map<refStart, map<sample, vector<tuple<start, end> > > >
                    {
                        chrStartSampleLociVecMap[chrTmp][refStart];
                    }

                    // �����sample���걻�ҵ��ˣ�����
                    if (sampleNoSynFindBoolMap[sampleTmp])  // map<sample, findBool>
                    {
                        continue;
                    }

                    // ��ʱ����
                    uint32_t noSynStartTmp = 1;
                    uint32_t noSynEndTmp = 1;

                    auto iter2Tmp = iter2;  // ��ʱ�ڵ�
                    
                    if (iter2 == iter1->second.begin())  // ����ǵ�һ���ڵ�
                    {
                        noSynStartTmp = 1;
                    }
                    else
                    {
                        iter2Tmp = iter2;
                        iter2Tmp--;  // �����ƶ�һλ

                        // �����һ���ڵ㻹��0���Ҳ�Ϊ��һ���ڵ㣬����ǰ��
                        while (iter2Tmp != iter1->second.begin() && get<1>(iter2Tmp->second[sampleTmp]) == 0)
                        {
                            iter2Tmp--;  // �����ƶ�һλ
                        }

                        // ��ֵǰһ������
                        if (get<1>(iter2Tmp->second[sampleTmp]) != 0)
                        {
                            noSynStartTmp = get<1>(iter2Tmp->second[sampleTmp]) + 1;
                        }
                    }

                    
                    iter2Tmp = iter2;  // ��ʼ��������
                    iter2Tmp++;  // �����ƶ�һλ

                    // �����һ���ڵ㻹��0���Ҳ�Ϊ���һ���ڵ㣬��������
                    while (iter2Tmp != iter1->second.end() && get<0>(iter2Tmp->second[sampleTmp]) == 0)
                    {
                        iter2Tmp++;  // �����ƶ�һλ
                    }

                    // ��������һ���ڵ㣬��ΪȾɫ�峤��
                    if (iter2Tmp == iter1->second.end())
                    {
                        // Ⱦɫ�峤��
                        noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);
                    }
                    else
                    {
                        // ��һ���ڵ����ʼ���� -1
                        noSynEndTmp = get<0>(iter2Tmp->second[sampleTmp]) - 1;
                    }

                    // ��ʼ����С����ֹ
                    if (noSynStartTmp < noSynEndTmp)
                    {
                        // ��ʼ���ֵ�
                        if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                        {
                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                        }
                        
                        chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(noSynStartTmp, noSynEndTmp));  // ��¼�ǹ���������

                        sampleNoSynFindBoolMap[sampleTmp] = true;  // ��¼����Ʒ���ҵ���
                    }
                }
                else
                {
                    /* �ж��Ƿ�Ϊ���һ���ڵ� */
                    auto iter2Tmp = iter2;  // ��ʱ�ڵ�
                    iter2Tmp++;  // �����ƶ�һλ

                    // ����ǵ�һ���ڵ㣬�ж��Ƿ�� 1 ��ʼ�ıȶԡ�����Ͳο��������ǹ����ԵĻ���������Ⱦɫ��ǰ�벿�֣�
                    if (iter2 == iter1->second.begin() && StartTmp == 0)
                    {
                        // ��ʼ���ֵ�
                        if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                        {
                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                        }

                        chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(1, StartTmp - 1));  // ��ֵ
                    }
                    else if (iter2Tmp == iter1->second.end())  // ��������һ���ڵ㣬�ж��Ƿ񵽴�Ⱦɫ���β
                    {
                        // Ⱦɫ�峤��
                        uint32_t noSynEndTmp = lengthClass.get_length(sampleTmp, chrTmp);

                        if (EndTmp < noSynEndTmp)
                        {
                            // ��ʼ���ֵ�
                            
                            if (chrStartSampleLociVecMap[chrTmp][refStart].find(sampleTmp) == chrStartSampleLociVecMap[chrTmp][refStart].end())
                            {
                                chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp];
                            }

                            chrStartSampleLociVecMap[chrTmp][refStart][sampleTmp].push_back(make_tuple(EndTmp + 1, noSynEndTmp));  // ��ֵ
                        }
                    }
                    
                    // ���¸���Ʒ��û�ҵ�
                    sampleNoSynFindBoolMap[sampleTmp] = false;
                }
            }
            // ���·ǹ����Ե�����
            preNoSynNum = noSynNum;
        }
    }


    // ���תΪ�ַ���
    for (auto iter1 : chrStartSampleLociVecMap)  // map<chr, map<refStart, map<sample, vector<tuple<start, end> > > > >
    {
        for (auto iter2 : iter1.second)  // map<refStart, map<sample, vector<tuple<start, end> > > >
        {
            // �յĽڵ�����
            if (iter2.second.size() == 0)
            {
                continue;
            }
            
            outTxt += iter1.first + "\t" + to_string(iter2.first);

            for (auto iter3 : iter2.second)  // map<sample, vector<tuple<start, end> > >
            {
                outTxt += "\t" + iter3.first + ":";

                for (size_t i = 0; i < iter3.second.size(); i++)  // vector<tuple<start, end> >
                {
                    if (i > 0)
                    {
                        outTxt += ";";
                    }
                    
                    outTxt += to_string(get<0>(iter3.second[i])) + "-" + to_string(get<1>(iter3.second[i]));
                }

                outTxt += "\t";
            }

            outTxt += "\n";
        }
    }
}