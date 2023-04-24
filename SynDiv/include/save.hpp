#ifndef save_hpp
#define save_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>
#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

namespace SAVE
{
    /**
     * @brief ��vcf�ļ�
     * 
     * @param outputFileName
     * 
     * @return
    **/
    class SAVE
    {
    private:
        // vcf�ļ�
        string outputFileName_;

        // ����ļ���
        ofstream fpO;
        gzFile gzfpO;
    public:
        SAVE() {}
        SAVE(string aliFileName) {
            outputFileName_ = aliFileName;

            if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos)
            {
                // ����ļ���
                gzfpO = gzopen(outputFileName_.c_str(), "wb");
                if(!gzfpO)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "'" << outputFileName_ << "': No such file or directory." << endl;
                    exit(1);
                }
            }
            else if (outputFileName_.size() > 0)
            {
                // ����ļ���
                fpO.open(outputFileName_, ios::out);
                if(!fpO)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "'" << outputFileName_ << "': No such file or directory." << endl;
                    exit(1);
                }
            }
        }
        ~SAVE() {
            if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos)
            {
                // �ر��ļ�
                gzclose(gzfpO);
            }
            else if (outputFileName_.size() > 0)
            {
                // �ر��ļ�
                fpO.close();
            }
        }
        int save(string & outTxt_);
    };


    /**
     * @brief �洢
     * 
     * @return int
    **/
    int SAVE::save(
        string & outTxt_
    )
    {
        outTxt_ = strip(outTxt_, '\n');  // ȥ�����з�

        if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos)
        {
            outTxt_ += "\n";
            gzwrite(gzfpO, outTxt_.c_str(), outTxt_.length());
        }
        else if (outputFileName_.size() > 0)
        {
            fpO << outTxt_ << endl;
        }
        else
        {
            cout << outTxt_ << endl;
        }
        
        return 0;
    }
} // namespace SAVE

#endif