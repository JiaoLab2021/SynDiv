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
     * @brief 打开vcf文件
     * 
     * @param outputFileName
     * 
     * @return
    **/
    class SAVE
    {
    private:
        // vcf文件
        string outputFileName_;

        // 输出文件流
        ofstream fpO;
        gzFile gzfpO;
    public:
        SAVE() {}
        SAVE(string aliFileName) {
            outputFileName_ = aliFileName;

            if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos)
            {
                // 输出文件流
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
                // 输出文件流
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
                // 关闭文件
                gzclose(gzfpO);
            }
            else if (outputFileName_.size() > 0)
            {
                // 关闭文件
                fpO.close();
            }
        }
        int save(string & outTxt_);
    };


    /**
     * @brief 存储
     * 
     * @return int
    **/
    int SAVE::save(
        string & outTxt_
    )
    {
        outTxt_ = strip(outTxt_, '\n');  // 去除换行符

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