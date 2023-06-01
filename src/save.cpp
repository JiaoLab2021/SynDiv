#include "../include/save.hpp"

using namespace std;


/**
 * @brief save result to file
 * 
 * @return int
**/
int SAVE::save(
    string & outTxt_
)
{
    outTxt_ = strip(outTxt_, '\n');  // remove '\n'

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