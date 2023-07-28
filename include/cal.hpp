#ifndef CAL_HPP
#define CAL_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "zlib.h"
#include <map>
#include <tuple>
#include <regex>
#include <memory>
#include <string.h>
#include <malloc.h>
#include <getopt.h>

#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "GzChunkReader.hpp"
#include "save.hpp"
#include "ThreadPool.hpp"


using namespace std;


// 调试代码
extern bool debugCal;

// 线程数
extern int threadsCal;

// 读取文件的缓存大小
extern int32_t readBuffer;

// define parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen))

void help_cal(char* argv[]);
int main_cal(int argc, char* argv[]);


namespace CALNAME
{
    // kseq.h 打开文件
    KSEQ_INIT(gzFile, gzread)

    /**
     * @brief 解析参数
     * 
     * @param alignsTitles     aligns 文件 title
     * @param alignsVec        aligns 路径Vec
     * @param syriConfigFileName   syri 输出的配置文件
     * 
     * @return make_tuple(alignsMap, sampleSampleSyriMap)       map<sampleName, alignsPath>, map<sample1, map<sample2, syriOutPath> >
    */
    tuple<map<string, string>, map<string, map<string, string> > > aligns_parameter(
        const vector<string> & alignsTitles, 
        const vector<string> & alignsVec, 
        const string & syriConfigFileName
    );



    struct CALSTRUCTURE
    {
        // sampleName
        string sampleName;

        // 记录位点共线性的数量
        map<string, map<uint32_t, uint32_t> > sampleSynOutMap;  // map<chr, map<refLoci, synNum> >

        // 记录位点共线性对应的样品名和位置
        map<string, map<uint32_t, string> > sampleSynOutMapTmp;  // map<chr, map<refLoci, sampleName> >  name2:pos1-pos2;name3:pos1-pos3
    };




    // 坐标转换
    class SYNCOOR
    {
    private:
        // 共线性坐标文件  'coor' 输出结果
        string fileName_;

        // 判断 'sampleNameVec' 中是否含有所有的sample
        bool sampleNameVecBool = false;
    public:
        // 共线性坐标
        map<string, map<int, map<string, tuple<uint32_t, uint32_t> > > > coorChrLociSampleLociMap;  // map<refChr, map<refStart, map<sample, tuple(start, end)> > >

        // 样品名
        vector<string> sampleNameVec;

        SYNCOOR() {}
        SYNCOOR(string fileName) {
            fileName_ = fileName;
        }
        ~SYNCOOR() {}

        // 初始化接口，用于继承类使用
        void init(string fileName) {
            fileName_ = fileName;
        }

        // 构建索引
        void open_coor()
        {
            // open file
            GzChunkReader GzChunkReaderClass(fileName_);

            // read line
            string line;
            while (GzChunkReaderClass.read_line(line))
            {
                // 跳过空行
                if (line.empty())
                {
                    continue;
                }

                // 制表符拆分
                std::istringstream iss(line);
                vector<string> lineVec(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

                // 先检查列数是不是三的倍数
                if (lineVec.size() % 3 != 0)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "'Error: The number of columns is not a multiple of three -> " << fileName_ << endl;
                    exit(1);
                }

                // 每三列记录一次
                int idxTmp = lineVec.size() / 3;
                // 临时键 
                string chrTmp = "";
                uint32_t refStartTmp = 0;
                for (size_t i = 0; i < idxTmp; i++)
                {
                    string sampleNameTmp = lineVec[3*i + 0];
                    uint32_t qryStartTmp = 0;
                    uint32_t qryEndTmp = 0;

                    // 判断是否是数字，不是一律为0
                    if (isdigit(lineVec[3*i + 1][0]))
                    {
                        qryStartTmp = stol(lineVec[3*i + 1]);
                    }
                    if (isdigit(lineVec[3*i + 2][0]))
                    {
                        qryEndTmp = stol(lineVec[3*i + 2]);
                    }
                    
                    // 染色体号
                    if (i == 0)
                    {
                        chrTmp = sampleNameTmp;
                        refStartTmp = stoi(lineVec[3*i + 1]);

                        // 初始化 map 
                        if (coorChrLociSampleLociMap.find(chrTmp) == coorChrLociSampleLociMap.end())
                        {
                            coorChrLociSampleLociMap[chrTmp];
                        }
                        if (coorChrLociSampleLociMap[chrTmp].find(refStartTmp) == coorChrLociSampleLociMap[chrTmp].end())
                        {
                            coorChrLociSampleLociMap[chrTmp][refStartTmp];
                        }
                        
                        sampleNameTmp = "reference";
                    }

                    // 记录 'sampleName'，只记录第一行的
                    if (!sampleNameVecBool)
                    {
                        sampleNameVec.push_back(sampleNameTmp);
                    }
                    
                    // 添加到总的哈希表中
                    coorChrLociSampleLociMap[chrTmp][refStartTmp][sampleNameTmp] = make_tuple(qryStartTmp, qryEndTmp);
                }

                // 只记录第一行的 'sampleName'
                sampleNameVecBool = true;
            }
        }
    };


    // 坐标转换
    class COORTRANS
    {
    private:
        // show-align 输出结果 'xx.aligns'
        string aliFileName_;

        // sampleName
        string sampleName;

        // open file
        unique_ptr<GzChunkReader> GzChunkReaderClass_;
        bool endBool = false;  // 记录文件是否遍历完

        // 记录是否要循环接下来的行
        bool chrBool = true;  // 记录该染色体对是否符合要求，目前是染色体号是否一致，不一致跳过
        bool aliBool = true;  // 记录该alignment的起始和终止是否是想要的
        bool findChrBool = true;  // 记录要找的染色体号和当前染色体号是否一致

        // 遍历文件时更改的值
        // 染色体号
        string chr = "";
        // alignment的方向、起始和终止
        string aliRefStrand = "";
        uint32_t aliRefStart = 0;
        uint32_t aliRefEnd = 0;
        string aliQryStrand = "";
        uint32_t aliQryStart = 0;
        uint32_t aliQryEnd = 0;
        // 文件流指针所在行的比对信息
        uint32_t refStart = 0;
        string refSeq = "";
        uint32_t refEnd = 0;
        uint32_t qryStart = 0;
        string qrySeq = "";
        uint32_t qryEnd = 0;
        // 记录遍历到了 '_ref_seq' 以及对应的坐标
        uint32_t refLoci = 0;
        int idxTmp = -1;
        uint32_t qryLoci = 0;

        tuple<string, uint32_t, uint32_t> get_alignment_loc(
            string infoTmp
        );
        int next_loci();  // 下一个比对段的坐标
    public:
        COORTRANS() {}

        explicit COORTRANS(string aliFileName) 
            : aliFileName_(aliFileName), GzChunkReaderClass_(std::make_unique<GzChunkReader>(aliFileName_, 1024 * 1024 * readBuffer)) {}

        COORTRANS(const COORTRANS& other) = delete;  // 禁止拷贝构造函数
        COORTRANS& operator=(const COORTRANS& other) = delete;  // 禁止拷贝赋值运算符

        // 移动构造函数
        COORTRANS(COORTRANS&& other) noexcept : 
            aliFileName_(std::move(other.aliFileName_)), 
            sampleName(std::move(other.sampleName)),
            GzChunkReaderClass_(std::move(other.GzChunkReaderClass_)),
            endBool(other.endBool),
            chrBool(other.chrBool),
            aliBool(other.aliBool),
            findChrBool(other.findChrBool),
            chr(std::move(other.chr)),
            aliRefStrand(std::move(other.aliRefStrand)),
            aliRefStart(other.aliRefStart),
            aliRefEnd(other.aliRefEnd),
            aliQryStrand(std::move(other.aliQryStrand)),
            aliQryStart(other.aliQryStart),
            aliQryEnd(other.aliQryEnd),
            refStart(other.refStart),
            refSeq(std::move(other.refSeq)),
            refEnd(other.refEnd),
            qryStart(other.qryStart),
            qrySeq(std::move(other.qrySeq)),
            qryEnd(other.qryEnd),
            refLoci(other.refLoci),
            idxTmp(other.idxTmp),
            qryLoci(other.qryLoci)
        {
            // 将原对象状态设为无效
            other.GzChunkReaderClass_ = nullptr;
        }

        // 移动赋值运算符
        COORTRANS& operator=(COORTRANS&& other) {
            if (this != &other) {
                aliFileName_ = std::move(other.aliFileName_);
                sampleName = std::move(other.sampleName);
                GzChunkReaderClass_ = std::move(other.GzChunkReaderClass_);
                endBool = other.endBool;
                chrBool = other.chrBool;
                aliBool = other.aliBool;
                findChrBool = other.findChrBool;
                chr = std::move(other.chr);
                aliRefStrand = std::move(other.aliRefStrand);
                aliRefStart = other.aliRefStart;
                aliRefEnd = other.aliRefEnd;
                aliQryStrand = std::move(other.aliQryStrand);
                aliQryStart = other.aliQryStart;
                aliQryEnd = other.aliQryEnd;
                refStart = other.refStart;
                refSeq = std::move(other.refSeq);
                refEnd = other.refEnd;
                qryStart = other.qryStart;
                qrySeq = std::move(other.qrySeq);
                qryEnd = other.qryEnd;
                refLoci = other.refLoci;
                idxTmp = other.idxTmp;
                qryLoci = other.qryLoci;
            }
            return *this;
        }
        
        // 找坐标
        uint32_t find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };


    // 获取共线性坐标
    class SYRIOUT
    {
    private:
        // syri.out
        string fileName_ = "";

        // open file
        unique_ptr<GzChunkReader> GzChunkReaderClass_;
        bool endBool = false;  // 记录文件是否遍历完

        // 记录共线性信息
        string chr = "";
        uint32_t refStart = 0;
        uint32_t refEnd = 0;
        uint32_t qryStart = 0;
        uint32_t qryEnd = 0;

        // 私有函数
        int next_loci();  // 下一个比对段的坐标
    public:
        SYRIOUT() : endBool(true) {}
        explicit SYRIOUT(string fileName) : 
            fileName_(fileName), GzChunkReaderClass_(std::make_unique<GzChunkReader>(fileName_, 1024 * 1024 * readBuffer)), endBool(false) {}

        SYRIOUT(const SYRIOUT& other) = delete;  // 禁止拷贝构造函数
        SYRIOUT& operator=(const SYRIOUT& other) = delete;  // 禁止拷贝赋值运算符

        // 移动构造函数
        SYRIOUT(SYRIOUT&& other) noexcept : 
            fileName_(std::move(other.fileName_)), 
            GzChunkReaderClass_(std::move(other.GzChunkReaderClass_)),
            endBool(other.endBool),
            chr(std::move(other.chr)),
            refStart(other.refStart),
            refEnd(other.refEnd),
            qryStart(other.qryStart),
            qryEnd(other.qryEnd) 
        {
            // 将原对象状态设为无效
            other.GzChunkReaderClass_ = nullptr;
        }

        // 移动赋值运算符
        SYRIOUT& operator=(SYRIOUT&& other) {
            if (this != &other) {
                fileName_ = std::move(other.fileName_);
                GzChunkReaderClass_ = std::move(other.GzChunkReaderClass_);
                endBool = other.endBool;
                chr = std::move(other.chr);
                refStart = other.refStart;
                refEnd = other.refEnd;
                qryStart = other.qryStart;
                qryEnd = other.qryEnd;
            }
            return *this;
        }

        int find_loci(
            string chr_, 
            uint32_t refLoci_
        );
    };


    /**
     * @brief 构建参考基因组索引
     * 
     * @param referenceFileName -> refgenome
     * 
     * @return refLenMap  map<chr, length>
    **/
    map<string, uint32_t> build_fasta_index(
        const string & referenceFileName
    );


    /* ************************************** calculate syntenic diversity ************************************** */

    /**
     * @brief 计算
     * 
     * @param sampleName               样品名
     * @param refLenMap                参考基因组长度
     * @param SynCoorTmp               共线性坐标
     * @param sampleSampleSyriMap      syri.out 输出路径字典
     * @param alignsMap                show-aligns 输出路径字典
     * 
     * @return CALSTRUCTURE
    **/
    CALSTRUCTURE calculate_run(
        const string sampleName, 
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap
    );


    /* ********************************************* memory-saving mode ********************************************* */

    /**
     * @brief 合并结果
     * 
     * @param calOutStr                 某一样品的所有计算结果
     * @param refLenMap                 染色体长度信息  map<string, length>
     * @param chrLociSynNumMap          保存最终的结果map<chr, vector<synNum> >
     * 
     * @return 0
    **/
    int merge(
        const CALSTRUCTURE& calOutStr,
        const map<string, uint32_t> & refLenMap,
        map<string, vector<uint32_t> >& chrSynNumVecMap
    );


    /**
     * @brief 计算
     * 
     * @param refLenMap                参考基因组长度
     * @param SynCoorTmp               共线性坐标
     * @param sampleSampleSyriMap      syri.out 输出路径字典
     * @param alignsMap                show-aligns 输出路径字典
     * @param allSynNum                所有组合数量
     * @param outputFileName           输出文件名
     * 
     * @return 0
    **/
    int calculate(
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap, 
        const uint32_t & allSynNum, 
        const string & outputFileName
    );



    /* ********************************************* quick mode ********************************************* */

    /**
     * @brief 合并结果
     * 
     * @param chromosome                         染色体号
     * @param chrLen                             染色体长度
     * @param CALSTRUCTUREVec                    多线程输出结果
     * 
     * @return tuple<chromosome, synOutVecTmp>   tuple<chr, vector<synNum> >
    **/
    tuple<string, vector<uint32_t> > merge_fast(
        string chromosome, 
        uint32_t chrLen, 
        const vector<CALSTRUCTURE> & CALSTRUCTUREVec
    );


    /**
     * @brief 计算
     * 
     * @param refLenMap                参考基因组长度
     * @param SynCoorTmp               共线性坐标
     * @param sampleSampleSyriMap      syri.out 输出路径字典
     * @param alignsMap                show-aligns 输出路径字典
     * @param allSynNum                所有组合数量
     * @param outputFileName           输出文件名
     * 
     * @return 0
    **/
    int calculate_fast(
        const map<string, uint32_t> & refLenMap, 
        const SYNCOOR & SynCoorTmp, 
        const map<string, map<string, string> > & sampleSampleSyriMap, 
        const map<string, string> & alignsMap, 
        const uint32_t & allSynNum, 
        const string & outputFileName
    );
}

#endif