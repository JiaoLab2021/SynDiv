#!/home/ltan/anaconda3/envs/syri_env/bin/python3
# coding=gb2312
# Created on 2023/4/12
# @author: Zezhen Du
# email: dzz0539@gmail.com or dzz0539@163.com


import os
import sys
import shutil
import datetime
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import logging
import pandas as pd
import syri_syn


# 创建目录
def makedir(
    path_dir: str
):
    """
    :param path_dir: 需要创建的文件夹路径
    :return: 0
    """
    # log
    logger = logging.getLogger('makedir')

    if os.path.isdir(path_dir):
        shutil.rmtree(path_dir)
        os.makedirs(path_dir)
        log = "[" + str(datetime.datetime.now()).split(".")[0] + \
                '] [makedir] \'{}\' already exists, clear and recreate.'.format(path_dir)
        logger.error(log)
    else:
        os.makedirs(path_dir)


class AlignmentLociClass:
    """
    需要比对的非共线性坐标
    """
    def __init__(
            self,
            ali_loci_filename: str, 
            configFileMap: dict, 
            args, 
            code_path
    ):
        """
        构造函数
        :param ali_loci_filename: 需要比对的坐标  no_syn输出
        :param configFileMap: getParser输出
        :param args: getParser输出
        """

        # 需要比对的坐标文件
        self._ali_loci_filename = ali_loci_filename

        # 基因组文件字典
        self._configFileMap = configFileMap  # dict{sample, {genome: "path", aligns: "path", syri_out: "path"}}

        # 输入参数
        self._args = args

        # 环境变量
        self._path = code_path

        # 索引结果
        self._ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >

        # 构建索引
        self._index()

        # syntenic result
        self._sample1Sample2SynPD = {}  # map<sample1, map<sample2, pd.DataFrame()> >


    def _index(
        self
    ):
        """
        构建坐标索引
        """
        with open(self._ali_loci_filename) as f:
            for infos in f.readlines():
                # 遇到空行跳过
                if len(infos) == 0:
                    continue

                # 去字符串并拆分
                infos_list = infos.strip().split()
                
                # 如果只有一个坐标，跳过该行
                if len(infos_list) <= 3:
                    continue

                # 染色体号
                chr_tmp = infos_list[0]

                # 临时键
                key_tmp = infos_list[0] + "_" + infos_list[1]

                # 初始化字典
                self._ali_loci_dict[key_tmp] = {}

                # 从第三列开始循环
                for idx1 in range(2, len(infos_list)):
                    info = infos_list[idx1]

                    info_list = info.strip().split(":")
                    sample_tmp = info_list[0]

                    loci_list = info_list[1].split(";")
                    
                    # 初始化字典
                    self._ali_loci_dict[key_tmp][sample_tmp] = []

                    for idx2 in range(len(loci_list)):
                        # 如果只有一个坐标，跳过
                        if len(loci_list[idx2].split("-")) != 2:
                            continue

                        # 记录坐标   map<sample, vector<chr:start-end>
                        self._ali_loci_dict[key_tmp][sample_tmp].append(chr_tmp + ":" + loci_list[idx2])

    
    def _alignment(
        self, 
        key1, 
        value1, 
        threads, 
        debug
    ):
        """
        提取序列并比对的多线程函数
        : param key1  参考基因组坐标, map<chr+"_"+start, map<sample, vector<chr:start-end> > >
        : param value1  map<sample, vector<chr:start-end> >
        : param threads  线程数
        : param debug    是否调试代码
        : return sample1Sample2SynPDTmp, map<sample1, map<sample2, pd.DataFrame()> >
        """
        logger = logging.getLogger('_alignment')

        # 临时存储共线性坐标
        sample1Sample2SynPDTmp = {}  # map<sample1, map<sample2, pd.DataFrame()> >

        # 存储提取的序列文件
        fileNameMap = {key2: os.path.abspath(f"{key1}_{key2}.fa") 
                    for key2 in value1.keys()}
        fileNameLocMap = {key2: {} for key2 in value1.keys()}
        sampleList = sorted(value1.keys())

        # 提取序列
        for key2, value2 in value1.items():  # map<sample, vector<chr:start-end> >
            # 提取的序列保存的地方
            outputFileName = fileNameMap[key2]

            # ########## 提取序列 ######### #
            for idx1, loci in enumerate(value2):
                # 检查染色体号是否符合规定  (chr1:257130-257225)
                if ":" not in loci or "-" not in loci:
                    raise ValueError(f"Error: incorrect '{loci}' found in '{self._ali_loci_filename}'.")
                    
                cmd = ""

                # 如果为第一个序列，保留染色体号
                if idx1 == 0:
                    # SynDiv 提交的配置文件有 ["genome"] 键
                    try:
                        cmd = f'echo \">{key1.split("_")[0]}\" > {outputFileName}; samtools faidx {self._configFileMap[key2]["genome"]} {loci} | grep -v \">\" >> {outputFileName}'
                    # SynDiv_ 提交的配置文件有没有 ["genome"] 键
                    except TypeError:
                        cmd = f'echo \">{key1.split("_")[0]}\" > {outputFileName}; samtools faidx {self._configFileMap[key2]} {loci} | grep -v \">\" >> {outputFileName}'
                else:  # 如果为之后的序列，不要染色体号
                    # SynDiv 提交的配置文件有 ["genome"] 键
                    try:
                        cmd = f'samtools faidx {self._configFileMap[key2]["genome"]} {loci} | grep -v \">\" >> {outputFileName}'
                    # SynDiv_ 提交的配置文件有没有 ["genome"] 键
                    except TypeError:
                        cmd = f'echo \">{key1.split("_")[0]}\" > {outputFileName}; samtools faidx {self._configFileMap[key2]} {loci} | grep -v \">\" >> {outputFileName}'

                # 提交任务
                cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=self._path)
                # 标准输出和错误输出
                stdout, stderr = cmd_out.communicate()
                exit_state = cmd_out.returncode
                if exit_state != 0 or "FileNotFoundError" in str(stderr) or \
                        "command not found" in str(stderr) or \
                        "error" in str(stderr) or \
                        "Error" in str(stderr):
                    logger.error(f"Error: Exit status is non-zero.")
                    logger.error(f"Error: CMD -> {cmd}")
                    logger.error(f"Error: stderr -> {stderr}")
                    sys.exit(1)

                # 记录坐标
                lociList = loci.split(":")[1].split("-")
                fileNameLocMap[key2][int(lociList[0])] = int(lociList[1])


        # ############### alignment ############### #
        # 样品名排序
        sampleList.sort()

        # 样品1
        for idx1, sample1 in enumerate(sampleList):
            filename1 = fileNameMap[sample1]

            # ref坐标和长度信息
            refStartTmpList = []
            refEndTmpList = []
            refLenTmpList = []
            for key3, value3 in fileNameLocMap[sample1].items():  # map<start, end>
                refLenTmp = abs(value3 - key3 + 1)
                refStartTmpList.append(key3)
                refEndTmpList.append(value3)
                refLenTmpList.append(refLenTmp)

            if debug:
                logger.error('\n', refStartTmpList)
                logger.error(refEndTmpList)
                logger.error(refLenTmpList)
            
            # 样品2
            for idx2 in range(idx1 + 1, len(sampleList)):
                sample2 = sampleList[idx2]
                filename2 = fileNameMap[sample2]
                
                # 输出文件名
                alignMentFileName = os.path.abspath(f'{key1}_{sample1}_{sample2}.sam')

                # minimap2 比对
                cmd = f'minimap2 -t {threads} -ax asm5 --eqx {filename1} {filename2} -o {alignMentFileName}'
                # 提交任务
                cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=self._path)
                # 标准输出和错误输出
                stdout, stderr = cmd_out.communicate()
                exit_state = cmd_out.returncode
                if exit_state != 0 or "FileNotFoundError" in str(stderr) or \
                        "command not found" in str(stderr) or \
                        "error" in str(stderr) or \
                        "Error" in str(stderr):
                    logger.error(f"Error: Exit status is non-zero.")
                    logger.error(f"Error: CMD -> {cmd}")
                    logger.error(f"Error: stderr -> {stderr}")
                    sys.exit(1)

                # syri 寻找共线性坐标
                # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
                synDatas = syri_syn.main(self._args, alignMentFileName)

                # remove file
                if os.path.exists(alignMentFileName) and debug==False:
                    os.remove(alignMentFileName)

                # 如果synDatas为空，下一个循环
                if synDatas.empty:
                    continue

                # 调试代码
                if debug:
                    logger.error(f"alignMentFileName:{alignMentFileName}\nsynDatas:{synDatas}")

                # 保存坐标到总表中
                # 按行打印每一行的元素
                # 用于判断坐标怎么计算
                idxTmp = 0  # 临时索引
                preLenTmp = 0  # 前边序列的长度(n)，同于判断是不是一个独立的比对段
                thisLenTmp = refLenTmpList[0]  # 前边序列的长度加上往下一个索引(n+1)，最大值代表下边没有了
                for index, row in synDatas.iterrows():
                    try:
                        # 临时的 DataFrame
                        synDataTmp1 = synDatas.iloc[[index]].copy()

                        if debug:
                            logger.error(f"synDataTmp1:{synDataTmp1}")

                        # 如果为空集或者反向比对，下一个循环
                        if synDataTmp1.empty or synDataTmp1.iloc[0, 7] == "-" or synDataTmp1.iloc[0, 8] == "-":
                            continue

                        # 判断坐标位置
                        while True and idxTmp < len(refLenTmpList):
                            # 如果起始大于当前所有序列的长度，更新索引和长度
                            while synDataTmp1.iloc[0, 0] > thisLenTmp and idxTmp < len(refLenTmpList):
                                idxTmp += 1
                                preLenTmp = sum(refLenTmpList[:idxTmp])

                                # 如果下边还有坐标，更新下一个长度
                                if idxTmp + 1 < len(refLenTmpList):
                                    thisLenTmp = sum(refLenTmpList[:idxTmp + 1])
                                else:  # 如果没有，赋值为最大值
                                    thisLenTmp = sum(refLenTmpList)

                            # 是独立的比对
                            if synDataTmp1.iloc[0, 0] >= preLenTmp and synDataTmp1.iloc[0, 1] <= thisLenTmp:
                                synDataTmp1.iloc[0, 0] += refStartTmpList[idxTmp] - 1 - preLenTmp
                                synDataTmp1.iloc[0, 1] += refStartTmpList[idxTmp] - 1 - preLenTmp
                                try:
                                    # 赋值
                                    sample1Sample2SynPDTmp[sample1][sample2] = pd.concat([sample1Sample2SynPDTmp[sample1][sample2], synDataTmp1], ignore_index=True)
                                except KeyError:
                                    # 初始化
                                    if sample1 not in sample1Sample2SynPDTmp:
                                        sample1Sample2SynPDTmp[sample1] = {}
                                    if sample2 not in sample1Sample2SynPDTmp[sample1]:
                                        sample1Sample2SynPDTmp[sample1][sample2] = pd.DataFrame()
                                    # 赋值
                                    sample1Sample2SynPDTmp[sample1][sample2] = pd.concat([sample1Sample2SynPDTmp[sample1][sample2], synDataTmp1], ignore_index=True)

                                if debug:
                                    logger.error(f"\n{synDataTmp1}\n")

                                break
                            # 跨越了两个比对序列
                            else:
                                # 将当前序列的比对坐标提取并保存
                                synDataTmp2 = synDataTmp1.copy()
                                synDataTmp2.iloc[0, 0] += refStartTmpList[idxTmp] - 1 - preLenTmp
                                synDataTmp2.iloc[0, 1] = refEndTmpList[idxTmp]
                                try:
                                    # 赋值
                                    sample1Sample2SynPDTmp[sample1][sample2] = pd.concat([sample1Sample2SynPDTmp[sample1][sample2], synDataTmp2], ignore_index=True)
                                except KeyError:
                                    # 初始化
                                    if sample1 not in sample1Sample2SynPDTmp:
                                        sample1Sample2SynPDTmp[sample1] = {}
                                    if sample2 not in sample1Sample2SynPDTmp[sample1]:
                                        sample1Sample2SynPDTmp[sample1][sample2] = pd.DataFrame()
                                    # 赋值
                                    sample1Sample2SynPDTmp[sample1][sample2] = pd.concat([sample1Sample2SynPDTmp[sample1][sample2], synDataTmp2], ignore_index=True)
                                if debug:
                                    logger.error(f"\n{synDataTmp2}\n")
                                # 更新synDataTmp1的坐标，计算剩下的坐标
                                synDataTmp1.iloc[0, 0] = thisLenTmp + 1
                    except KeyError as e:  # 处理异常
                        logger.error(f"KeyError occurred:{e}")
                        continue
                    except IndexError as e:  # 处理异常
                        logger.error(f"IndexError occurred:{e}")
                        if debug:
                            logger.error(f'Error refStartTmpList:{refStartTmpList}')
                            logger.error(f"Error refEndTmpList:{refEndTmpList}")
                            logger.error(f"Error refLenTmpList:{refLenTmpList}")
                            logger.error(f"Error synDataTmp1:{synDataTmp1}")
                            logger.error(f"Error index:{index}", index)
                            logger.error(f"Error alignMentFileName:{alignMentFileName}")
                        continue

        # remove file
        for key2, value2 in fileNameMap.items():
            if os.path.exists(value2) and debug==False:
                os.remove(value2)

        # 返回当前位置的所有样品组合的共线性坐标
        return sample1Sample2SynPDTmp


def main(args, configFileMap, no_synPath, workDir, code_path):
    """
    :param args:  getParser输出参数
    :param configFileMap:   getParser输出每个样品所有需要样品的路径
    :param no_synPath:  Nexus_c no_syn 输出的非共线性在qurey上的坐标
    :param workDir:  工作路径
    :param code_path: 环境变量
    :return: coorPath
    """
    # 创建工作路径
    makedir(workDir)
    workDir = workDir + os.sep
    # 切换工作路径
    os.chdir(workDir)

    logger = logging.getLogger('no_syn_alignment')

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Running.')

    # 构造类
    ali_loci_class = AlignmentLociClass(no_synPath, configFileMap, args, code_path)
    # 构建索引
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Building index.')
    ali_loci_class._index()

    # 多线程结果
    output = []

    # 线程池提交任务
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Alignment.')
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        # 存储所有的 future 对象
        results = []
        # 遍历 ali_loci_class._ali_loci_dict 字典，并多线程执行 ali_loci_class._alignmentRun 函数
        for key1, value1 in ali_loci_class._ali_loci_dict.items():  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >
            future = executor.submit(ali_loci_class._alignment, key1, value1, args.threads, args.debug)
            results.append(future)
            # 控制多线程执行间隔
            time.sleep(0.00005)

        # 等待所有 future 对象执行完毕并获取结果
        for future in as_completed(results):
            result = future.result()
            output.append(result)

    # 合并结果
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Merge result.')
    for sample1Sample2SynPDTmp in output:
        for key2, value2 in sample1Sample2SynPDTmp.items():  # map<sample1, map<sample2, pd.DataFrame()> >
            for key3, value3 in value2.items():  # map<sample2, pd.DataFrame()>
                try:
                    # 赋值
                    ali_loci_class._sample1Sample2SynPD[key2][key3] = pd.concat([ali_loci_class._sample1Sample2SynPD[key2][key3], value3], ignore_index=True)
                except KeyError:
                    # 初始化
                    if key2 not in ali_loci_class._sample1Sample2SynPD:
                        ali_loci_class._sample1Sample2SynPD[key2] = {}
                    if key3 not in ali_loci_class._sample1Sample2SynPD[key2]:
                        ali_loci_class._sample1Sample2SynPD[key2][key3] = pd.DataFrame()
                    # 赋值
                    ali_loci_class._sample1Sample2SynPD[key2][key3] = pd.concat([ali_loci_class._sample1Sample2SynPD[key2][key3], value3], ignore_index=True)


    # 对 _sample1Sample2SynPD 进行排序并输出
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Sort and save the results.')
    # All output file information
    allSyriOutInfo = ""
    for key1, value1 in ali_loci_class._sample1Sample2SynPD.items():  # map<sample1, map<sample2, pd.DataFrame()> >
        for key2, value2 in value1.items():  # map<sample2, pd.DataFrame()>
            # pd.DataFrame() 排序
            value2.sort_values(by=['aChr', 'aStart', 'aEnd'], ascending=[True, True, True], inplace=True)

            # 保存输出字符串
            outTxt = ""

            # 将结果写入输出字符串中
            for row in value2.itertuples():  # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
                outTxt += '\t'.join(
                    [row.aChr, 
                    str(row.aStart), 
                    str(row.aEnd), 
                    "-", 
                    "-", 
                    row.aChr, 
                    str(row.bStart), 
                    str(row.bEnd), 
                    "SYNAL", 
                    "SYN", 
                    "SYNAL", 
                    "-"]) + "\n"

            # 输出到文件
            outputFilePath = os.path.abspath("{}_{}.syn.out".format(key1, key2))
            with open(outputFilePath, 'w', encoding="utf-8") as outputFile:
                outputFile.write(outTxt)

            # 输出文件的信息
            allSyriOutInfo += '\t'.join([key1, key2, outputFilePath]) + "\n"

    # 输出文件信息
    no_syn_alignmentPath = os.path.join(os.getcwd() + os.sep, f"{args.prefix}no_syn_alignment.out")
    with open(no_syn_alignmentPath, 'w', encoding="utf-8") as f:
        f.write(allSyriOutInfo)

    # 关闭线程池
    executor.shutdown(wait=True)

    return no_syn_alignmentPath
