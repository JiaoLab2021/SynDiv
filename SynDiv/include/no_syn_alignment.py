#!/home/ltan/anaconda3/envs/syri_env/bin/python3
# coding=gb2312
# Created on 2023/4/12
# @author: Zezhen Du
# email: dzz0539@gmail.com or dzz0539@163.com


import os
import shutil
import subprocess
import multiprocessing
import logging
import pandas as pd
import syri_syn


# 定义日志格式
logger = logging.getLogger('no_syn_alignment')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # 输出到控制台
handler.setFormatter(formatter)
logger.addHandler(handler)


# 创建目录
def makedir(
    path_dir: str
):
    """
    :param path_dir: 需要创建的文件夹路径
    :return: 0
    """
    if os.path.isdir(path_dir):
        shutil.rmtree(path_dir)
        os.makedirs(path_dir)
        log = '[makedir] \'{}\' already exists, clear and recreate.'.format(path_dir)
        logger.error(log)
    else:
        os.makedirs(path_dir)


# 运行多线程任务
def run_command(command, envPath):
    # 提交任务
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=envPath)

    # 获取退出状态
    exit_state = proc.wait()

    stdout, stderr = proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()

    # 标准输出和错误输出
    if exit_state != 0 or "FileNotFoundError" in str(stderr) or \
            "command not found" in str(stderr) or \
            "error" in str(stderr) or \
            "Error" in str(stderr):
        logger.error(f"Error: Exit status is {exit_state}")
        logger.error(f"Error: CMD -> {command}")
        logger.error(f"Error: stderr -> {stderr}")

    return exit_state, stdout, stderr


class AlignmentLociClass:
    """
    需要比对的非共线性坐标
    """
    def __init__(
            self,
            ali_loci_filename: str, 
            config_file_map: dict, 
            args, 
            code_path
    ):
        """
        构造函数
        :param ali_loci_filename: 需要比对的坐标  no_syn输出
        :param config_file_map: getParser输出
        :param args: getParser输出
        """
        # 需要比对的坐标文件
        self._ali_loci_filename = ali_loci_filename

        # 基因组文件字典
        self._config_file_map = config_file_map  # dict{sample, {genome: "path", aligns: "path", syri_out: "path"}}

        # 输入参数
        self._args = args

        # 环境变量
        self._path = code_path

        # 需要比对的坐标
        self._long_ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >， 列数过多的行
        self._short_ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >，列数较少的行

        # syntenic result
        self._sample1_sample2_synpd = {}  # map<sample1, map<sample2, pd.DataFrame()> >


    def _index(
        self
    ):
        """
        构建坐标索引
        """
        # 所有需要比对的坐标
        ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >

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
                ali_loci_dict[key_tmp] = {}

                # 从第三列开始循环
                for idx1 in range(2, len(infos_list)):
                    info = infos_list[idx1]

                    info_list = info.strip().split(":")
                    sample_tmp = info_list[0]

                    loci_list = info_list[1].split(";")
                    
                    # 初始化字典
                    ali_loci_dict[key_tmp][sample_tmp] = []

                    for idx2 in range(len(loci_list)):
                        # 如果只有一个坐标，跳过
                        if len(loci_list[idx2].split("-")) != 2:
                            continue

                        # 记录坐标   map<sample, vector<chr:start-end>
                        ali_loci_dict[key_tmp][sample_tmp].append(chr_tmp + ":" + loci_list[idx2])

        # 使用字典推导式和条件表达式将 ali_loci_dict 拆分为列数长和短的两个字典
        import math
        threshold_number = (-1 + math.sqrt(1 + 8 * max(1, self._args.jobs * 100))) / 2  # y = (x * (x + 1)) / 2    x = (-1 + sqrt(1 + 8 * y)) / 2
        self._long_ali_loci_dict = {k: v for k, v in ali_loci_dict.items() if len(v) > threshold_number}
        self._short_ali_loci_dict = {k: v for k, v in ali_loci_dict.items() if len(v) <= threshold_number}
        # 打印数量
        logger.error(f'Number of tasks running in parallel by column: {len(self._long_ali_loci_dict)}')
        logger.error(f'Number of tasks running in parallel by row: {len(self._short_ali_loci_dict)}')



    def _alignment_line(
        self, 
        key_value
    ):
        """
        提取序列并比对的多线程函数
        : param key_value   参考基因组坐标_map<sample, vector<chr:start-end> >   map<sample, vector<chr:start-end> >
        : param key1        参考基因组坐标, chr+"_"+start
        : param value1      map<sample, vector<chr:start-end> >
        :
        : return key1, sample1_sample2_synpd, map<sample1, map<sample2, pd.DataFrame()> >
        """
        key1, value1 = key_value
        # 临时存储共线性坐标
        sample1_sample2_synpd = {}  # map<sample1, map<sample2, pd.DataFrame()> >

        # 存储提取的序列文件
        file_name_map = {key2: os.path.abspath(f"{key1}_{key2}.fa") 
                    for key2 in value1.keys()}
        file_name_loc_map = {key2: {} for key2 in value1.keys()}
        sample_list = sorted(value1.keys())

        # 染色体号
        chromosome = ""

        # ############### subseq ############### #
        for key2, value2 in value1.items():  # map<sample, vector<chr:start-end> >
            # 提取的序列保存的地方
            output_file_name = file_name_map[key2]

            # 提交任务
            chromosome, file_name_loc_map_tmp = samtools(
                self._path, 
                self._ali_loci_filename,
                key1, 
                key2, 
                value2, 
                self._config_file_map, 
                output_file_name
            )

            # 记录该样品的坐标信息
            for sample2_tmp, loci_map in file_name_loc_map_tmp.items():
                file_name_loc_map[sample2_tmp] = loci_map


        # ############### minimap2 + syri ############### #
        # 记录和sample1共线性的sample，用于判断之后的sample之间是否需要继续相互比对
        syn_sample_dict = {}  # map<sample1, map<sample2, 0/1> >   0->noSyn  1->Syn

        # 样品名排序
        sample_list.sort()

        # 样品1
        for idx1, sample1 in enumerate(sample_list):
            filename1 = file_name_map[sample1]

            # 记录和sample1共线性的sample，用于判断之后的sample之间是否需要继续相互比对
            synSampleList = []
            # 记录和sample1非共线性的sample，用于判断之后的sample之间是否需要继续相互比对
            noSynSampleList = []

            # ref坐标和长度信息
            ref_start_tmp_list = []
            ref_end_tmp_list = []
            ref_len_tmp_list = []
            for key3, value3 in file_name_loc_map[sample1].items():  # map<start, end>
                refLenTmp = abs(value3 - key3 + 1)
                ref_start_tmp_list.append(key3)
                ref_end_tmp_list.append(value3)
                ref_len_tmp_list.append(refLenTmp)

            if self._args.debug:
                logger.error('\n', ref_start_tmp_list)
                logger.error(ref_end_tmp_list)
                logger.error(ref_len_tmp_list)
            
            # 样品2
            for idx2 in range(idx1 + 1, len(sample_list)):
                sample2 = sample_list[idx2]
                filename2 = file_name_map[sample2]

                # qry坐标和长度信息
                # qryStartTmpList = []
                # qryEndTmpList = []
                # qryLenTmpList = []
                qryLenSum = 0
                for key3, value3 in file_name_loc_map[sample2].items():  # map<start, end>
                    qryLenTmp = abs(value3 - key3 + 1)
                    # qryStartTmpList.append(key3)
                    # qryEndTmpList.append(value3)
                    # qryLenTmpList.append(qryLenTmp)
                    qryLenSum += qryLenTmp

                # 判断是否需要比对和计算共线性
                if sample1 in syn_sample_dict:
                    if sample2 in syn_sample_dict[sample1]:
                        # 1代表 self._args.synRatio 为共线性，直接添加
                        if syn_sample_dict[sample1][sample2] == 1:
                            listSize = len(ref_start_tmp_list)
                            # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
                            SynPDTmp = pd.DataFrame({'aStart': ref_start_tmp_list, 
                                                    'aEnd': ref_end_tmp_list, 
                                                    'bStart': [1]*listSize, 
                                                    'bEnd': [1]*listSize, 
                                                    'aLen': ref_len_tmp_list, 
                                                    'bLen': [1]*listSize, 
                                                    'iden': [100]*listSize, 
                                                    'aDir': [1]*listSize, 
                                                    'bDir': [1]*listSize, 
                                                    'aChr': [chromosome]*listSize, 
                                                    'bChr': [chromosome]*listSize})
                            try:
                                # 赋值
                                sample1_sample2_synpd[sample1][sample2] = SynPDTmp
                            except KeyError:
                                # 初始化
                                if sample1 not in sample1_sample2_synpd:
                                    sample1_sample2_synpd[sample1] = {}
                                if sample2 not in sample1_sample2_synpd[sample1]:
                                    sample1_sample2_synpd[sample1][sample2] = pd.DataFrame()
                                # 赋值
                                sample1_sample2_synpd[sample1][sample2] = SynPDTmp
                            continue
                        # 0代表为非共线性，同样不比对，直接跳过该sample2
                        else:
                            continue
                
                # 输出文件名
                alignment_file_name = os.path.abspath(f'{key1}_{sample1}_{sample2}.sam')

                # 提交任务
                # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
                synpd = minimap2_syn(self._args, self._path, filename1, filename2, alignment_file_name, ref_start_tmp_list, ref_end_tmp_list, ref_len_tmp_list)

                # 如果为空值，直接下一个循环
                if synpd.empty:
                    # 记录为非共线性的样本
                    noSynSampleList.append(sample2)
                    continue

                try:
                    # 赋值
                    sample1_sample2_synpd[sample1][sample2] = synpd
                except KeyError:
                    # 初始化
                    if sample1 not in sample1_sample2_synpd:
                        sample1_sample2_synpd[sample1] = {}
                    if sample2 not in sample1_sample2_synpd[sample1]:
                        sample1_sample2_synpd[sample1][sample2] = pd.DataFrame()
                    # 赋值
                    sample1_sample2_synpd[sample1][sample2] = synpd

                # 在这里添加一个列表，把是共线性的样品全部收集起来，后续比对跳过这些样品  要大于阈值 共线性长度占 self._args.synRatio 以上
                qryTotal = synpd['bLen'].sum()
                alignRatio = float(qryTotal)/max(qryLenSum, 1)
                if alignRatio >= self._args.synRatio:
                    synSampleList.append(sample2)
                elif alignRatio <= self._args.nosynRatio:
                    # 记录为非共线性的样本
                    noSynSampleList.append(sample2)

            # 将共线性的样品添加到总的dict中，用于之后比对查询
            # 样品1
            for idx2, sample2 in enumerate(synSampleList):
                # 初始化 dict
                if sample2 not in syn_sample_dict:
                    syn_sample_dict[sample2] = {}
                # 样品2
                for idx3 in range(idx2 + 1, len(synSampleList)):
                    syn_sample_dict[sample2][synSampleList[idx3]] = 1
            # 将非共线性的样品添加到总的dict中，用于之后比对查询
            # 样品1
            for sample2 in synSampleList:
                # 初始化 dict
                if sample2 not in syn_sample_dict:
                    syn_sample_dict[sample2] = {}
                # 样品2
                for sample3 in noSynSampleList:
                    syn_sample_dict[sample2][sample3] = 0

        # remove file (fasta)
        for key2, value2 in file_name_map.items():
            if os.path.exists(value2) and self._args.debug==False:
                os.remove(value2)

        # 返回当前位置的所有样品组合的共线性坐标
        return key1, sample1_sample2_synpd
    

    def _alignment_column(
        self, 
        key1, 
        value1, 
        pool
    ):
        """
        提取序列并比对的多线程函数
        : param key1    参考基因组坐标, map<chr+"_"+start, map<sample, vector<chr:start-end> > >
        : param value1  map<sample, vector<chr:start-end> >
        : param pool    进程池
        :
        : return key1, sample1_sample2_synpd_tmp, map<sample1, map<sample2, pd.DataFrame()> >
        """
        # 临时存储共线性坐标
        sample1_sample2_synpd_tmp = {}  # map<sample1, map<sample2, pd.DataFrame()> >

        # 存储提取的序列文件
        file_name_map = {key2: os.path.abspath(f"{key1}_{key2}.fa") 
                    for key2 in value1.keys()}
        file_name_loc_map = {key2: {} for key2 in value1.keys()}
        sample_list = sorted(value1.keys())

        # ############### subseq ############### #
        results = []
        for key2, value2 in value1.items():  # map<sample, vector<chr:start-end> >
            # 提取的序列保存的地方
            output_file_name = file_name_map[key2]

            # 提交任务
            result = pool.apply_async(
                samtools, 
                args=(
                    self._path, 
                    self._ali_loci_filename,
                    key1, 
                    key2, 
                    value2, 
                    self._config_file_map, 
                    output_file_name,
                )
            )
            results.append(result)

        # 获取函数的返回值
        for result in results:
            chromosome, file_name_loc_map_tmp = result.get()

            # 记录该样品的坐标信息
            for sample2_tmp, loci_map in file_name_loc_map_tmp.items():
                file_name_loc_map[sample2_tmp] = loci_map


        # ############### minimap2 + syri ############### #
        # 样品名排序
        sample_list.sort()

        # 存储比对的最终结果
        results = []
        
        # 样品1
        for idx1, sample1 in enumerate(sample_list):
            filename1 = file_name_map[sample1]
            
            # ref坐标和长度信息
            ref_start_tmp_list = []
            ref_end_tmp_list = []
            ref_len_tmp_list = []
            for key3, value3 in file_name_loc_map[sample1].items():  # map<start, end>
                refLenTmp = abs(value3 - key3 + 1)
                ref_start_tmp_list.append(key3)
                ref_end_tmp_list.append(value3)
                ref_len_tmp_list.append(refLenTmp)

            if self._args.debug:
                logger.error('\n', ref_start_tmp_list)
                logger.error(ref_end_tmp_list)
                logger.error(ref_len_tmp_list)
            
            # 样品2
            for idx2 in range(idx1 + 1, len(sample_list)):
                sample2 = sample_list[idx2]
                filename2 = file_name_map[sample2]
                
                # 输出文件名
                alignment_file_name = os.path.abspath(f'{key1}_{sample1}_{sample2}.sam')

                # minimap2 + syn
                result = pool.apply_async(
                    minimap2_syn, 
                    args=(
                        self._args, 
                        self._path, 
                        filename1, 
                        filename2, 
                        alignment_file_name, 
                        ref_start_tmp_list, 
                        ref_end_tmp_list, 
                        ref_len_tmp_list,
                    )
                )
                results.append((sample1, sample2, result))

        # 获取函数的返回值
        for result in results:
            # 获取输出结果
            sample1, sample2, synpd = result
            synpd = synpd.get()
            
            # 如果为空值，直接下一个循环
            if synpd.empty:
                continue

            try:
                # 赋值
                sample1_sample2_synpd_tmp[sample1][sample2] = synpd
            except KeyError:
                # 初始化
                if sample1 not in sample1_sample2_synpd_tmp:
                    sample1_sample2_synpd_tmp[sample1] = {}
                if sample2 not in sample1_sample2_synpd_tmp[sample1]:
                    sample1_sample2_synpd_tmp[sample1][sample2] = pd.DataFrame()
                # 赋值
                sample1_sample2_synpd_tmp[sample1][sample2] = synpd

        # remove file
        for key2, value2 in file_name_map.items():
            if os.path.exists(value2) and self._args.debug==False:
                os.remove(value2)

        # 返回当前位置的所有样品组合的共线性坐标
        return key1, sample1_sample2_synpd_tmp
    

def samtools(
    path, 
    ali_loci_filename,
    key1, 
    key2, 
    value2, 
    config_file_map, 
    output_file_name
):
    """
    :description:                 提取序列
    :param path:                  环境变量
    :param ali_loci_filename:     SynDic_c no_syn 输出的需要比对的坐标
    :param key1:                  chr+"_"+start
    :param key2:                  sample
    :param value2:                vector<chr:start-end>
    :param config_file_map:         输入的配置文件dict
    :param output_file_name:        输出文件名
    :return file_name_loc_map_tmp:    map<sample2, map<start, end> >
    """
    # ########## 提取序列 ######### #
    # 记录sample2的坐标信息
    file_name_loc_map_tmp = {}
    file_name_loc_map_tmp[key2] = {}  # map<sample2, map<start, end> >

    for idx1, loci in enumerate(value2):
        # 检查染色体号是否符合规定  (chr1:257130-257225)
        if ":" not in loci or "-" not in loci:
            raise ValueError(f"Error: incorrect '{loci}' found in '{ali_loci_filename}'.")
        
        # 提取染色体号
        chromosome = loci.split(":")[0]

        cmd = ""

        # 如果为第一个序列，保留染色体号
        if idx1 == 0:
            # SynDiv 提交的配置文件有 ["genome"] 键
            try:
                cmd = f'echo \">{key1.split("_")[0]}\" > {output_file_name}; samtools faidx {config_file_map[key2]["genome"]} {loci} | grep -v \">\" >> {output_file_name}'
            # SynDiv_ 提交的配置文件有没有 ["genome"] 键
            except TypeError:
                cmd = f'echo \">{key1.split("_")[0]}\" > {output_file_name}; samtools faidx {config_file_map[key2]} {loci} | grep -v \">\" >> {output_file_name}'
        else:  # 如果为之后的序列，不要染色体号
            # SynDiv 提交的配置文件有 ["genome"] 键
            try:
                cmd = f'samtools faidx {config_file_map[key2]["genome"]} {loci} | grep -v \">\" >> {output_file_name}'
            # SynDiv_ 提交的配置文件有没有 ["genome"] 键
            except TypeError:
                cmd = f'samtools faidx {config_file_map[key2]} {loci} | grep -v \">\" >> {output_file_name}'

        # 记录坐标
        lociList = loci.split(":")[1].split("-")
        file_name_loc_map_tmp[key2][int(lociList[0])] = int(lociList[1])

        # 提交任务
        run_command(cmd, path)

    return chromosome, file_name_loc_map_tmp


def minimap2_syn(
    args, 
    path, 
    filename1, 
    filename2, 
    alignment_file_name, 
    ref_start_tmp_list, 
    ref_end_tmp_list, 
    ref_len_tmp_list
):
    """
    :description:                 比对并寻找共线性坐标
    :param args:                  输入的参数信息
    :param path:                  环境变量
    :param filename1:             样品1名字
    :param filename2:             样品2名字
    :param alignment_file_name:     比对文件输出路径
    :param ref_start_tmp_list:       样品1的起始坐标列表
    :param ref_end_tmp_list:         样品1的终止坐标列表
    :param ref_len_tmp_list:         样品1的长度列表
    """
    # minimap2 比对
    cmd = f'minimap2 -ax asm5 --eqx {filename1} {filename2} -o {alignment_file_name}'
    # 提交任务
    run_command(cmd, path)

    # 如果文件不存在，直接返回
    if not os.path.exists(alignment_file_name):
        return pd.DataFrame()

    # Syn
    synpd = syri(args, ref_start_tmp_list, ref_end_tmp_list, ref_len_tmp_list, alignment_file_name)

    return synpd


def syri(
        args, 
        ref_start_tmp_list, 
        ref_end_tmp_list, 
        ref_len_tmp_list, 
        alignment_file_name
    ):
        """
        :description:                           解析minimap2结果, 寻找共线性坐标
        :param args:                            输入的参数信息
        :param ref_start_tmp_list:                 ref坐标和长度信息
        :param ref_end_tmp_list:                   ref坐标和长度信息
        :param ref_len_tmp_list:                   ref坐标和长度信息
        :param alignment_file_name:               minimap2 比对的结果文件
        :return: synpd                          synpd -> pd.DataFrame()
        """
        # 临时存储共线性坐标
        # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
        synpd = pd.DataFrame()

        # syri 寻找共线性信息
        synDatas = pd.DataFrame()
        try:
            synDatas = syri_syn.main(args, alignment_file_name)
        except:
            return synpd

        # remove file
        if os.path.exists(alignment_file_name) and args.debug==False:
            os.remove(alignment_file_name)

        # 如果synDatas为空，直接返回空值
        if synDatas.empty:
            return synpd
        
        # 调试代码
        if args.debug:
            logger.error(f"alignment_file_name:{alignment_file_name}\nsynDatas:{synDatas}")

        # 保存坐标到总表中
        # 按行打印每一行的元素
        # 用于判断坐标怎么计算
        idxTmp = 0  # 临时索引
        preLenTmp = 0  # 前边序列的长度(n)，同于判断是不是一个独立的比对段
        thisLenTmp = ref_len_tmp_list[0]  # 前边序列的长度加上往下一个索引(n+1)，最大值代表下边没有了
        for index, row in synDatas.iterrows():
            try:
                # 临时的 DataFrame
                synDataTmp1 = synDatas.iloc[[index]].copy()

                if args.debug:
                    logger.error(f"synDataTmp1:{synDataTmp1}")

                # 如果为空集或者反向比对，下一个循环
                if synDataTmp1.empty or synDataTmp1.iloc[0, 7] == -1 or synDataTmp1.iloc[0, 8] == -1:
                    continue

                # 判断坐标位置
                while True and idxTmp < len(ref_len_tmp_list):
                    # 如果起始大于当前所有序列的长度，更新索引和长度
                    while synDataTmp1.iloc[0, 0] > thisLenTmp and idxTmp < len(ref_len_tmp_list):
                        idxTmp += 1
                        preLenTmp = sum(ref_len_tmp_list[:idxTmp])

                        # 如果下边还有坐标，更新下一个长度
                        if idxTmp + 1 < len(ref_len_tmp_list):
                            thisLenTmp = sum(ref_len_tmp_list[:idxTmp + 1])
                        else:  # 如果没有，赋值为最大值
                            thisLenTmp = sum(ref_len_tmp_list)

                    # 是独立的比对
                    if synDataTmp1.iloc[0, 0] >= preLenTmp and synDataTmp1.iloc[0, 1] <= thisLenTmp:
                        synDataTmp1.iloc[0, 0] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                        synDataTmp1.iloc[0, 1] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                        try:
                            # 赋值
                            synpd = pd.concat([synpd, synDataTmp1], ignore_index=True)
                        except KeyError:
                            # 赋值
                            synpd = pd.concat([synpd, synDataTmp1], ignore_index=True)

                        if args.debug:
                            logger.error(f"\n{synDataTmp1}\n")

                        break
                    # 跨越了两个比对序列
                    else:
                        # 将当前序列的比对坐标提取并保存
                        synDataTmp2 = synDataTmp1.copy()
                        synDataTmp2.iloc[0, 0] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                        synDataTmp2.iloc[0, 1] = ref_end_tmp_list[idxTmp]
                        try:
                            # 赋值
                            synpd = pd.concat([synpd, synDataTmp2], ignore_index=True)
                        except KeyError:
                            # 赋值
                            synpd = pd.concat([synpd, synDataTmp2], ignore_index=True)
                        if args.debug:
                            logger.error(f"\n{synDataTmp2}\n")
                        # 更新synDataTmp1的坐标，计算剩下的坐标
                        synDataTmp1.iloc[0, 0] = thisLenTmp + 1
            except KeyError as e:  # 处理异常
                if args.debug:
                    logger.error(f"KeyError occurred: {e}")
                continue
            except IndexError as e:  # 处理异常
                if args.debug:
                    logger.error(f"IndexError occurred: {e}")
                    logger.error(f'Error ref_start_tmp_list: {ref_start_tmp_list}')
                    logger.error(f"Error ref_end_tmp_list: {ref_end_tmp_list}")
                    logger.error(f"Error ref_len_tmp_list: {ref_len_tmp_list}")
                    logger.error(f"Error index:{index}", index)
                    logger.error(f"Error alignment_file_name: {alignment_file_name}")
                continue
        
        return synpd


# 保存结果到文件中
def save_result(sample1, sample2, synpd):
    """
    :param sample1:     样品1名字
    :param sample2:     样品2名字
    :param synpd:      共线性的坐标   pd.DataFrame()
    :
    :return: 0
    """
    # 保存输出字符串
    outTxt = ""

    # 将结果写入输出字符串中
    for row in synpd.itertuples():  # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
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
    outputFilePath = os.path.abspath("{}_{}.syn.out".format(sample1, sample2))
    with open(outputFilePath, 'a', encoding="utf-8") as outputFile:
        outputFile.write(outTxt)

    return 0


# 保存结果到文件中
def sort_result(filePath):
    """
    :param filePath:    需要排序的文件路径
    :
    :return: 0
    """
    # 打开文件并读取每一行
    with open(filePath, 'r') as f:
        lines = f.readlines()
    
    # 按染色体和起始位置排序并输出
    sorted_lines = sorted(lines, key=lambda x: (x.split()[0], int(x.split()[1])))
    
    # 将排序后的行写入新的文件中
    with open(filePath, 'w') as f:
        f.writelines(sorted_lines)

    return 0


def main(args, config_file_map, no_synPath, workDir, code_path):
    """
    :param args:            getParser输出参数
    :param config_file_map:   getParser输出每个样品所有需要样品的路径
    :param no_synPath:      no_syn 输出的需要比对的非共线性坐标
    :param workDir:         工作路径
    :param code_path:       环境变量
    :
    :return: coorPath
    """
    # 创建工作路径
    makedir(workDir)
    workDir = workDir + os.sep

    # 切换工作路径
    os.chdir(workDir)

    
    # 构造类
    ali_loci_class = AlignmentLociClass(no_synPath, config_file_map, args, code_path)


    # #################################### 构建索引 #################################### #
    logger.error(f'Building index.')
    ali_loci_class._index()


    # #################################### Alignment #################################### #
    logger.error(f'Alignment started.')

    # 创建进程池，指定最大进程数为 args.jobs
    pool = multiprocessing.Pool(processes=args.jobs)

    # 任务计数器初始化
    task_count = 0
    task_total = len(ali_loci_class._long_ali_loci_dict) + len(ali_loci_class._short_ali_loci_dict)

    # 遍历 ali_loci_class._short_ali_loci_dict 字典，多进程
    if len(ali_loci_class._short_ali_loci_dict) > 0:
        results = pool.imap_unordered(ali_loci_class._alignment_line, ali_loci_class._short_ali_loci_dict.items())

        # #################################### Save #################################### #
        for loci, result in results:
            # 每个任务完成后，将任务计数器加1
            task_count += 1
            # 打印进度条
            progress = task_count / task_total * 100
            logger.error(f'Alignment Progress: {progress:.2f}% | {task_count}/{task_total} | {loci}')

            # 保存结果
            for key1, value1 in result.items():
                for key2, value2 in value1.items():
                    save_result(key1, key2, value2)

    # 遍历 ali_loci_class._long_ali_loci_dict 字典，多进程
    for key1, value1 in ali_loci_class._long_ali_loci_dict.items():  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >
        loci, result = ali_loci_class._alignment_column(key1, value1, pool)

        # 每个任务完成后，将任务计数器加1
        task_count += 1
        # 打印进度条
        progress = task_count / task_total * 100
        logger.error(f'Alignment Progress: {progress:.2f}% | {task_count}/{task_total} | {loci}')

        # #################################### Save #################################### #
        for key1, value1 in result.items():
            for key2, value2 in value1.items():
                save_result(key1, key2, value2)
    

    # #################################### Sort #################################### #
    import glob
    synOutFilePathList = glob.glob("*.syn.out")
    pool.map_async(sort_result, synOutFilePathList).get()

    # 关闭进程池，等待所有任务完成
    pool.close()
    pool.join()


    # #################################### 输出配置文件信息 #################################### #
    no_syn_alignmentPath = os.path.abspath(f"{args.prefix}no_syn_alignment.out")
    with open(no_syn_alignmentPath, 'w', encoding="utf-8") as f:
        for synOutFilePath in synOutFilePathList:
            fileNameList = os.path.basename(synOutFilePath).replace(".syn.out", "").split("_")
            f.write('\t'.join([fileNameList[0], fileNameList[1], os.path.abspath(synOutFilePath)]) + "\n")

    return no_syn_alignmentPath
