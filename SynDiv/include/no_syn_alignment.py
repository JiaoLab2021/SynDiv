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


# ����Ŀ¼
def makedir(
    path_dir: str
):
    """
    :param path_dir: ��Ҫ�������ļ���·��
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
    ��Ҫ�ȶԵķǹ���������
    """
    def __init__(
            self,
            ali_loci_filename: str, 
            configFileMap: dict, 
            args, 
            code_path
    ):
        """
        ���캯��
        :param ali_loci_filename: ��Ҫ�ȶԵ�����  no_syn���
        :param configFileMap: getParser���
        :param args: getParser���
        """

        # ��Ҫ�ȶԵ������ļ�
        self._ali_loci_filename = ali_loci_filename

        # �������ļ��ֵ�
        self._configFileMap = configFileMap  # dict{sample, {genome: "path", aligns: "path", syri_out: "path"}}

        # �������
        self._args = args

        # ��������
        self._path = code_path

        # �������
        self._ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >

        # ��������
        self._index()

        # syntenic result
        self._sample1Sample2SynPD = {}  # map<sample1, map<sample2, pd.DataFrame()> >


    def _index(
        self
    ):
        """
        ������������
        """
        with open(self._ali_loci_filename) as f:
            for infos in f.readlines():
                # ������������
                if len(infos) == 0:
                    continue

                # ȥ�ַ��������
                infos_list = infos.strip().split()
                
                # ���ֻ��һ�����꣬��������
                if len(infos_list) <= 3:
                    continue

                # Ⱦɫ���
                chr_tmp = infos_list[0]

                # ��ʱ��
                key_tmp = infos_list[0] + "_" + infos_list[1]

                # ��ʼ���ֵ�
                self._ali_loci_dict[key_tmp] = {}

                # �ӵ����п�ʼѭ��
                for idx1 in range(2, len(infos_list)):
                    info = infos_list[idx1]

                    info_list = info.strip().split(":")
                    sample_tmp = info_list[0]

                    loci_list = info_list[1].split(";")
                    
                    # ��ʼ���ֵ�
                    self._ali_loci_dict[key_tmp][sample_tmp] = []

                    for idx2 in range(len(loci_list)):
                        # ���ֻ��һ�����꣬����
                        if len(loci_list[idx2].split("-")) != 2:
                            continue

                        # ��¼����   map<sample, vector<chr:start-end>
                        self._ali_loci_dict[key_tmp][sample_tmp].append(chr_tmp + ":" + loci_list[idx2])

    
    def _alignment(
        self, 
        key1, 
        value1, 
        threads, 
        debug
    ):
        """
        ��ȡ���в��ȶԵĶ��̺߳���
        : param key1  �ο�����������, map<chr+"_"+start, map<sample, vector<chr:start-end> > >
        : param value1  map<sample, vector<chr:start-end> >
        : param threads  �߳���
        : param debug    �Ƿ���Դ���
        : return sample1Sample2SynPDTmp, map<sample1, map<sample2, pd.DataFrame()> >
        """
        logger = logging.getLogger('_alignment')

        # ��ʱ�洢����������
        sample1Sample2SynPDTmp = {}  # map<sample1, map<sample2, pd.DataFrame()> >

        # �洢��ȡ�������ļ�
        fileNameMap = {key2: os.path.abspath(f"{key1}_{key2}.fa") 
                    for key2 in value1.keys()}
        fileNameLocMap = {key2: {} for key2 in value1.keys()}
        sampleList = sorted(value1.keys())

        # ��ȡ����
        for key2, value2 in value1.items():  # map<sample, vector<chr:start-end> >
            # ��ȡ�����б���ĵط�
            outputFileName = fileNameMap[key2]

            # ########## ��ȡ���� ######### #
            for idx1, loci in enumerate(value2):
                # ���Ⱦɫ����Ƿ���Ϲ涨  (chr1:257130-257225)
                if ":" not in loci or "-" not in loci:
                    raise ValueError(f"Error: incorrect '{loci}' found in '{self._ali_loci_filename}'.")
                    
                cmd = ""

                # ���Ϊ��һ�����У�����Ⱦɫ���
                if idx1 == 0:
                    # SynDiv �ύ�������ļ��� ["genome"] ��
                    try:
                        cmd = f'echo \">{key1.split("_")[0]}\" > {outputFileName}; samtools faidx {self._configFileMap[key2]["genome"]} {loci} | grep -v \">\" >> {outputFileName}'
                    # SynDiv_ �ύ�������ļ���û�� ["genome"] ��
                    except TypeError:
                        cmd = f'echo \">{key1.split("_")[0]}\" > {outputFileName}; samtools faidx {self._configFileMap[key2]} {loci} | grep -v \">\" >> {outputFileName}'
                else:  # ���Ϊ֮������У���ҪȾɫ���
                    # SynDiv �ύ�������ļ��� ["genome"] ��
                    try:
                        cmd = f'samtools faidx {self._configFileMap[key2]["genome"]} {loci} | grep -v \">\" >> {outputFileName}'
                    # SynDiv_ �ύ�������ļ���û�� ["genome"] ��
                    except TypeError:
                        cmd = f'echo \">{key1.split("_")[0]}\" > {outputFileName}; samtools faidx {self._configFileMap[key2]} {loci} | grep -v \">\" >> {outputFileName}'

                # �ύ����
                cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=self._path)
                # ��׼����ʹ������
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

                # ��¼����
                lociList = loci.split(":")[1].split("-")
                fileNameLocMap[key2][int(lociList[0])] = int(lociList[1])


        # ############### alignment ############### #
        # ��Ʒ������
        sampleList.sort()

        # ��Ʒ1
        for idx1, sample1 in enumerate(sampleList):
            filename1 = fileNameMap[sample1]

            # ref����ͳ�����Ϣ
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
            
            # ��Ʒ2
            for idx2 in range(idx1 + 1, len(sampleList)):
                sample2 = sampleList[idx2]
                filename2 = fileNameMap[sample2]
                
                # ����ļ���
                alignMentFileName = os.path.abspath(f'{key1}_{sample1}_{sample2}.sam')

                # minimap2 �ȶ�
                cmd = f'minimap2 -t {threads} -ax asm5 --eqx {filename1} {filename2} -o {alignMentFileName}'
                # �ύ����
                cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=self._path)
                # ��׼����ʹ������
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

                # syri Ѱ�ҹ���������
                # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
                synDatas = syri_syn.main(self._args, alignMentFileName)

                # remove file
                if os.path.exists(alignMentFileName) and debug==False:
                    os.remove(alignMentFileName)

                # ���synDatasΪ�գ���һ��ѭ��
                if synDatas.empty:
                    continue

                # ���Դ���
                if debug:
                    logger.error(f"alignMentFileName:{alignMentFileName}\nsynDatas:{synDatas}")

                # �������굽�ܱ���
                # ���д�ӡÿһ�е�Ԫ��
                # �����ж�������ô����
                idxTmp = 0  # ��ʱ����
                preLenTmp = 0  # ǰ�����еĳ���(n)��ͬ���ж��ǲ���һ�������ıȶԶ�
                thisLenTmp = refLenTmpList[0]  # ǰ�����еĳ��ȼ�������һ������(n+1)�����ֵ�����±�û����
                for index, row in synDatas.iterrows():
                    try:
                        # ��ʱ�� DataFrame
                        synDataTmp1 = synDatas.iloc[[index]].copy()

                        if debug:
                            logger.error(f"synDataTmp1:{synDataTmp1}")

                        # ���Ϊ�ռ����߷���ȶԣ���һ��ѭ��
                        if synDataTmp1.empty or synDataTmp1.iloc[0, 7] == "-" or synDataTmp1.iloc[0, 8] == "-":
                            continue

                        # �ж�����λ��
                        while True and idxTmp < len(refLenTmpList):
                            # �����ʼ���ڵ�ǰ�������еĳ��ȣ����������ͳ���
                            while synDataTmp1.iloc[0, 0] > thisLenTmp and idxTmp < len(refLenTmpList):
                                idxTmp += 1
                                preLenTmp = sum(refLenTmpList[:idxTmp])

                                # ����±߻������꣬������һ������
                                if idxTmp + 1 < len(refLenTmpList):
                                    thisLenTmp = sum(refLenTmpList[:idxTmp + 1])
                                else:  # ���û�У���ֵΪ���ֵ
                                    thisLenTmp = sum(refLenTmpList)

                            # �Ƕ����ıȶ�
                            if synDataTmp1.iloc[0, 0] >= preLenTmp and synDataTmp1.iloc[0, 1] <= thisLenTmp:
                                synDataTmp1.iloc[0, 0] += refStartTmpList[idxTmp] - 1 - preLenTmp
                                synDataTmp1.iloc[0, 1] += refStartTmpList[idxTmp] - 1 - preLenTmp
                                try:
                                    # ��ֵ
                                    sample1Sample2SynPDTmp[sample1][sample2] = pd.concat([sample1Sample2SynPDTmp[sample1][sample2], synDataTmp1], ignore_index=True)
                                except KeyError:
                                    # ��ʼ��
                                    if sample1 not in sample1Sample2SynPDTmp:
                                        sample1Sample2SynPDTmp[sample1] = {}
                                    if sample2 not in sample1Sample2SynPDTmp[sample1]:
                                        sample1Sample2SynPDTmp[sample1][sample2] = pd.DataFrame()
                                    # ��ֵ
                                    sample1Sample2SynPDTmp[sample1][sample2] = pd.concat([sample1Sample2SynPDTmp[sample1][sample2], synDataTmp1], ignore_index=True)

                                if debug:
                                    logger.error(f"\n{synDataTmp1}\n")

                                break
                            # ��Խ�������ȶ�����
                            else:
                                # ����ǰ���еıȶ�������ȡ������
                                synDataTmp2 = synDataTmp1.copy()
                                synDataTmp2.iloc[0, 0] += refStartTmpList[idxTmp] - 1 - preLenTmp
                                synDataTmp2.iloc[0, 1] = refEndTmpList[idxTmp]
                                try:
                                    # ��ֵ
                                    sample1Sample2SynPDTmp[sample1][sample2] = pd.concat([sample1Sample2SynPDTmp[sample1][sample2], synDataTmp2], ignore_index=True)
                                except KeyError:
                                    # ��ʼ��
                                    if sample1 not in sample1Sample2SynPDTmp:
                                        sample1Sample2SynPDTmp[sample1] = {}
                                    if sample2 not in sample1Sample2SynPDTmp[sample1]:
                                        sample1Sample2SynPDTmp[sample1][sample2] = pd.DataFrame()
                                    # ��ֵ
                                    sample1Sample2SynPDTmp[sample1][sample2] = pd.concat([sample1Sample2SynPDTmp[sample1][sample2], synDataTmp2], ignore_index=True)
                                if debug:
                                    logger.error(f"\n{synDataTmp2}\n")
                                # ����synDataTmp1�����꣬����ʣ�µ�����
                                synDataTmp1.iloc[0, 0] = thisLenTmp + 1
                    except KeyError as e:  # �����쳣
                        logger.error(f"KeyError occurred:{e}")
                        continue
                    except IndexError as e:  # �����쳣
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

        # ���ص�ǰλ�õ�������Ʒ��ϵĹ���������
        return sample1Sample2SynPDTmp


def main(args, configFileMap, no_synPath, workDir, code_path):
    """
    :param args:  getParser�������
    :param configFileMap:   getParser���ÿ����Ʒ������Ҫ��Ʒ��·��
    :param no_synPath:  Nexus_c no_syn ����ķǹ�������qurey�ϵ�����
    :param workDir:  ����·��
    :param code_path: ��������
    :return: coorPath
    """
    # ��������·��
    makedir(workDir)
    workDir = workDir + os.sep
    # �л�����·��
    os.chdir(workDir)

    logger = logging.getLogger('no_syn_alignment')

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Running.')

    # ������
    ali_loci_class = AlignmentLociClass(no_synPath, configFileMap, args, code_path)
    # ��������
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Building index.')
    ali_loci_class._index()

    # ���߳̽��
    output = []

    # �̳߳��ύ����
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Alignment.')
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        # �洢���е� future ����
        results = []
        # ���� ali_loci_class._ali_loci_dict �ֵ䣬�����߳�ִ�� ali_loci_class._alignmentRun ����
        for key1, value1 in ali_loci_class._ali_loci_dict.items():  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >
            future = executor.submit(ali_loci_class._alignment, key1, value1, args.threads, args.debug)
            results.append(future)
            # ���ƶ��߳�ִ�м��
            time.sleep(0.00005)

        # �ȴ����� future ����ִ����ϲ���ȡ���
        for future in as_completed(results):
            result = future.result()
            output.append(result)

    # �ϲ����
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Merge result.')
    for sample1Sample2SynPDTmp in output:
        for key2, value2 in sample1Sample2SynPDTmp.items():  # map<sample1, map<sample2, pd.DataFrame()> >
            for key3, value3 in value2.items():  # map<sample2, pd.DataFrame()>
                try:
                    # ��ֵ
                    ali_loci_class._sample1Sample2SynPD[key2][key3] = pd.concat([ali_loci_class._sample1Sample2SynPD[key2][key3], value3], ignore_index=True)
                except KeyError:
                    # ��ʼ��
                    if key2 not in ali_loci_class._sample1Sample2SynPD:
                        ali_loci_class._sample1Sample2SynPD[key2] = {}
                    if key3 not in ali_loci_class._sample1Sample2SynPD[key2]:
                        ali_loci_class._sample1Sample2SynPD[key2][key3] = pd.DataFrame()
                    # ��ֵ
                    ali_loci_class._sample1Sample2SynPD[key2][key3] = pd.concat([ali_loci_class._sample1Sample2SynPD[key2][key3], value3], ignore_index=True)


    # �� _sample1Sample2SynPD �����������
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Sort and save the results.')
    # All output file information
    allSyriOutInfo = ""
    for key1, value1 in ali_loci_class._sample1Sample2SynPD.items():  # map<sample1, map<sample2, pd.DataFrame()> >
        for key2, value2 in value1.items():  # map<sample2, pd.DataFrame()>
            # pd.DataFrame() ����
            value2.sort_values(by=['aChr', 'aStart', 'aEnd'], ascending=[True, True, True], inplace=True)

            # ��������ַ���
            outTxt = ""

            # �����д������ַ�����
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

            # ������ļ�
            outputFilePath = os.path.abspath("{}_{}.syn.out".format(key1, key2))
            with open(outputFilePath, 'w', encoding="utf-8") as outputFile:
                outputFile.write(outTxt)

            # ����ļ�����Ϣ
            allSyriOutInfo += '\t'.join([key1, key2, outputFilePath]) + "\n"

    # ����ļ���Ϣ
    no_syn_alignmentPath = os.path.join(os.getcwd() + os.sep, f"{args.prefix}no_syn_alignment.out")
    with open(no_syn_alignmentPath, 'w', encoding="utf-8") as f:
        f.write(allSyriOutInfo)

    # �ر��̳߳�
    executor.shutdown(wait=True)

    return no_syn_alignmentPath
