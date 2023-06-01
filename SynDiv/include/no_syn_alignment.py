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


# ������־��ʽ
logger = logging.getLogger('no_syn_alignment')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # ���������̨
handler.setFormatter(formatter)
logger.addHandler(handler)


# ����Ŀ¼
def makedir(
    path_dir: str
):
    """
    :param path_dir: ��Ҫ�������ļ���·��
    :return: 0
    """
    if os.path.isdir(path_dir):
        shutil.rmtree(path_dir)
        os.makedirs(path_dir)
        log = '[makedir] \'{}\' already exists, clear and recreate.'.format(path_dir)
        logger.error(log)
    else:
        os.makedirs(path_dir)


# ���ж��߳�����
def run_command(command, envPath):
    # �ύ����
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=envPath)

    # ��ȡ�˳�״̬
    exit_state = proc.wait()

    stdout, stderr = proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()

    # ��׼����ʹ������
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
    ��Ҫ�ȶԵķǹ���������
    """
    def __init__(
            self,
            ali_loci_filename: str, 
            config_file_map: dict, 
            args, 
            code_path
    ):
        """
        ���캯��
        :param ali_loci_filename: ��Ҫ�ȶԵ�����  no_syn���
        :param config_file_map: getParser���
        :param args: getParser���
        """
        # ��Ҫ�ȶԵ������ļ�
        self._ali_loci_filename = ali_loci_filename

        # �������ļ��ֵ�
        self._config_file_map = config_file_map  # dict{sample, {genome: "path", aligns: "path", syri_out: "path"}}

        # �������
        self._args = args

        # ��������
        self._path = code_path

        # ��Ҫ�ȶԵ�����
        self._long_ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >�� �����������
        self._short_ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >���������ٵ���

        # syntenic result
        self._sample1_sample2_synpd = {}  # map<sample1, map<sample2, pd.DataFrame()> >


    def _index(
        self
    ):
        """
        ������������
        """
        # ������Ҫ�ȶԵ�����
        ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >

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
                ali_loci_dict[key_tmp] = {}

                # �ӵ����п�ʼѭ��
                for idx1 in range(2, len(infos_list)):
                    info = infos_list[idx1]

                    info_list = info.strip().split(":")
                    sample_tmp = info_list[0]

                    loci_list = info_list[1].split(";")
                    
                    # ��ʼ���ֵ�
                    ali_loci_dict[key_tmp][sample_tmp] = []

                    for idx2 in range(len(loci_list)):
                        # ���ֻ��һ�����꣬����
                        if len(loci_list[idx2].split("-")) != 2:
                            continue

                        # ��¼����   map<sample, vector<chr:start-end>
                        ali_loci_dict[key_tmp][sample_tmp].append(chr_tmp + ":" + loci_list[idx2])

        # ʹ���ֵ��Ƶ�ʽ���������ʽ�� ali_loci_dict ���Ϊ�������Ͷ̵������ֵ�
        import math
        threshold_number = (-1 + math.sqrt(1 + 8 * max(1, self._args.jobs * 100))) / 2  # y = (x * (x + 1)) / 2    x = (-1 + sqrt(1 + 8 * y)) / 2
        self._long_ali_loci_dict = {k: v for k, v in ali_loci_dict.items() if len(v) > threshold_number}
        self._short_ali_loci_dict = {k: v for k, v in ali_loci_dict.items() if len(v) <= threshold_number}
        # ��ӡ����
        logger.error(f'Number of tasks running in parallel by column: {len(self._long_ali_loci_dict)}')
        logger.error(f'Number of tasks running in parallel by row: {len(self._short_ali_loci_dict)}')



    def _alignment_line(
        self, 
        key_value
    ):
        """
        ��ȡ���в��ȶԵĶ��̺߳���
        : param key_value   �ο�����������_map<sample, vector<chr:start-end> >   map<sample, vector<chr:start-end> >
        : param key1        �ο�����������, chr+"_"+start
        : param value1      map<sample, vector<chr:start-end> >
        :
        : return key1, sample1_sample2_synpd, map<sample1, map<sample2, pd.DataFrame()> >
        """
        key1, value1 = key_value
        # ��ʱ�洢����������
        sample1_sample2_synpd = {}  # map<sample1, map<sample2, pd.DataFrame()> >

        # �洢��ȡ�������ļ�
        file_name_map = {key2: os.path.abspath(f"{key1}_{key2}.fa") 
                    for key2 in value1.keys()}
        file_name_loc_map = {key2: {} for key2 in value1.keys()}
        sample_list = sorted(value1.keys())

        # Ⱦɫ���
        chromosome = ""

        # ############### subseq ############### #
        for key2, value2 in value1.items():  # map<sample, vector<chr:start-end> >
            # ��ȡ�����б���ĵط�
            output_file_name = file_name_map[key2]

            # �ύ����
            chromosome, file_name_loc_map_tmp = samtools(
                self._path, 
                self._ali_loci_filename,
                key1, 
                key2, 
                value2, 
                self._config_file_map, 
                output_file_name
            )

            # ��¼����Ʒ��������Ϣ
            for sample2_tmp, loci_map in file_name_loc_map_tmp.items():
                file_name_loc_map[sample2_tmp] = loci_map


        # ############### minimap2 + syri ############### #
        # ��¼��sample1�����Ե�sample�������ж�֮���sample֮���Ƿ���Ҫ�����໥�ȶ�
        syn_sample_dict = {}  # map<sample1, map<sample2, 0/1> >   0->noSyn  1->Syn

        # ��Ʒ������
        sample_list.sort()

        # ��Ʒ1
        for idx1, sample1 in enumerate(sample_list):
            filename1 = file_name_map[sample1]

            # ��¼��sample1�����Ե�sample�������ж�֮���sample֮���Ƿ���Ҫ�����໥�ȶ�
            synSampleList = []
            # ��¼��sample1�ǹ����Ե�sample�������ж�֮���sample֮���Ƿ���Ҫ�����໥�ȶ�
            noSynSampleList = []

            # ref����ͳ�����Ϣ
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
            
            # ��Ʒ2
            for idx2 in range(idx1 + 1, len(sample_list)):
                sample2 = sample_list[idx2]
                filename2 = file_name_map[sample2]

                # qry����ͳ�����Ϣ
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

                # �ж��Ƿ���Ҫ�ȶԺͼ��㹲����
                if sample1 in syn_sample_dict:
                    if sample2 in syn_sample_dict[sample1]:
                        # 1���� self._args.synRatio Ϊ�����ԣ�ֱ�����
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
                                # ��ֵ
                                sample1_sample2_synpd[sample1][sample2] = SynPDTmp
                            except KeyError:
                                # ��ʼ��
                                if sample1 not in sample1_sample2_synpd:
                                    sample1_sample2_synpd[sample1] = {}
                                if sample2 not in sample1_sample2_synpd[sample1]:
                                    sample1_sample2_synpd[sample1][sample2] = pd.DataFrame()
                                # ��ֵ
                                sample1_sample2_synpd[sample1][sample2] = SynPDTmp
                            continue
                        # 0����Ϊ�ǹ����ԣ�ͬ�����ȶԣ�ֱ��������sample2
                        else:
                            continue
                
                # ����ļ���
                alignment_file_name = os.path.abspath(f'{key1}_{sample1}_{sample2}.sam')

                # �ύ����
                # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
                synpd = minimap2_syn(self._args, self._path, filename1, filename2, alignment_file_name, ref_start_tmp_list, ref_end_tmp_list, ref_len_tmp_list)

                # ���Ϊ��ֵ��ֱ����һ��ѭ��
                if synpd.empty:
                    # ��¼Ϊ�ǹ����Ե�����
                    noSynSampleList.append(sample2)
                    continue

                try:
                    # ��ֵ
                    sample1_sample2_synpd[sample1][sample2] = synpd
                except KeyError:
                    # ��ʼ��
                    if sample1 not in sample1_sample2_synpd:
                        sample1_sample2_synpd[sample1] = {}
                    if sample2 not in sample1_sample2_synpd[sample1]:
                        sample1_sample2_synpd[sample1][sample2] = pd.DataFrame()
                    # ��ֵ
                    sample1_sample2_synpd[sample1][sample2] = synpd

                # ���������һ���б����ǹ����Ե���Ʒȫ���ռ������������ȶ�������Щ��Ʒ  Ҫ������ֵ �����Գ���ռ self._args.synRatio ����
                qryTotal = synpd['bLen'].sum()
                alignRatio = float(qryTotal)/max(qryLenSum, 1)
                if alignRatio >= self._args.synRatio:
                    synSampleList.append(sample2)
                elif alignRatio <= self._args.nosynRatio:
                    # ��¼Ϊ�ǹ����Ե�����
                    noSynSampleList.append(sample2)

            # �������Ե���Ʒ��ӵ��ܵ�dict�У�����֮��ȶԲ�ѯ
            # ��Ʒ1
            for idx2, sample2 in enumerate(synSampleList):
                # ��ʼ�� dict
                if sample2 not in syn_sample_dict:
                    syn_sample_dict[sample2] = {}
                # ��Ʒ2
                for idx3 in range(idx2 + 1, len(synSampleList)):
                    syn_sample_dict[sample2][synSampleList[idx3]] = 1
            # ���ǹ����Ե���Ʒ��ӵ��ܵ�dict�У�����֮��ȶԲ�ѯ
            # ��Ʒ1
            for sample2 in synSampleList:
                # ��ʼ�� dict
                if sample2 not in syn_sample_dict:
                    syn_sample_dict[sample2] = {}
                # ��Ʒ2
                for sample3 in noSynSampleList:
                    syn_sample_dict[sample2][sample3] = 0

        # remove file (fasta)
        for key2, value2 in file_name_map.items():
            if os.path.exists(value2) and self._args.debug==False:
                os.remove(value2)

        # ���ص�ǰλ�õ�������Ʒ��ϵĹ���������
        return key1, sample1_sample2_synpd
    

    def _alignment_column(
        self, 
        key1, 
        value1, 
        pool
    ):
        """
        ��ȡ���в��ȶԵĶ��̺߳���
        : param key1    �ο�����������, map<chr+"_"+start, map<sample, vector<chr:start-end> > >
        : param value1  map<sample, vector<chr:start-end> >
        : param pool    ���̳�
        :
        : return key1, sample1_sample2_synpd_tmp, map<sample1, map<sample2, pd.DataFrame()> >
        """
        # ��ʱ�洢����������
        sample1_sample2_synpd_tmp = {}  # map<sample1, map<sample2, pd.DataFrame()> >

        # �洢��ȡ�������ļ�
        file_name_map = {key2: os.path.abspath(f"{key1}_{key2}.fa") 
                    for key2 in value1.keys()}
        file_name_loc_map = {key2: {} for key2 in value1.keys()}
        sample_list = sorted(value1.keys())

        # ############### subseq ############### #
        results = []
        for key2, value2 in value1.items():  # map<sample, vector<chr:start-end> >
            # ��ȡ�����б���ĵط�
            output_file_name = file_name_map[key2]

            # �ύ����
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

        # ��ȡ�����ķ���ֵ
        for result in results:
            chromosome, file_name_loc_map_tmp = result.get()

            # ��¼����Ʒ��������Ϣ
            for sample2_tmp, loci_map in file_name_loc_map_tmp.items():
                file_name_loc_map[sample2_tmp] = loci_map


        # ############### minimap2 + syri ############### #
        # ��Ʒ������
        sample_list.sort()

        # �洢�ȶԵ����ս��
        results = []
        
        # ��Ʒ1
        for idx1, sample1 in enumerate(sample_list):
            filename1 = file_name_map[sample1]
            
            # ref����ͳ�����Ϣ
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
            
            # ��Ʒ2
            for idx2 in range(idx1 + 1, len(sample_list)):
                sample2 = sample_list[idx2]
                filename2 = file_name_map[sample2]
                
                # ����ļ���
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

        # ��ȡ�����ķ���ֵ
        for result in results:
            # ��ȡ������
            sample1, sample2, synpd = result
            synpd = synpd.get()
            
            # ���Ϊ��ֵ��ֱ����һ��ѭ��
            if synpd.empty:
                continue

            try:
                # ��ֵ
                sample1_sample2_synpd_tmp[sample1][sample2] = synpd
            except KeyError:
                # ��ʼ��
                if sample1 not in sample1_sample2_synpd_tmp:
                    sample1_sample2_synpd_tmp[sample1] = {}
                if sample2 not in sample1_sample2_synpd_tmp[sample1]:
                    sample1_sample2_synpd_tmp[sample1][sample2] = pd.DataFrame()
                # ��ֵ
                sample1_sample2_synpd_tmp[sample1][sample2] = synpd

        # remove file
        for key2, value2 in file_name_map.items():
            if os.path.exists(value2) and self._args.debug==False:
                os.remove(value2)

        # ���ص�ǰλ�õ�������Ʒ��ϵĹ���������
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
    :description:                 ��ȡ����
    :param path:                  ��������
    :param ali_loci_filename:     SynDic_c no_syn �������Ҫ�ȶԵ�����
    :param key1:                  chr+"_"+start
    :param key2:                  sample
    :param value2:                vector<chr:start-end>
    :param config_file_map:         ����������ļ�dict
    :param output_file_name:        ����ļ���
    :return file_name_loc_map_tmp:    map<sample2, map<start, end> >
    """
    # ########## ��ȡ���� ######### #
    # ��¼sample2��������Ϣ
    file_name_loc_map_tmp = {}
    file_name_loc_map_tmp[key2] = {}  # map<sample2, map<start, end> >

    for idx1, loci in enumerate(value2):
        # ���Ⱦɫ����Ƿ���Ϲ涨  (chr1:257130-257225)
        if ":" not in loci or "-" not in loci:
            raise ValueError(f"Error: incorrect '{loci}' found in '{ali_loci_filename}'.")
        
        # ��ȡȾɫ���
        chromosome = loci.split(":")[0]

        cmd = ""

        # ���Ϊ��һ�����У�����Ⱦɫ���
        if idx1 == 0:
            # SynDiv �ύ�������ļ��� ["genome"] ��
            try:
                cmd = f'echo \">{key1.split("_")[0]}\" > {output_file_name}; samtools faidx {config_file_map[key2]["genome"]} {loci} | grep -v \">\" >> {output_file_name}'
            # SynDiv_ �ύ�������ļ���û�� ["genome"] ��
            except TypeError:
                cmd = f'echo \">{key1.split("_")[0]}\" > {output_file_name}; samtools faidx {config_file_map[key2]} {loci} | grep -v \">\" >> {output_file_name}'
        else:  # ���Ϊ֮������У���ҪȾɫ���
            # SynDiv �ύ�������ļ��� ["genome"] ��
            try:
                cmd = f'samtools faidx {config_file_map[key2]["genome"]} {loci} | grep -v \">\" >> {output_file_name}'
            # SynDiv_ �ύ�������ļ���û�� ["genome"] ��
            except TypeError:
                cmd = f'samtools faidx {config_file_map[key2]} {loci} | grep -v \">\" >> {output_file_name}'

        # ��¼����
        lociList = loci.split(":")[1].split("-")
        file_name_loc_map_tmp[key2][int(lociList[0])] = int(lociList[1])

        # �ύ����
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
    :description:                 �ȶԲ�Ѱ�ҹ���������
    :param args:                  ����Ĳ�����Ϣ
    :param path:                  ��������
    :param filename1:             ��Ʒ1����
    :param filename2:             ��Ʒ2����
    :param alignment_file_name:     �ȶ��ļ����·��
    :param ref_start_tmp_list:       ��Ʒ1����ʼ�����б�
    :param ref_end_tmp_list:         ��Ʒ1����ֹ�����б�
    :param ref_len_tmp_list:         ��Ʒ1�ĳ����б�
    """
    # minimap2 �ȶ�
    cmd = f'minimap2 -ax asm5 --eqx {filename1} {filename2} -o {alignment_file_name}'
    # �ύ����
    run_command(cmd, path)

    # ����ļ������ڣ�ֱ�ӷ���
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
        :description:                           ����minimap2���, Ѱ�ҹ���������
        :param args:                            ����Ĳ�����Ϣ
        :param ref_start_tmp_list:                 ref����ͳ�����Ϣ
        :param ref_end_tmp_list:                   ref����ͳ�����Ϣ
        :param ref_len_tmp_list:                   ref����ͳ�����Ϣ
        :param alignment_file_name:               minimap2 �ȶԵĽ���ļ�
        :return: synpd                          synpd -> pd.DataFrame()
        """
        # ��ʱ�洢����������
        # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
        synpd = pd.DataFrame()

        # syri Ѱ�ҹ�������Ϣ
        synDatas = pd.DataFrame()
        try:
            synDatas = syri_syn.main(args, alignment_file_name)
        except:
            return synpd

        # remove file
        if os.path.exists(alignment_file_name) and args.debug==False:
            os.remove(alignment_file_name)

        # ���synDatasΪ�գ�ֱ�ӷ��ؿ�ֵ
        if synDatas.empty:
            return synpd
        
        # ���Դ���
        if args.debug:
            logger.error(f"alignment_file_name:{alignment_file_name}\nsynDatas:{synDatas}")

        # �������굽�ܱ���
        # ���д�ӡÿһ�е�Ԫ��
        # �����ж�������ô����
        idxTmp = 0  # ��ʱ����
        preLenTmp = 0  # ǰ�����еĳ���(n)��ͬ���ж��ǲ���һ�������ıȶԶ�
        thisLenTmp = ref_len_tmp_list[0]  # ǰ�����еĳ��ȼ�������һ������(n+1)�����ֵ�����±�û����
        for index, row in synDatas.iterrows():
            try:
                # ��ʱ�� DataFrame
                synDataTmp1 = synDatas.iloc[[index]].copy()

                if args.debug:
                    logger.error(f"synDataTmp1:{synDataTmp1}")

                # ���Ϊ�ռ����߷���ȶԣ���һ��ѭ��
                if synDataTmp1.empty or synDataTmp1.iloc[0, 7] == -1 or synDataTmp1.iloc[0, 8] == -1:
                    continue

                # �ж�����λ��
                while True and idxTmp < len(ref_len_tmp_list):
                    # �����ʼ���ڵ�ǰ�������еĳ��ȣ����������ͳ���
                    while synDataTmp1.iloc[0, 0] > thisLenTmp and idxTmp < len(ref_len_tmp_list):
                        idxTmp += 1
                        preLenTmp = sum(ref_len_tmp_list[:idxTmp])

                        # ����±߻������꣬������һ������
                        if idxTmp + 1 < len(ref_len_tmp_list):
                            thisLenTmp = sum(ref_len_tmp_list[:idxTmp + 1])
                        else:  # ���û�У���ֵΪ���ֵ
                            thisLenTmp = sum(ref_len_tmp_list)

                    # �Ƕ����ıȶ�
                    if synDataTmp1.iloc[0, 0] >= preLenTmp and synDataTmp1.iloc[0, 1] <= thisLenTmp:
                        synDataTmp1.iloc[0, 0] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                        synDataTmp1.iloc[0, 1] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                        try:
                            # ��ֵ
                            synpd = pd.concat([synpd, synDataTmp1], ignore_index=True)
                        except KeyError:
                            # ��ֵ
                            synpd = pd.concat([synpd, synDataTmp1], ignore_index=True)

                        if args.debug:
                            logger.error(f"\n{synDataTmp1}\n")

                        break
                    # ��Խ�������ȶ�����
                    else:
                        # ����ǰ���еıȶ�������ȡ������
                        synDataTmp2 = synDataTmp1.copy()
                        synDataTmp2.iloc[0, 0] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                        synDataTmp2.iloc[0, 1] = ref_end_tmp_list[idxTmp]
                        try:
                            # ��ֵ
                            synpd = pd.concat([synpd, synDataTmp2], ignore_index=True)
                        except KeyError:
                            # ��ֵ
                            synpd = pd.concat([synpd, synDataTmp2], ignore_index=True)
                        if args.debug:
                            logger.error(f"\n{synDataTmp2}\n")
                        # ����synDataTmp1�����꣬����ʣ�µ�����
                        synDataTmp1.iloc[0, 0] = thisLenTmp + 1
            except KeyError as e:  # �����쳣
                if args.debug:
                    logger.error(f"KeyError occurred: {e}")
                continue
            except IndexError as e:  # �����쳣
                if args.debug:
                    logger.error(f"IndexError occurred: {e}")
                    logger.error(f'Error ref_start_tmp_list: {ref_start_tmp_list}')
                    logger.error(f"Error ref_end_tmp_list: {ref_end_tmp_list}")
                    logger.error(f"Error ref_len_tmp_list: {ref_len_tmp_list}")
                    logger.error(f"Error index:{index}", index)
                    logger.error(f"Error alignment_file_name: {alignment_file_name}")
                continue
        
        return synpd


# ���������ļ���
def save_result(sample1, sample2, synpd):
    """
    :param sample1:     ��Ʒ1����
    :param sample2:     ��Ʒ2����
    :param synpd:      �����Ե�����   pd.DataFrame()
    :
    :return: 0
    """
    # ��������ַ���
    outTxt = ""

    # �����д������ַ�����
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

    # ������ļ�
    outputFilePath = os.path.abspath("{}_{}.syn.out".format(sample1, sample2))
    with open(outputFilePath, 'a', encoding="utf-8") as outputFile:
        outputFile.write(outTxt)

    return 0


# ���������ļ���
def sort_result(filePath):
    """
    :param filePath:    ��Ҫ������ļ�·��
    :
    :return: 0
    """
    # ���ļ�����ȡÿһ��
    with open(filePath, 'r') as f:
        lines = f.readlines()
    
    # ��Ⱦɫ�����ʼλ���������
    sorted_lines = sorted(lines, key=lambda x: (x.split()[0], int(x.split()[1])))
    
    # ����������д���µ��ļ���
    with open(filePath, 'w') as f:
        f.writelines(sorted_lines)

    return 0


def main(args, config_file_map, no_synPath, workDir, code_path):
    """
    :param args:            getParser�������
    :param config_file_map:   getParser���ÿ����Ʒ������Ҫ��Ʒ��·��
    :param no_synPath:      no_syn �������Ҫ�ȶԵķǹ���������
    :param workDir:         ����·��
    :param code_path:       ��������
    :
    :return: coorPath
    """
    # ��������·��
    makedir(workDir)
    workDir = workDir + os.sep

    # �л�����·��
    os.chdir(workDir)

    
    # ������
    ali_loci_class = AlignmentLociClass(no_synPath, config_file_map, args, code_path)


    # #################################### �������� #################################### #
    logger.error(f'Building index.')
    ali_loci_class._index()


    # #################################### Alignment #################################### #
    logger.error(f'Alignment started.')

    # �������̳أ�ָ����������Ϊ args.jobs
    pool = multiprocessing.Pool(processes=args.jobs)

    # �����������ʼ��
    task_count = 0
    task_total = len(ali_loci_class._long_ali_loci_dict) + len(ali_loci_class._short_ali_loci_dict)

    # ���� ali_loci_class._short_ali_loci_dict �ֵ䣬�����
    if len(ali_loci_class._short_ali_loci_dict) > 0:
        results = pool.imap_unordered(ali_loci_class._alignment_line, ali_loci_class._short_ali_loci_dict.items())

        # #################################### Save #################################### #
        for loci, result in results:
            # ÿ��������ɺ󣬽������������1
            task_count += 1
            # ��ӡ������
            progress = task_count / task_total * 100
            logger.error(f'Alignment Progress: {progress:.2f}% | {task_count}/{task_total} | {loci}')

            # ������
            for key1, value1 in result.items():
                for key2, value2 in value1.items():
                    save_result(key1, key2, value2)

    # ���� ali_loci_class._long_ali_loci_dict �ֵ䣬�����
    for key1, value1 in ali_loci_class._long_ali_loci_dict.items():  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >
        loci, result = ali_loci_class._alignment_column(key1, value1, pool)

        # ÿ��������ɺ󣬽������������1
        task_count += 1
        # ��ӡ������
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

    # �رս��̳أ��ȴ������������
    pool.close()
    pool.join()


    # #################################### ��������ļ���Ϣ #################################### #
    no_syn_alignmentPath = os.path.abspath(f"{args.prefix}no_syn_alignment.out")
    with open(no_syn_alignmentPath, 'w', encoding="utf-8") as f:
        for synOutFilePath in synOutFilePathList:
            fileNameList = os.path.basename(synOutFilePath).replace(".syn.out", "").split("_")
            f.write('\t'.join([fileNameList[0], fileNameList[1], os.path.abspath(synOutFilePath)]) + "\n")

    return no_syn_alignmentPath
