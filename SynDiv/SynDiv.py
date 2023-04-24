#!/home/ltan/anaconda3/envs/syri_env/bin/python3
# coding=gb2312
# __name__ = "SynDiv"

__data__ = "2023/04/24"
__version__ = "1.0.1"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
import datetime
import argparse
import sys
import shutil
import subprocess
import logging

# include
include_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
sys.path.append(os.path.join(include_dir, "include/"))
import no_syn_alignment

code_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep

# ��������
code_path = {'PATH': os.environ.get('PATH')}


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


# ��ȡfastaÿ�����еĳ���
def getLength(fastaFilePath, outputFilePath):
    """
    :param fastaFilePath: ��Ҫ�����fasta�ļ�·��
    :param outputFilePath:  ����������ļ�
    :return: 0
    """
    with open(fastaFilePath, "r") as f:
        seq = ""
        name = ""
        length_dict = {}
        for line in f:
            if line.startswith('>'):
                if seq != "":
                    # ͳ�����г��Ȳ����浽�ֵ���
                    length_dict[name] = len(seq)
                    seq = ""
                name = line.strip()[1:]
            else:
                seq += line.strip()
        # �������һ������
        length_dict[name] = len(seq)

    # ���ÿ�����е�ID�ͳ���
    with open(outputFilePath, "w") as f:
        for seq_id, seq_len in length_dict.items():
            f.write(f"{seq_id}\t{seq_len}\n")

    return 0


# ��������
def getParser():
    # log
    logger = logging.getLogger('getParser')

    logger.error(f"\ndata: {__data__}")
    logger.error(f"version: {__version__}")
    logger.error(f"author: {__author__}")
    logger.error(f"\nIf you encounter any issues related to the code, please don't hesitate to contact us via email at {__email__}.\n")

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("Input Files")
    required.add_argument("-r, --reference", dest="reference", help="Input FASTA reference", type=argparse.FileType('r'), required=True)
    required.add_argument("-c, --config", dest="config", help="Configuration file (format: sample\tgenome_path\taligns_path\tsyri_out_path)", type=argparse.FileType('r'), required=True)

    other = parser.add_argument_group("Additional arguments")
    other.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="S", choices=['T', 'S', 'B', 'P'])
    other.add_argument('-f', dest='f', help='Filter out low quality and small alignments. Use this parameter to use the full list of alignments without any filtering.', default=True, action='store_false')
    other.add_argument('--dir', dest='dir', help="Path to working directory (if not current directory). All files must be in this directory.", action='store')
    other.add_argument("--prefix", dest="prefix", help="Prefix to add before the output file Names", type=str, default="SynDiv.")
    other.add_argument("--no-chrmatch", dest='chrmatch', help="Don't allow automatic matching chromosome ids between the two genomes if they are not equal", default=False, action='store_true')
    other.add_argument('--mode', dest="mode", help="enabling quick mode will increase memory consumption", default="normal", choices=["fast", "normal"])
    other.add_argument('-w, --window', dest="window", help="window size", type=int, default=5000)
    other.add_argument('-s, --step', dest="step", help="step size", type=int, default=1000)
    other.add_argument("--debug", dest='debug', help="debug code", default=False, action='store_true')

    # ���в���
    parallel_option = parser.add_argument_group(title="parallel")
    # ����
    parallel_option.add_argument("--jobs", dest="jobs", help="Run n jobs in parallel", type=int, default=100)
    # �߳�
    parallel_option.add_argument("--threads", dest="threads", help="Number of threads used per job", type=int, default=3)

    args = parser.parse_args()

    # ���òο�������Ϊȫ·��
    referencePath = os.path.abspath(args.reference.name)

    # Set CWD and check if it exists
    if args.dir is None:
        args.dir = os.getcwd() + os.sep
    else:
        if os.path.isdir(args.dir):
            args.dir = os.path.abspath(args.dir + os.sep)
            # �л�����·��
            os.chdir(os.path.abspath(args.dir))
        else:
            logger.error(args.dir + ' is not a valid folder. Exiting.')
            sys.exit(1)

    # set the dict of all genomes path
    configFileMap = {}  # dict{sample, {genome: "path", aligns: "path", syri_out: "path"}}
    with open(args.config.name, "r") as f:
        for info in f.readlines():
            # ��������
            if len(info) == 0 or "#" in info:
                continue

            info = info.strip()
            infoList = info.split()

            # ��������Ƿ���Ϲ涨
            if len(infoList) != 4:
                logger.error("Error: '" + args.genomes.name + f"' is not four columns -> {info}")
                sys.exit(1)

            # ��¼����Ʒ�Ļ�����ȫ·��
            genomeFilePath = os.path.abspath(infoList[1])
            alignsFilePath = os.path.abspath(infoList[2])
            syriOutFilePath = os.path.abspath(infoList[3])
            configFileMap[infoList[0]] = {
                "genome": genomeFilePath, 
                "aligns": alignsFilePath, 
                "syri_out": syriOutFilePath
            }
            # ����ļ��Ƿ�����Ҵ���0
            for filePath in [genomeFilePath, alignsFilePath, syriOutFilePath]:
                if os.path.isfile(filePath) and os.path.getsize(filePath) > 0:
                    pass
                else:
                    logger.error(f"Error: '{filePath}': File does not exist or is empty.")
                    sys.exit(1)

    # ������Դ��룬���̳�����Ϊ1
    if args.debug:
        args.jobs = 1

    return args, referencePath, configFileMap


# ���� SynDiv_c multiinter
def multiinter(configFileMap, prefix, workDir):
    """
    :param configFileMap:   getParser���ÿ����Ʒ������Ҫ��Ʒ��·��
    :param prefix:  ���ǰ׺
    :param workDir:  ����·��
    :return: multiinterPath
    """
    logger = logging.getLogger('multiinter')

    # ��������·��
    makedir(workDir)
    workDir = workDir + os.sep
    # �л�����·��
    os.chdir(workDir)

    # multiinter ����
    cmd = f"{code_dir}SynDiv_c multiinter "
    filePathList = []
    sampleNameList = []

    for sampleName, value in configFileMap.items():
        filePathList.append(value["syri_out"])
        sampleNameList.append(sampleName)
    
    multiinterPath = os.path.join(workDir, f"{prefix}multiinter.out")
    cmd += "-i " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {multiinterPath}"

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] CMD: {cmd}')

    # �ύ����
    cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=code_path)
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

    return multiinterPath


# ���� SynDiv_c coor
def coor(configFileMap, multiinterPath, threads, prefix, workDir):
    """
    :param configFileMap:   getParser���ÿ����Ʒ������Ҫ��Ʒ��·��
    :param multiinterPath:  SynDiv_c multiinter ����Ĺ������ཻ���
    :param threads:  �߳���
    :param prefix:  ���ǰ׺
    :param workDir:  ����·��
    :return: coorPath
    """
    logger = logging.getLogger('coor')

    # ��������·��
    makedir(workDir)
    workDir = workDir + os.sep
    # �л�����·��
    os.chdir(workDir)

    # coor ����
    cmd = f"{code_dir}SynDiv_c coor -s {multiinterPath} -t {threads} "
    filePathList = []
    sampleNameList = []

    for sampleName, value in configFileMap.items():
        filePathList.append(value["aligns"])
        sampleNameList.append(sampleName)
    
    coorPath = os.path.join(workDir, f"{prefix}coor.out")
    cmd += "-i " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {coorPath}"

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] CMD: {cmd}')

    # �ύ����
    cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=code_path)
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

    return coorPath


# ���� SynDiv_c no_syn
def no_syn(configFileMap, coorPath, prefix, workDir):
    """
    :param configFileMap:   getParser���ÿ����Ʒ������Ҫ��Ʒ��·��
    :param coorPath:  SynDiv_c coor����Ĺ�������qurey�ϵ�����
    :param prefix:  ���ǰ׺
    :param workDir:  ����·��
    :return: coorPath
    """
    logger = logging.getLogger('no_syn')

    # ��������·��
    makedir(workDir)
    workDir = workDir + os.sep
    # �л�����·��
    os.chdir(workDir)

    # ����Ⱦɫ�峤����Ϣ
    filePathList = []
    sampleNameList = []
    for key, value in configFileMap.items():
        sampleNameList.append(key)
        filePath = os.path.join(workDir, f"{key}.length")
        filePathList.append(filePath)
        # ��ȡȾɫ�峤���ļ�
        getLength(value["genome"], filePath)

    # no_syn ����
    cmd = f"{code_dir}SynDiv_c no_syn --coor {coorPath} "
    
    no_synPath = os.path.join(workDir, f"{prefix}coor.out")
    cmd += "--lengths " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {no_synPath}"

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] CMD: {cmd}')

    # �ύ����
    cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=code_path)
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

    return no_synPath


# ���� SynDiv_c cal
def cal(args, referencePath, configFileMap, coorPath, no_syn_alignmentPath, workDir):
    """
    :param args:  getParser�������
    :param referencePath: �ο�������·��
    :param configFileMap:   getParser���ÿ����Ʒ������Ҫ��Ʒ��·��
    :param coorPath:  SynDiv_c coor����Ĺ�������qurey�ϵ�����
    :param no_syn_alignmentPath: no_syn_alignment ����������ļ�·��
    :param workDir:  ����·��
    :return: calPath
    """
    logger = logging.getLogger('cal')

    # ��������·��
    makedir(workDir)
    workDir = workDir + os.sep
    # �л�����·��
    os.chdir(workDir)

    # ������Ϣ
    filePathList = []
    sampleNameList = []
    for sampleName, value in configFileMap.items():
        filePathList.append(value["aligns"])
        sampleNameList.append(sampleName)

    # no_syn ����
    cmd = f"{code_dir}SynDiv_c cal -t {args.threads} -r {referencePath} --coor {coorPath} --syri_outs {no_syn_alignmentPath} "

    if args.mode == "fast":
        cmd += "--fast "
    
    calPath = os.path.join(workDir, f"{args.prefix}coor.out")
    cmd += "--aligns " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {calPath}"

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] CMD: {cmd}')

    # �ύ����
    cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=code_path)
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

    return calPath


# ���� SynDiv_c cal
def window(args, referencePath, calPath):
    """
    :param args:  getParser�������
    :param referencePath: �ο�������·��
    :param calPath:   SynDiv_c cal �����ÿ������� syntenic diversity
    :return: windowPath
    """
    logger = logging.getLogger('window')

    windowPath = os.path.abspath(f"{args.prefix}SynDiv.out")

    # no_syn ����
    cmd = f"{code_dir}SynDiv_c window -r {referencePath} -i {calPath} -w {args.window} -s {args.step} -o {windowPath}"

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] CMD: {cmd}')

    # �ύ����
    cmd_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=code_path)
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

    return windowPath


def main():
    logger = logging.getLogger('SynDiv')

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Running.')

    # ################################################# �������� ################################################# #
    args, referencePath, configFileMap = getParser()

    # ################################################# SynDiv_c multiinter ################################################# #
    os.chdir(args.dir)
    multiinterPath = multiinter(configFileMap, args.prefix, os.path.join(args.dir, "multiinter"))

    # ################################################# SynDiv_c coor ################################################# #
    os.chdir(args.dir)
    coorPath = coor(configFileMap, multiinterPath, args.jobs, args.prefix, os.path.join(args.dir, "coor"))

    # ################################################# SynDiv_c no_syn ################################################# #
    os.chdir(args.dir)
    no_synPath = no_syn(configFileMap, coorPath, args.prefix, os.path.join(args.dir, "no_syn"))

    # ################################################# no_syn_alignment ################################################# #
    os.chdir(args.dir)
    no_syn_alignmentPath = no_syn_alignment.main(args, configFileMap, no_synPath, os.path.join(args.dir, "no_syn_alignment"), code_path)

    # ################################################# no_syn_alignment ################################################# #
    os.chdir(args.dir)
    calPath = cal(args, referencePath, configFileMap, coorPath, no_syn_alignmentPath, os.path.join(args.dir, "cal"))
    
    # ################################################# no_syn_alignment ################################################# #
    os.chdir(args.dir)
    winPath = window(args, referencePath, calPath)

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Result: {calPath}')
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Result: {winPath}')
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Done.')


if __name__ == '__main__':
    main()
