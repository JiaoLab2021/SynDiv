#!/home/ltan/anaconda3/envs/syri_env/bin/python3
# coding=gb2312
# __name__ = "SynDiv"

__data__ = "2023/04/24"
__version__ = "1.0.1"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
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

# 环境变量
code_path = {'PATH': os.environ.get('PATH')}

# log
logger = logging.getLogger('SynDiv')
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
        sys.exit(1)

    return exit_state, stdout, stderr


# 获取fasta每条序列的长度
def getLength(fastaFilePath, outputFilePath):
    """
    :param fastaFilePath: 需要计算的fasta文件路径
    :param outputFilePath:  长度输出到文件
    :return: 0
    """
    with open(fastaFilePath, "r") as f:
        seq = ""
        name = ""
        length_dict = {}
        for line in f:
            if line.startswith('>'):
                if seq != "":
                    # 统计序列长度并保存到字典中
                    length_dict[name] = len(seq)
                    seq = ""
                name = line.strip()[1:]
            else:
                seq += line.strip()
        # 处理最后一条序列
        length_dict[name] = len(seq)

    # 输出每条序列的ID和长度
    with open(outputFilePath, "w") as f:
        for seq_id, seq_len in length_dict.items():
            f.write(f"{seq_id}\t{seq_len}\n")

    return 0


# 解析参数
def getParser():
    # log
    logger = logging.getLogger('getParser')

    logger.error(f"data: {__data__}")
    logger.error(f"version: {__version__}")
    logger.error(f"author: {__author__}")
    logger.error(f"\nIf you encounter any issues related to the code, please don't hesitate to contact us via email at {__email__}.\n")

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("Input Files")
    required.add_argument("-r, --reference", dest="reference", help="Input FASTA reference", type=str, required=True)
    required.add_argument("-c, --config", dest="config", help="Configuration file (format: sample\tgenome_path\taligns_path\tsyri_out_path)", type=str, required=True)

    other = parser.add_argument_group("Additional arguments")
    other.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="S", choices=['T', 'S', 'B', 'P'])
    other.add_argument('-f', dest='f', help='Filter out low quality and small alignments. Use this parameter to use the full list of alignments without any filtering.', default=True, action='store_false')
    other.add_argument('--dir', dest='dir', help="Path to working directory (if not current directory). All files must be in this directory.", action='store')
    other.add_argument("--prefix", dest="prefix", help="Prefix to add before the output file Names", type=str, default="SynDiv.")
    other.add_argument("--no-chrmatch", dest='chrmatch', help="Don't allow automatic matching chromosome ids between the two genomes if they are not equal", default=False, action='store_true')
    other.add_argument('--mode', dest="mode", help="Enabling quick mode will increase memory consumption", default="normal", choices=["fast", "normal"])
    other.add_argument('--synRatio', dest="synRatio", help="Threshold for complete synteny detection. Lower values increase alignment speed (0,1].", type=float, default=0.8)
    other.add_argument('--nosynRatio', dest="nosynRatio", help="Threshold for complete no-synteny detection. Higher values increase alignment speed (0,1].", type=float, default=0.05)
    other.add_argument('--window', dest="window", help="Window size", type=int, default=5000)
    other.add_argument('--step', dest="step", help="Step size", type=int, default=1000)
    other.add_argument("--debug", dest='debug', help="Debug code", default=False, action='store_true')

    # 并行参数
    parallel_option = parser.add_argument_group(title="parallel")
    # 进程
    parallel_option.add_argument("--jobs", dest="jobs", help="Run n jobs in parallel", type=int, default=150)
    # 线程
    parallel_option.add_argument("--threads", dest="threads", help="Number of threads used by SynDiv_c", type=int, default=10)
    
    args = parser.parse_args()

    # 设置配置文件为全路径
    args.reference= os.path.abspath(args.reference)

    # 检查文件是否存在
    if not os.path.exists(args.reference):
        logger.error(f"'{args.reference}' No such file or directory.")
        sys.exit(1)
    if not os.path.exists(args.config):
        logger.error(f"'{args.config}' No such file or directory.")
        sys.exit(1)

    # 检查参数是否正确
    if args.synRatio > 1 or args.synRatio < 0:
        logger.error(f"The synRatio value should be between 0 and 1: {args.ratio}.")
        sys.exit(1)
    if args.nosynRatio > 1 or args.nosynRatio < 0:
        logger.error(f"The nosynRatio value should be between 0 and 1: {args.ratio}.")
        sys.exit(1)

    # Set CWD and check if it exists
    if args.dir is None:
        args.dir = os.getcwd() + os.sep
    else:
        if os.path.isdir(args.dir):
            args.dir = os.path.abspath(args.dir + os.sep)
            # 切换工作路径
            os.chdir(os.path.abspath(args.dir))
        else:
            logger.error(args.dir + ' is not a valid folder. Exiting.')
            sys.exit(1)

    # set the dict of all genomes path
    config_file_map = {}  # dict{sample, {genome: "path", aligns: "path", syri_out: "path"}}
    with open(args.config, "r") as f:
        for info in f.readlines():
            # 跳过空行
            if len(info) == 0 or "#" in info:
                continue

            info = info.strip()
            infoList = info.split()

            # 检查列数是否符合规定
            if len(infoList) != 4:
                logger.error("Error: '" + args.config + f"' is not four columns -> {info}")
                sys.exit(1)

            # 记录该样品的基因组全路径
            genomeFilePath = os.path.abspath(infoList[1])
            alignsFilePath = os.path.abspath(infoList[2])
            syriOutFilePath = os.path.abspath(infoList[3])
            config_file_map[infoList[0]] = {
                "genome": genomeFilePath, 
                "aligns": alignsFilePath, 
                "syri_out": syriOutFilePath
            }
            # 检查文件是否存在且大于0
            for filePath in [genomeFilePath, alignsFilePath, syriOutFilePath]:
                if os.path.isfile(filePath) and os.path.getsize(filePath) > 0:
                    pass
                else:
                    logger.error(f"Error: '{filePath}': File does not exist or is empty.")
                    sys.exit(1)

    # 如果调试代码，进程池数量为1
    if args.debug:
        args.jobs = 1

    return args, config_file_map


# 运行 SynDiv_c multiinter
def multiinter(config_file_map, prefix, workDir):
    """
    :param config_file_map:   getParser输出每个样品所有需要样品的路径
    :param prefix:  输出前缀
    :param workDir:  工作路径
    :return: multiinterPath
    """
    # 创建工作路径
    makedir(workDir)
    workDir = workDir + os.sep
    # 切换工作路径
    os.chdir(workDir)

    # multiinter 命令
    cmd = f"{code_dir}SynDiv_c multiinter "
    filePathList = []
    sampleNameList = []

    for sampleName, value in config_file_map.items():
        filePathList.append(value["syri_out"])
        sampleNameList.append(sampleName)
    
    multiinterPath = os.path.join(workDir, f"{prefix}multiinter.out")
    cmd += "-i " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {multiinterPath}"

    logger.error(f'CMD: {cmd}')

    # 提交任务
    run_command(cmd, code_path)

    return multiinterPath


# 运行 SynDiv_c coor
def coor(config_file_map, multiinterPath, threads, prefix, workDir):
    """
    :param config_file_map:   getParser输出每个样品所有需要样品的路径
    :param multiinterPath:  SynDiv_c multiinter 输出的共线性相交结果
    :param threads:  线程数
    :param prefix:  输出前缀
    :param workDir:  工作路径
    :return: coorPath
    """
    # 创建工作路径
    makedir(workDir)
    workDir = workDir + os.sep
    # 切换工作路径
    os.chdir(workDir)

    # coor 命令
    cmd = f"{code_dir}SynDiv_c coor -s {multiinterPath} -t {threads} "
    filePathList = []
    sampleNameList = []

    for sampleName, value in config_file_map.items():
        filePathList.append(value["aligns"])
        sampleNameList.append(sampleName)
    
    coorPath = os.path.join(workDir, f"{prefix}coor.out")
    cmd += "-i " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {coorPath}"

    logger.error(f'CMD: {cmd}')

    # 提交任务
    run_command(cmd, code_path)

    return coorPath


# 运行 SynDiv_c no_syn
def no_syn(config_file_map, coorPath, prefix, workDir):
    """
    :param config_file_map:   getParser输出每个样品所有需要样品的路径
    :param coorPath:  SynDiv_c coor输出的共线性在qurey上的坐标
    :param prefix:  输出前缀
    :param workDir:  工作路径
    :return: coorPath
    """
    # 创建工作路径
    makedir(workDir)
    workDir = workDir + os.sep
    # 切换工作路径
    os.chdir(workDir)

    # 生成染色体长度信息
    filePathList = []
    sampleNameList = []
    for key, value in config_file_map.items():
        sampleNameList.append(key)
        filePath = os.path.join(workDir, f"{key}.length")
        filePathList.append(filePath)
        # 获取染色体长度文件
        getLength(value["genome"], filePath)

    # no_syn 命令
    cmd = f"{code_dir}SynDiv_c no_syn --coor {coorPath} "
    
    no_synPath = os.path.join(workDir, f"{prefix}coor.out")
    cmd += "--lengths " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {no_synPath}"

    logger.error(f'CMD: {cmd}')

    # 提交任务
    run_command(cmd, code_path)

    return no_synPath


# 运行 SynDiv_c cal
def cal(args, config_file_map, coorPath, no_syn_alignmentPath, workDir):
    """
    :param args:  getParser输出参数
    :param config_file_map:   getParser输出每个样品所有需要样品的路径
    :param coorPath:  SynDiv_c coor输出的共线性在qurey上的坐标
    :param no_syn_alignmentPath: no_syn_alignment 输出的配置文件路径
    :param workDir:  工作路径
    :return: calPath
    """
    # 创建工作路径
    makedir(workDir)
    workDir = workDir + os.sep
    # 切换工作路径
    os.chdir(workDir)

    # 参数信息
    filePathList = []
    sampleNameList = []
    for sampleName, value in config_file_map.items():
        filePathList.append(value["aligns"])
        sampleNameList.append(sampleName)

    # no_syn 命令
    cmd = f"{code_dir}SynDiv_c cal -t {args.threads} -r {args.reference} --coor {coorPath} --syri_outs {no_syn_alignmentPath} "

    if args.mode == "fast":
        cmd += "--fast "
    
    calPath = os.path.join(workDir, f"{args.prefix}coor.out")
    cmd += "--aligns " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {calPath}"

    logger.error(f'CMD: {cmd}')

    # 提交任务
    run_command(cmd, code_path)

    return calPath


# 运行 SynDiv_c cal
def window(args, calPath):
    """
    :param args:  getParser输出参数
    :param calPath:   SynDiv_c cal 计算的每个坐标的 syntenic diversity
    :return: windowPath
    """
    # 输出文件
    windowPath = os.path.abspath(f"{args.prefix}SynDiv.out")

    # no_syn 命令
    cmd = f"{code_dir}SynDiv_c window -r {args.reference} -i {calPath} -w {args.window} -s {args.step} -o {windowPath}"

    logger.error(f'CMD: {cmd}')

    # 提交任务
    run_command(cmd, code_path)

    return windowPath


def main():
    # ################################################# 解析参数 ################################################# #
    args, config_file_map = getParser()

    logger.error(f'Running.')

    # ################################################# SynDiv_c multiinter ################################################# #
    os.chdir(args.dir)
    multiinterPath = multiinter(config_file_map, args.prefix, os.path.join(args.dir, "multiinter"))

    # ################################################# SynDiv_c coor ################################################# #
    os.chdir(args.dir)
    coorPath = coor(config_file_map, multiinterPath, args.threads, args.prefix, os.path.join(args.dir, "coor"))

    # ################################################# SynDiv_c no_syn ################################################# #
    os.chdir(args.dir)
    no_synPath = no_syn(config_file_map, coorPath, args.prefix, os.path.join(args.dir, "no_syn"))

    # ################################################# no_syn_alignment ################################################# #
    os.chdir(args.dir)
    no_syn_alignmentPath = no_syn_alignment.main(args, config_file_map, no_synPath, os.path.join(args.dir, "no_syn_alignment"), code_path)

    # ################################################# no_syn_alignment ################################################# #
    os.chdir(args.dir)
    calPath = cal(args, config_file_map, coorPath, no_syn_alignmentPath, os.path.join(args.dir, "cal"))
    
    # ################################################# no_syn_alignment ################################################# #
    os.chdir(args.dir)
    winPath = window(args, calPath)

    logger.error(f'Result: {calPath}')
    logger.error(f'Result: {winPath}')
    logger.error(f'Done.')


if __name__ == '__main__':
    main()
