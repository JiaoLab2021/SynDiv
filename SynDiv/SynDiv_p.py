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
import logging

# include
include_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
sys.path.append(os.path.join(include_dir, "include/"))
import no_syn_alignment

# 环境变量
code_path = {'PATH': os.environ.get('PATH')}

# log
logger = logging.getLogger('SynDiv_p')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # 输出到控制台
handler.setFormatter(formatter)
logger.addHandler(handler)


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
    required.add_argument("-i, --input", dest="input", help="File containing non-syntenic coordinates, output by SynDiv no_syn", type=str, required=True)
    required.add_argument("-c, --config", dest="config", help="Genome configuration file (format: sample\tpath)", type=str, required=True)

    other = parser.add_argument_group("Additional arguments")
    other.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="S", choices=['T', 'S', 'B', 'P'])
    other.add_argument('-f', dest='f', help='Filter out low quality and small alignments. Use this parameter to use the full list of alignments without any filtering.', default=True, action='store_false')
    other.add_argument('--dir', dest='dir', help="Path to working directory (if not current directory). All files must be in this directory.", action='store')
    other.add_argument("--prefix", dest="prefix", help="Prefix to add before the output file Names", type=str, default="SynDiv.")
    other.add_argument("--no-chrmatch", dest='chrmatch', help="Don't allow automatic matching chromosome ids between the two genomes if they are not equal", default=False, action='store_true')
    other.add_argument('--synRatio', dest="synRatio", help="Threshold for complete synteny detection. Lower values increase alignment speed (0,1].", type=float, default=0.8)
    other.add_argument('--nosynRatio', dest="nosynRatio", help="Threshold for complete no-synteny detection. Higher values increase alignment speed (0,1].", type=float, default=0.05)
    other.add_argument("--debug", dest='debug', help="Debug code", default=False, action='store_true')

    # 并行参数
    parallel_option = parser.add_argument_group(title="parallel")
    # 进程
    parallel_option.add_argument("--jobs", dest="jobs", help="Run n jobs in parallel", type=int, default=150)
    # 线程
    parallel_option.add_argument("--threads", dest="threads", help="Number of threads used by SynDiv_c", type=int, default=10)

    args = parser.parse_args()

    # 设置配置文件为全路径
    args.input = os.path.abspath(args.input)

    # 检查文件是否存在
    if not os.path.exists(args.input):
        logger.error(f"'{args.input}' No such file or directory.")
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
            args.dir = args.dir + os.sep
            # 切换工作路径
            os.chdir(args.dir)
        else:
            logger.error(args.dir + ' is not a valid folder. Exiting.')
            sys.exit(1)

    # set the dict of all genomes
    genomesFileMap = {}  # map<sample, genomeFileName>
    with open(args.config, "r") as f:
        for info in f.readlines():
            # 跳过空行
            if len(info) == 0 or "#" in info:
                continue

            info = info.strip()
            infoList = info.split()

            # 检查列数是否符合规定
            if len(infoList) != 2:
                logger.error("Error: '" + args.config + f"' is not two columns -> {info}")
                sys.exit(1)

            # 记录该样品的基因组全路径
            genomeFilePath = os.path.abspath(infoList[1])
            genomesFileMap[infoList[0]] = genomeFilePath
            # 检查文件是否存在且大于0
            if os.path.isfile(genomeFilePath) and os.path.getsize(genomeFilePath) > 0:
                pass
            else:
                logger.error(f"Error: '{genomeFilePath}': File does not exist or is empty.")
                sys.exit(1)

    # 如果调试代码，进程池数量为1
    if args.debug:
        args.jobs = 1

    return args, genomesFileMap


def main():
    # 解析参数
    args, genomesFileMap = getParser()

    logger.error(f'Running.')

    # ################################################# no_syn_alignment ################################################# #
    os.chdir(args.dir)
    no_syn_alignmentPath = no_syn_alignment.main(args, genomesFileMap, args.input, os.path.join(args.dir, "no_syn_alignment"), code_path)

    logger.error(f'Result: {no_syn_alignmentPath}')
    logger.error(f'Done.')


if __name__ == '__main__':
    main()
