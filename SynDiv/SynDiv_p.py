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
import logging

# include
include_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
sys.path.append(os.path.join(include_dir, "include/"))
import no_syn_alignment

# ��������
code_path = {'PATH': os.environ.get('PATH')}


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
    required.add_argument("-i", dest="input", help="File containing non-syntenic coordinates, output by SynDiv no_syn", type=argparse.FileType('r'), required=True)
    required.add_argument("-g", dest="genomes", help="Genome configuration file (format: sample\tpath)", type=argparse.FileType('r'), required=True)

    other = parser.add_argument_group("Additional arguments")
    other.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="S", choices=['T', 'S', 'B', 'P'])
    other.add_argument('-f', dest='f', help='Filter out low quality and small alignments. Use this parameter to use the full list of alignments without any filtering.', default=True, action='store_false')
    other.add_argument('--dir', dest='dir', help="Path to working directory (if not current directory). All files must be in this directory.", action='store')
    other.add_argument("--prefix", dest="prefix", help="Prefix to add before the output file Names", type=str, default="SynDiv.")
    other.add_argument("--no-chrmatch", dest='chrmatch', help="Don't allow automatic matching chromosome ids between the two genomes if they are not equal", default=False, action='store_true')
    other.add_argument("--debug", dest='debug', help="debug code", default=False, action='store_true')

    # ���в���
    parallel_option = parser.add_argument_group(title="parallel")
    # ����
    parallel_option.add_argument("--jobs", dest="jobs", help="Run n jobs in parallel", type=int, default=100)
    # �߳�
    parallel_option.add_argument("--threads", dest="threads", help="Number of threads used per job", type=int, default=3)

    args = parser.parse_args()

    # ���������ļ�Ϊȫ·��
    noSynFilePath = os.path.abspath(args.input.name)

    # Set CWD and check if it exists
    if args.dir is None:
        args.dir = os.getcwd() + os.sep
    else:
        if os.path.isdir(args.dir):
            args.dir = args.dir + os.sep
            # �л�����·��
            os.chdir(args.dir)
        else:
            logger.error(args.dir + ' is not a valid folder. Exiting.')
            sys.exit(1)

    # set the dict of all genomes
    genomesFileMap = {}  # map<sample, genomeFileName>
    with open(args.genomes.name, "r") as f:
        for info in f.readlines():
            # ��������
            if len(info) == 0 or "#" in info:
                continue

            info = info.strip()
            infoList = info.split()

            # ��������Ƿ���Ϲ涨
            if len(infoList) != 2:
                logger.error("Error: '" + args.genomes.name + f"' is not two columns -> {info}")
                sys.exit(1)

            # ��¼����Ʒ�Ļ�����ȫ·��
            genomeFilePath = os.path.abspath(infoList[1])
            genomesFileMap[infoList[0]] = genomeFilePath
            # ����ļ��Ƿ�����Ҵ���0
            if os.path.isfile(genomeFilePath) and os.path.getsize(genomeFilePath) > 0:
                pass
            else:
                logger.error(f"Error: '{genomeFilePath}': File does not exist or is empty.")
                sys.exit(1)

    # ������Դ��룬���̳�����Ϊ1
    if args.debug:
        args.jobs = 1

    return args, genomesFileMap, noSynFilePath


def main():
    logger = logging.getLogger('SynDiv')

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Running.')

    # ��������
    args, genomesFileMap, noSynFilePath = getParser()

    # ################################################# no_syn_alignment ################################################# #
    os.chdir(args.dir)
    no_syn_alignmentPath = no_syn_alignment.main(args, genomesFileMap, noSynFilePath, os.path.join(args.dir, "no_syn_alignment"), code_path)

    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Result: {no_syn_alignmentPath}')
    logger.error(f'[{str(datetime.datetime.now()).split(".")[0]}] Done.')


if __name__ == '__main__':
    main()
