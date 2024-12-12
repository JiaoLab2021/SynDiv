#!/usr/bin/env python3

# -*- coding: utf-8 -*-

__data__ = "2024/12/11"
__version__ = "1.1.3"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
import argparse
import sys
import shutil
import logging

# include
code_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(code_dir, "include/"))
import no_syn_alignment
from file_open import MyOpener

# environment variable
code_env = {'PATH': os.environ.get('PATH')}


# get parser
class MyParser:
    def __init__(self):
        self.setup_logging()
        self.initialize_parser()

        # Parse arguments
        self.args = self.parser.parse_args()
        self.setup_paths()

    def setup_logging(self):
        """Configure the logging"""
        self.logger = logging.getLogger('MyParser')
        # log
        self.logger.error(f"data: {__data__}")
        self.logger.error(f"version: {__version__}")
        self.logger.error(f"author: {__author__}")
        self.logger.error(f"If you encounter any issues related to the code, please don't hesitate to contact us via email at {__email__}.\n")
        formatter = logging.Formatter('[%(asctime)s] %(message)s')
        handler = logging.StreamHandler()  # output to the console
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def initialize_parser(self):
        """Setup command line argument parser"""
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        self.optional = self.parser._action_groups.pop()
        # Input options
        self.required = self.parser.add_argument_group("Input Files")
        self.required.add_argument("-i, --input", dest="input", help="File containing non-syntenic coordinates, output by 'SynDiv no_syn'", type=str, required=True)
        self.required.add_argument("-c, --config", dest="config", help="Genome configuration file (format: sample\tpath)", type=str, required=True)
        # Alignment options
        self.AlignmentFun = self.parser.add_argument_group("Alignment arguments")
        self.AlignmentFun.add_argument('-F', dest="ftype", help=argparse.SUPPRESS, default="S", choices=['T', 'S', 'B', 'P'])
        self.AlignmentFun.add_argument('-f', dest='f', help=argparse.SUPPRESS, default=True, action='store_false')
        self.AlignmentFun.add_argument("--no-chrmatch", dest='chrmatch', help=argparse.SUPPRESS, default=False, action='store_true')
        self.AlignmentFun.add_argument('--synRatio', dest="synRatio", help="Threshold for complete synteny detection. Lower values increase alignment speed (0,1].", type=float, default=0.8)
        self.AlignmentFun.add_argument('--nosynRatio', dest="nosynRatio", help="Threshold for complete no-synteny detection. Higher values increase alignment speed (0,1].", type=float, default=0.05)
        # Parallel options
        self.parallel_option = self.parser.add_argument_group(title="parallel")
        self.parallel_option.add_argument("--jobs", dest="jobs", help="Run n jobs in parallel", type=int, default=30)
        self.parallel_option.add_argument("--threads", dest="threads", help="Number of threads used by SynDiv_c", type=int, default=10)
        # Other options
        self.other = self.parser.add_argument_group("Additional arguments")
        self.other.add_argument('--dir', dest='dir', help="Path to working directory (if not current directory)", action='store')
        self.other.add_argument("--prefix", dest="prefix", help="Prefix to add before the output file Names", type=str, default="SynDiv")
        self.other.add_argument("--debug", dest='debug', help="Debug code", default=False, action='store_true')

    def setup_paths(self):
        """Setup paths based on parsed arguments"""
        # input
        self.input = os.path.abspath(self.args.input)
        self.config = os.path.abspath(self.args.config)
        # Check if the file exists
        if not os.path.exists(self.input):
            self.logger.error(f"'{self.input}' No such file or directory.")
            sys.exit(1)
        if not os.path.exists(self.config):
            self.logger.error(f"'{self.config}' No such file or directory.")
            sys.exit(1)

        # Check that the parameters are correct
        if self.args.synRatio > 1 or self.args.synRatio < 0:
            self.logger.error(f"The 'synRatio' value must be between 0 and 1. Current value: {self.args.ratio}.")
            sys.exit(1)
        if self.args.nosynRatio > 1 or self.args.nosynRatio < 0:
            self.logger.error(f"The 'nosynRatio' value must be between 0 and 1. Current value: {self.args.ratio}.")
            sys.exit(1)

        # path
        self.workDir = os.getcwd() if self.args.dir is None else os.path.abspath(self.args.dir)
        # Check if the directory exists, if not create it
        if not os.path.isdir(self.workDir):
            os.makedirs(self.workDir)

        # debug
        if self.args.debug:
            self.args.jobs = 1
            self.args.threads = 1

# SynDiv
class MySynDiv(MyParser):
    def __init__(self):
        super().__init__()

        # config file
        self.genomesFileMap = {}  # dict{sample: genomeFileName}

        # No_syn_alignment output
        self.no_syn_alignmentFilePath = ""

    def set_file_path_dict(self):
        self.logger.error(f"Running ...")

        fileOpener = MyOpener(self.config)

        line = [None]
        while fileOpener.read_line(line):
            info = line[0].strip()

            # skip empty lines and comments
            if not info.strip() or info.strip().startswith("#"):
                continue

            infoList = info.split()

            # Check that the number of columns meets the specification
            if len(infoList) != 2:
                self.logger.error("Error: '" + self.args.config + f"' is not two columns -> {info}")
                sys.exit(1)

            # Record the full genome path of the sample
            genomeFilePath = os.path.abspath(infoList[1])
            self.genomesFileMap[infoList[0]] = genomeFilePath
            # Check if the file exists and is greater than 0
            if not os.path.isfile(genomeFilePath) or os.path.getsize(genomeFilePath) == 0:
                self.logger.error(f"Error: The file '{genomeFilePath}' does not exist or is empty.")
                sys.exit(1)

    # Create a directory
    def makedir(self, path_dir: str):
        if os.path.isdir(path_dir):
            shutil.rmtree(path_dir)
            os.makedirs(path_dir)
            log = f'\'{path_dir}\' already exists, clear and recreate.'
            self.logger.error(log)
        else:
            os.makedirs(path_dir)

    # no_syn_alignment
    def no_syn_alignment(self):
        # mkdir
        no_syn_alignmentWorkDir = os.path.join(self.workDir, "no_syn_alignment")
        self.no_syn_alignmentFilePath = no_syn_alignment.main(self.args, self.genomesFileMap, self.input, no_syn_alignmentWorkDir, code_env)

    # Print Result
    def print_result(self):
        self.logger.error(f'Result: {self.no_syn_alignmentFilePath}')
        self.logger.error(f'All Done ...')


def main():
    # SynDiv
    SynDivClass = MySynDiv()
    # set_file_path_dict
    SynDivClass.set_file_path_dict()
    # no_syn_alignment
    SynDivClass.no_syn_alignment()
    # print_result
    SynDivClass.print_result()


if __name__ == '__main__':
    main()
