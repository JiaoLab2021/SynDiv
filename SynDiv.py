#!/usr/bin/env python3

# -*- coding: utf-8 -*-

__data__ = "2024/02/23"
__version__ = "1.1.1"
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
import no_syn_alignment, run_cmd, fasta_length
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
        self.required.add_argument("-r, --reference", dest="reference", help="Input FASTA reference", type=str, required=True)
        self.required.add_argument("-c, --config", dest="config", help="Configuration file (format: sample\tgenome_path\taligns_path\tsyri_out_path)", type=str, required=True)
        # Alignment options
        self.AlignmentFun = self.parser.add_argument_group("Alignment arguments")
        self.AlignmentFun.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="S", choices=['T', 'S', 'B', 'P'])
        self.AlignmentFun.add_argument('-f', dest='f', help='Filter out low quality and small alignments. Use this parameter to use the full list of alignments without any filtering.', default=True, action='store_false')
        self.AlignmentFun.add_argument("--no-chrmatch", dest='chrmatch', help="Don't allow automatic matching chromosome ids between the two genomes if they are not equal", default=False, action='store_true')
        self.AlignmentFun.add_argument('--synRatio', dest="synRatio", help="Threshold for complete synteny detection. Lower values increase alignment speed (0,1].", type=float, default=0.8)
        self.AlignmentFun.add_argument('--nosynRatio', dest="nosynRatio", help="Threshold for complete no-synteny detection. Higher values increase alignment speed (0,1].", type=float, default=0.05)
        # SynDic_c cal options
        self.calFun = self.parser.add_argument_group("SynDic_c cal arguments")
        self.calFun.add_argument('--mode', dest="mode", help="Enabling quick mode will increase memory consumption", default="normal", choices=["fast", "normal"])
        self.calFun.add_argument('--buffer', dest="buffer", help="Buffer size for file reading, measured in MB [1].", type=int, default=1)
        # SynDic_c window options
        self.windowFun = self.parser.add_argument_group("SynDic_c window arguments")
        self.windowFun.add_argument('--window', dest="window", help="Window size", type=int, default=5000)
        self.windowFun.add_argument('--step', dest="step", help="Step size", type=int, default=1000)
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
        self.reference = os.path.abspath(self.args.reference)
        self.config = os.path.abspath(self.args.config)
        # Check if the file exists
        if not os.path.exists(self.reference):
            self.logger.error(f"'{self.reference}' No such file or directory.")
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
        self.configFileMap = {}  # dict{sample, {genome: "path", aligns: "path", syri: "path"}}

        # Multiinter output
        self.multiinterFilePath = ""

        # Coor output
        self.coorFilePath = ""

        # No_syn output
        self.no_synFilePath = ""

        # No_syn_alignment output
        self.no_syn_alignmentFilePath = ""

        # Cal output
        self.calFilePath = ""

        # Window output
        self.winFilePath = ""

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
            if len(infoList) != 4:
                self.logger.error(f"Error: The configuration file '{self.args.config}' is expected to have four columns, but the following line was found -> {info}")
                sys.exit(1)

            # Record the full genome path of the sample
            genomeFilePath = os.path.abspath(infoList[1])
            alignsFilePath = os.path.abspath(infoList[2])
            syriOutFilePath = os.path.abspath(infoList[3])
            self.configFileMap[infoList[0]] = {
                "genome": genomeFilePath, 
                "aligns": alignsFilePath, 
                "syri": syriOutFilePath
            }

            # Check if the file exists and is greater than 0
            for filePath in [genomeFilePath, alignsFilePath, syriOutFilePath]:
                if not (os.path.isfile(filePath) and os.path.getsize(filePath) > 0):
                    self.logger.error(f"Error: The file '{filePath}' either does not exist or is empty.")
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

    # Run SynDiv_c multiinter
    def multiinter(self):
        # mkdir
        multiinterWorkDir = os.path.join(self.workDir, "multiinter")
        self.makedir(multiinterWorkDir)
        os.chdir(multiinterWorkDir)

        # multiinter command
        cmd = os.path.join(code_dir, "SynDiv_c multiinter ")

        filePathList = [value["syri"] for sampleName, value in self.configFileMap.items()]
        sampleNameList = list(self.configFileMap.keys())
        
        self.multiinterFilePath = os.path.join(multiinterWorkDir, f"{self.args.prefix}.multiinter.out")
        cmd += "-i " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {self.multiinterFilePath}"

        self.logger.error(f'CMD: {cmd}')

        # Run the command
        stdout, stderr, log = run_cmd.run(cmd, code_env)
        if log:
            raise SystemExit(log)
    
    # Run SynDiv_c coor
    def coor(self):
        # mkdir
        coorWorkDir = os.path.join(self.workDir, "coor")
        self.makedir(coorWorkDir)
        os.chdir(coorWorkDir)

        # coor command
        cmd = f"{os.path.join(code_dir, 'SynDiv_c coor')} -s {self.multiinterFilePath} -t {self.args.threads} "

        filePathList = [value["aligns"] for sampleName, value in self.configFileMap.items()]
        sampleNameList = list(self.configFileMap.keys())
        
        self.coorFilePath = os.path.join(coorWorkDir, f"{self.args.prefix}.coor.out")
        cmd += "-i " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {self.coorFilePath}"

        self.logger.error(f'CMD: {cmd}')

        # Run the command
        stdout, stderr, log = run_cmd.run(cmd, code_env)
        if log:
            raise SystemExit(log)

    # Run SynDiv_c no_syn
    def no_syn(self):
        # mkdir
        no_synWorkDir = os.path.join(self.workDir, "no_syn")
        self.makedir(no_synWorkDir)
        os.chdir(no_synWorkDir)

        filePathList = []
        sampleNameList = list(self.configFileMap.keys())

        for key, value in self.configFileMap.items():
            filePath = os.path.join(no_synWorkDir, f"{key}.length")
            filePathList.append(filePath)
            # Get Chromosome Length File
            fasta_length.getLength(value["genome"], filePath)

        # no_syn command
        cmd = f"{os.path.join(code_dir, 'SynDiv_c no_syn')} --coor {self.coorFilePath} "
        
        self.no_synFilePath = os.path.join(no_synWorkDir, f"{self.args.prefix}.no_syn.out")
        cmd += "--lengths " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {self.no_synFilePath}"

        self.logger.error(f'CMD: {cmd}')

        # Run the command
        stdout, stderr, log = run_cmd.run(cmd, code_env)
        if log:
            raise SystemExit(log)

    # no_syn_alignment
    def no_syn_alignment(self):
        # mkdir
        no_syn_alignmentWorkDir = os.path.join(self.workDir, "no_syn_alignment")
        self.no_syn_alignmentFilePath = no_syn_alignment.main(self.args, self.configFileMap, self.no_synFilePath, no_syn_alignmentWorkDir, code_env)

    # Run SynDiv_c cal
    def cal(self):
        # mkdir
        calWorkDir = os.path.join(self.workDir, "cal")
        self.makedir(calWorkDir)
        os.chdir(calWorkDir)

        filePathList = [value["aligns"] for sampleName, value in self.configFileMap.items()]
        sampleNameList = list(self.configFileMap.keys())

        # no_syn command
        cmd = f"{os.path.join(code_dir, 'SynDiv_c cal')} -t {self.args.threads} -r {self.reference} --coor {self.coorFilePath} --syri_outs {self.no_syn_alignmentFilePath} --buffer {self.args.buffer} "

        if self.args.mode == "fast":
            cmd += "--fast "

        self.calFilePath = os.path.join(calWorkDir, f"{self.args.prefix}.cal.out")
        cmd += "--aligns " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {self.calFilePath}"

        self.logger.error(f'CMD: {cmd}')

        # Run the command
        stdout, stderr, log = run_cmd.run(cmd, code_env)
        if log:
            raise SystemExit(log)
    
    # Run SynDiv_c window
    def window(self):
        # chdir
        os.chdir(self.workDir)

        # output file
        self.winFilePath = os.path.join(self.workDir, f"{self.args.prefix}.win.out")

        # no_syn command
        cmd = f"{os.path.join(code_dir, 'SynDiv_c window')} -r {self.reference} -i {self.calFilePath} -w {self.args.window} -s {self.args.step} -o {self.winFilePath}"

        self.logger.error(f'CMD: {cmd}')

        # Run the command
        stdout, stderr, log = run_cmd.run(cmd, code_env)
        if log:
            raise SystemExit(log)

    # Print Result
    def print_result(self):
        self.logger.error(f'Result: {self.calFilePath}')
        self.logger.error(f'Result: {self.winFilePath}')
        self.logger.error(f'All Done ...')


def main():
    # SynDiv
    SynDivClass = MySynDiv()
    # set_file_path_dict
    SynDivClass.set_file_path_dict()
    # multiinter
    SynDivClass.multiinter()
    # coor
    SynDivClass.coor()
    # no_syn
    SynDivClass.no_syn()
    # no_syn_alignment
    SynDivClass.no_syn_alignment()
    # cal
    SynDivClass.cal()
    # window
    SynDivClass.window()
    # print_result
    SynDivClass.print_result()


if __name__ == '__main__':
    main()