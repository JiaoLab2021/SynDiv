#!/usr/bin/env python3

# -*- coding: utf-8 -*-

__data__ = "2024/06/12"
__version__ = "1.1.2"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
import argparse
import sys
import shutil
import logging
import multiprocessing

# include
code_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(code_dir, "include/"))
from file_open import MyOpener
import run_cmd

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
        self.required.add_argument("-r, --refgenome", dest="refgenome", help="Input FASTA reference", type=str, required=True)
        self.required.add_argument("-c, --config", dest="config", help="Genome configuration file (format: sample\tpath)", type=str, required=True)
        self.required.add_argument("-o, --out", dest="out", help="Output file name", type=str, required=False, default="SynDiv.config")
        # Parallel options
        self.parallel_option = self.parser.add_argument_group(title="parallel")
        self.parallel_option.add_argument("--jobs", dest="jobs", help="Run n jobs in parallel", type=int, default=3)
        self.parallel_option.add_argument("--threads", dest="threads", help="Number of threads used by nucmer", type=int, default=10)

    def setup_paths(self):
        """Setup paths based on parsed arguments"""
        # input
        self.refgenome = os.path.abspath(self.args.refgenome)
        self.config = os.path.abspath(self.args.config)
        self.output = os.path.abspath(self.args.out)
        # Check if the file exists
        if not os.path.exists(self.refgenome):
            self.logger.error(f"'{self.refgenome}' No such file or directory.")
            sys.exit(1)
        if not os.path.exists(self.config):
            self.logger.error(f"'{self.config}' No such file or directory.")
            sys.exit(1)

        # work dir
        self.workDir = os.getcwd()


# genome2SynDiv
class MyGenome2SynDiv(MyParser):
    def __init__(self):
        super().__init__()

        # Create a process pool, specify the maximum number of processes as self.args.jobs
        self.pool = multiprocessing.Pool(processes=self.args.jobs)

        # chromosome list
        self.chromosomes = []

        # output dict
        self.outputFileMap = {}  # dict{sample: {genome: genomeFileName, aligns: alignsFileName, syri: syriFileName}}

    # chromosome list
    def get_chromosome_list(self):
        with open(self.refgenome, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    # Assuming the chromosome number follows '>'
                    chromosome = line[1:].strip().split()[0]
                    self.chromosomes.append(chromosome)

    # query genomes
    def set_genome_path(self):
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
            self.outputFileMap[infoList[0]] = {
                "genome": genomeFilePath, 
                "aligns": "", 
                "syri": ""}
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

    # alignment
    def alignment_syri(self):
        # mkdir
        workDir = os.path.join(self.workDir, "nucmer_syri")
        self.makedir(workDir)

        results = []
        for key, value in self.outputFileMap.items():
            sampleName = key
            sampleGenomePath = value["genome"]

            # sbumit task
            result = self.pool.apply_async(
                nucmer_syri, 
                args=(
                    self.refgenome, 
                    sampleName,
                    sampleGenomePath, 
                    self.chromosomes, 
                    workDir, 
                    self.args.threads, 
                    code_env
                )
            )
            results.append(result)

        # Get the return value of the function
        for result in results:
            try:
                log, sampleName, sampleAligns, sampleSyri = result.get()
                if log:
                    self.logger.error(log)
                    self.logger.error(f"Error: Processing of '{sampleName}' failed.")
                else:
                    self.logger.error(f"Processing of '{sampleName}' completed successfully.")
                self.outputFileMap[sampleName]["aligns"] = sampleAligns
                self.outputFileMap[sampleName]["syri"] = sampleSyri
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")
    
    # save result
    def save_result(self):
        with open(f"{self.output}", "w") as f:
            lines = [f"{key}\t{value['genome']}\t{value['aligns']}\t{value['syri']}\n" for key, value in self.outputFileMap.items()]
            f.writelines(lines)

    # Print Result
    def print_result(self):
        self.logger.error(f'Result: {os.path.abspath(self.output)}')
        self.logger.error(f'All Done ...')


# nucmer
def nucmer_syri(refgenomePath, sampleName, sampleGenomePath, chromosomes, workDir, thread, envPath):
    os.chdir(workDir)

    commands = [
        f"nucmer -l 40 -g 90 -c 100 -b 200 -t {thread} -p {sampleName} {refgenomePath} {sampleGenomePath}",
        f"delta-filter -1 -i 90 -l 100 {sampleName}.delta 1>{sampleName}.filter.delta",
        f"show-coords -THrd {sampleName}.filter.delta 1>{sampleName}.coords"
    ]

    for chromosome in chromosomes:
        commands.append(f"show-aligns -r {sampleName}.filter.delta {chromosome} {chromosome} 1>{sampleName}.{chromosome}.aligns")

    aligns_paths = [f"{sampleName}.{chromosome}.aligns" for chromosome in chromosomes]
    commands.append(f"cat {' '.join(aligns_paths)} > {sampleName}.aligns")

    commands.append(f"syri -c {sampleName}.coords -d {sampleName}.filter.delta -r {refgenomePath} -q {sampleGenomePath} --nc 5 --all --prefix {sampleName}.")

    for cmd in commands:
        stdout, stderr, log = run_cmd.run(cmd, envPath)
        if log:
            return log, sampleName, f"{sampleName}.aligns", f"{sampleName}.syri.out"

    return log, sampleName, os.path.abspath(f"{sampleName}.aligns"), os.path.abspath(f"{sampleName}.syri.out")


def main():
    MyGenome2SynDivClass = MyGenome2SynDiv()
    MyGenome2SynDivClass.get_chromosome_list()
    MyGenome2SynDivClass.set_genome_path()
    MyGenome2SynDivClass.alignment_syri()
    MyGenome2SynDivClass.print_result()
    MyGenome2SynDivClass.save_result()


if __name__ == '__main__':
    main()
