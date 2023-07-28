#!/usr/bin/env python3

# -*- coding: utf-8 -*-

__data__ = "2023/07/28"
__version__ = "1.0.5"
__author__ = "Zezhen Du"
__email__ = "dzz0539@gmail.com or dzz0539@163.com"

import os
import argparse
import sys
import shutil
import subprocess
import tempfile
import logging

# include
include_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
sys.path.append(os.path.join(include_dir, "include/"))
import no_syn_alignment

code_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep

# environment variable
code_path = {'PATH': os.environ.get('PATH')}

# log
logger = logging.getLogger('SynDiv')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # output to the console
handler.setFormatter(formatter)
logger.addHandler(handler)


# Create a directory
def makedir(
    path_dir: str
):
    """
    :param path_dir: The folder path to be created
    :return: 0
    """
    if os.path.isdir(path_dir):
        shutil.rmtree(path_dir)
        os.makedirs(path_dir)
        log = '[makedir] \'{}\' already exists, clear and recreate.'.format(path_dir)
        logger.error(log)
    else:
        os.makedirs(path_dir)


# Run multithreaded tasks
def run_command(command, envPath):
    # Create output file
    stdout_file = tempfile.NamedTemporaryFile(delete=False)
    stderr_file = tempfile.NamedTemporaryFile(delete=False)

    # submit task
    proc = subprocess.Popen(command, shell=True, stdout=stdout_file, stderr=stderr_file, env=envPath)

    # Wait for the command to complete
    proc.wait()

    # Read the output file and close the file handle
    with open(stdout_file.name, 'rb') as f:
        stdout_data = f.read()
    os.unlink(stdout_file.name)

    with open(stderr_file.name, 'rb') as f:
        stderr_data = f.read()
    os.unlink(stderr_file.name)

    exit_state = proc.returncode
    stdout = stdout_data.strip().decode('utf-8')
    stderr = stderr_data.strip().decode('utf-8')

    # standard output and error output
    if exit_state != 0 or "FileNotFoundError" in stderr or \
            "command not found" in stderr or \
            "error" in stderr or \
            "Error" in stderr:
        logger.error(f"Error: Exit status is {exit_state}")
        logger.error(f"Error: CMD -> {command}")
        logger.error(f"Error: stderr -> {stderr}")
        sys.exit(1)

    return exit_state, stdout, stderr


# Get the length of each sequence in fasta
def getLength(fastaFilePath, outputFilePath):
    """
    :param fastaFilePath:   The fasta file path that needs to be calculated
    :param outputFilePath:  output length to file
    :return: 0
    """
    with open(fastaFilePath, "r") as f:
        seq = ""
        name = ""
        length_dict = {}
        for line in f:
            if line.startswith('>'):
                if seq != "":
                    # Count the length of the sequence and save it in the dictionary
                    length_dict[name] = len(seq)
                    seq = ""
                name = line.strip()[1:]
            else:
                seq += line.strip()
        # process the last sequence
        length_dict[name] = len(seq)

    # Output the ID and length of each sequence
    with open(outputFilePath, "w") as f:
        for seq_id, seq_len in length_dict.items():
            f.write(f"{seq_id}\t{seq_len}\n")

    return 0


# parsing parameters
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

    AlignmentFun = parser.add_argument_group("Alignment arguments")
    AlignmentFun.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="S", choices=['T', 'S', 'B', 'P'])
    AlignmentFun.add_argument('-f', dest='f', help='Filter out low quality and small alignments. Use this parameter to use the full list of alignments without any filtering.', default=True, action='store_false')
    AlignmentFun.add_argument("--no-chrmatch", dest='chrmatch', help="Don't allow automatic matching chromosome ids between the two genomes if they are not equal", default=False, action='store_true')
    AlignmentFun.add_argument('--synRatio', dest="synRatio", help="Threshold for complete synteny detection. Lower values increase alignment speed (0,1].", type=float, default=0.8)
    AlignmentFun.add_argument('--nosynRatio', dest="nosynRatio", help="Threshold for complete no-synteny detection. Higher values increase alignment speed (0,1].", type=float, default=0.05)

    calFun = parser.add_argument_group("SynDic_c cal arguments")
    calFun.add_argument('--mode', dest="mode", help="Enabling quick mode will increase memory consumption", default="normal", choices=["fast", "normal"])
    calFun.add_argument('--buffer', dest="buffer", help="Buffer size for file reading, measured in MB [1].", type=int, default=1)

    windowFun = parser.add_argument_group("SynDic_c window arguments")
    windowFun.add_argument('--window', dest="window", help="Window size", type=int, default=5000)
    windowFun.add_argument('--step', dest="step", help="Step size", type=int, default=1000)
    
    # parallel parameter
    parallel_option = parser.add_argument_group(title="parallel")
    # jobs
    parallel_option.add_argument("--jobs", dest="jobs", help="Run n jobs in parallel", type=int, default=30)
    # threads
    parallel_option.add_argument("--threads", dest="threads", help="Number of threads used by SynDiv_c", type=int, default=10)

    # other
    other = parser.add_argument_group("Additional arguments")
    other.add_argument('--dir', dest='dir', help="Path to working directory (if not current directory). All files must be in this directory.", action='store')
    other.add_argument("--prefix", dest="prefix", help="Prefix to add before the output file Names", type=str, default="SynDiv.")
    other.add_argument("--debug", dest='debug', help="Debug code", default=False, action='store_true')
    
    args = parser.parse_args()

    # Set the configuration file to the full path
    args.reference= os.path.abspath(args.reference)
    args.config = os.path.abspath(args.config)

    # Check if the file exists
    if not os.path.exists(args.reference):
        logger.error(f"'{args.reference}' No such file or directory.")
        sys.exit(1)
    if not os.path.exists(args.config):
        logger.error(f"'{args.config}' No such file or directory.")
        sys.exit(1)

    # Check that the parameters are correct
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
            # cd
            os.chdir(os.path.abspath(args.dir))
        else:
            logger.error(args.dir + ' is not a valid folder. Exiting.')
            sys.exit(1)

    # set the dict of all genomes path
    config_file_map = {}  # dict{sample, {genome: "path", aligns: "path", syri_out: "path"}}
    with open(args.config, "r") as f:
        for info in f.readlines():
            # skip empty lines
            if len(info) == 0 or "#" in info:
                continue

            info = info.strip()
            infoList = info.split()

            # Check that the number of columns meets the specification
            if len(infoList) != 4:
                logger.error("Error: '" + args.config + f"' is not four columns -> {info}")
                sys.exit(1)

            # Record the full genome path of the sample
            genomeFilePath = os.path.abspath(infoList[1])
            alignsFilePath = os.path.abspath(infoList[2])
            syriOutFilePath = os.path.abspath(infoList[3])
            config_file_map[infoList[0]] = {
                "genome": genomeFilePath, 
                "aligns": alignsFilePath, 
                "syri_out": syriOutFilePath
            }
            # Check if the file exists and is greater than 0
            for filePath in [genomeFilePath, alignsFilePath, syriOutFilePath]:
                if os.path.isfile(filePath) and os.path.getsize(filePath) > 0:
                    pass
                else:
                    logger.error(f"Error: '{filePath}': File does not exist or is empty.")
                    sys.exit(1)

    # If you debug the code, the number of process pools is 1
    if args.debug:
        args.jobs = 1

    return args, config_file_map


# Run SynDiv_c multiinter
def multiinter(config_file_map, prefix, workDir):
    """
    :param config_file_map:   getParser outputs paths to all required samples for each sample
    :param prefix:            output prefix
    :param workDir:           work path
    :return:                  multiinterPath
    """
    # mkdir
    makedir(workDir)
    workDir = workDir + os.sep
    # cd
    os.chdir(workDir)

    # multiinter command
    cmd = f"{code_dir}SynDiv_c multiinter "
    filePathList = []
    sampleNameList = []

    for sampleName, value in config_file_map.items():
        filePathList.append(value["syri_out"])
        sampleNameList.append(sampleName)
    
    multiinterPath = os.path.join(workDir, f"{prefix}multiinter.out")
    cmd += "-i " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {multiinterPath}"

    logger.error(f'CMD: {cmd}')

    # submit task
    run_command(cmd, code_path)

    return multiinterPath


# Run SynDiv_c coor
def coor(config_file_map, multiinterPath, threads, prefix, workDir):
    """
    :param config_file_map:   getParser outputs paths to all required samples for each sample
    :param multiinterPath:  Collinear intersection result of SynDiv_c multiinter output
    :param threads:  threads
    :param prefix:  output prefix
    :param workDir:  work path
    :return: coorPath
    """
    # mkdir
    makedir(workDir)
    workDir = workDir + os.sep
    # cd
    os.chdir(workDir)

    # coor command
    cmd = f"{code_dir}SynDiv_c coor -s {multiinterPath} -t {threads} "
    filePathList = []
    sampleNameList = []

    for sampleName, value in config_file_map.items():
        filePathList.append(value["aligns"])
        sampleNameList.append(sampleName)
    
    coorPath = os.path.join(workDir, f"{prefix}coor.out")
    cmd += "-i " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {coorPath}"

    logger.error(f'CMD: {cmd}')

    # submit task
    run_command(cmd, code_path)

    return coorPath


# Run SynDiv_c no_syn
def no_syn(config_file_map, coorPath, prefix, workDir):
    """
    :param config_file_map:   getParser outputs paths to all required samples for each sample
    :param coorPath:  The coordinates of the collinearity output by SynDiv_c coor on qurey
    :param prefix:  output prefix
    :param workDir:  work path
    :return: coorPath
    """
    # mkdir
    makedir(workDir)
    workDir = workDir + os.sep
    # cd
    os.chdir(workDir)

    # Generating Chromosome Length Information
    filePathList = []
    sampleNameList = []
    for key, value in config_file_map.items():
        sampleNameList.append(key)
        filePath = os.path.join(workDir, f"{key}.length")
        filePathList.append(filePath)
        # Get Chromosome Length File
        getLength(value["genome"], filePath)

    # no_syn command
    cmd = f"{code_dir}SynDiv_c no_syn --coor {coorPath} "
    
    no_synPath = os.path.join(workDir, f"{prefix}no_syn.out")
    cmd += "--lengths " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {no_synPath}"

    logger.error(f'CMD: {cmd}')

    # submit task
    run_command(cmd, code_path)

    return no_synPath


# Run SynDiv_c cal
def cal(args, config_file_map, coorPath, no_syn_alignmentPath, workDir):
    """
    :param args:                    getParser output parameter
    :param config_file_map:         getParser outputs paths to all required samples for each sample
    :param coorPath:                The coordinates of the collinearity output by SynDiv_c coor on qurey
    :param no_syn_alignmentPath:    config file path for no_syn_alignment output
    :param workDir:                 work path
    :return: calPath
    """
    # mkdir
    makedir(workDir)
    workDir = workDir + os.sep
    # cd
    os.chdir(workDir)

    # Parameter information
    filePathList = []
    sampleNameList = []
    for sampleName, value in config_file_map.items():
        filePathList.append(value["aligns"])
        sampleNameList.append(sampleName)

    # no_syn command
    cmd = f"{code_dir}SynDiv_c cal -t {args.threads} -r {args.reference} --coor {coorPath} --syri_outs {no_syn_alignmentPath} --buffer {args.buffer} "

    if args.mode == "fast":
        cmd += "--fast "
    
    calPath = os.path.join(workDir, f"{args.prefix}cal.out")
    cmd += "--aligns " + " ".join(filePathList) + " -n " + " ".join(sampleNameList) + f" -o {calPath}"

    logger.error(f'CMD: {cmd}')

    # submit task
    run_command(cmd, code_path)

    return calPath


# Run SynDiv_c cal
def window(args, calPath):
    """
    :param args:  getParser output parameter
    :param calPath:  syntenic diversity for each coordinate calculated by SynDiv_c cal
    :return: windowPath
    """
    # output file
    windowPath = os.path.abspath(f"{args.prefix}SynDiv.out")

    # no_syn command
    cmd = f"{code_dir}SynDiv_c window -r {args.reference} -i {calPath} -w {args.window} -s {args.step} -o {windowPath}"

    logger.error(f'CMD: {cmd}')

    # submit task
    run_command(cmd, code_path)

    return windowPath


def main():
    # ################################################# parsing parameters ################################################# #
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

    # ################################################# cal ################################################# #
    os.chdir(args.dir)
    calPath = cal(args, config_file_map, coorPath, no_syn_alignmentPath, os.path.join(args.dir, "cal"))
    
    # ################################################# window ################################################# #
    os.chdir(args.dir)
    winPath = window(args, calPath)

    logger.error(f'Result: {calPath}')
    logger.error(f'Result: {winPath}')
    logger.error(f'Done.')


if __name__ == '__main__':
    main()
