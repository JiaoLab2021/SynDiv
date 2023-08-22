#!/usr/bin/env python3

# -*- coding: utf-8 -*-


import os
import shutil
import subprocess
import multiprocessing
import logging
import pandas as pd
import syri_syn


# Define the log format
logger = logging.getLogger('no_syn_alignment')
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
    # submit task
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=envPath)

    # get exit status
    exit_state = proc.wait()

    stdout, stderr = proc.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()

    # standard output and error output
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
    Non-collinear coordinates that need to be compared
    """
    def __init__(
            self,
            ali_loci_filename: str, 
            config_file_map: dict, 
            args, 
            code_path
    ):
        """
        Constructor
        :param ali_loci_filename: Coordinates to be compared, the output of SynDic_c no_syn
        :param config_file_map: the out put of getParser
        :param args: the output of getParser
        """
        # Coordinate files to be compared
        self._ali_loci_filename = ali_loci_filename

        # Genome File Dictionary
        self._config_file_map = config_file_map  # dict{sample, {genome: "path", aligns: "path", syri_out: "path"}}

        # Input parameters
        self._args = args

        # environment variable
        self._path = code_path

        # Coordinates to be compared
        self._long_ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >, row with too many columns
        self._short_ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >, rows with fewer columns

        # syntenic result
        self._sample1_sample2_synpd = {}  # map<sample1, map<sample2, pd.DataFrame()> >


    def _index(
        self
    ):
        """
        Build coordinate index
        """
        # All coordinates that need to be compared
        ali_loci_dict = {}  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >

        # Number of samples
        sample_num = 0

        with open(self._ali_loci_filename) as f:
            for infos in f.readlines():
                # skip blank line
                if len(infos) == 0:
                    continue

                # go string and split
                infos_list = infos.strip().split()

                if len(infos_list) - 2 > sample_num:
                    sample_num = len(infos_list) - 2
                
                # If there is only one coordinate, skip that line
                if len(infos_list) <= 3:
                    continue

                # chromosome
                chr_tmp = infos_list[0]

                # temporary key
                key_tmp = f"{infos_list[0]}_{infos_list[1]}"

                # initialize dictionary
                ali_loci_dict[key_tmp] = {}

                # Loop from the third column
                for idx1 in range(2, len(infos_list)):
                    info = infos_list[idx1]

                    info_list = info.strip().split(":")
                    sample_tmp = info_list[0]

                    loci_list = info_list[1].split(";")
                    
                    # initialize dictionary
                    ali_loci_dict[key_tmp][sample_tmp] = []

                    # Skip if only one coordinate. record coordinates   map<sample, vector<chr:start-end>
                    ali_loci_dict[key_tmp][sample_tmp] = [f"{chr_tmp}:{loci}" for loci in loci_list if len(loci.split("-")) == 2]

        # Split ali_loci_dict into two dictionaries with long and short columns using dictionary comprehensions and conditional expressions
        threshold = sample_num * 0.8
        self._long_ali_loci_dict = {k: v for k, v in ali_loci_dict.items() if len(v) > threshold}
        self._short_ali_loci_dict = {k: v for k, v in ali_loci_dict.items() if len(v) <= threshold}
        # print quantity
        logger.error(f'Number of tasks running in parallel by row: {len(self._short_ali_loci_dict)}')
        logger.error(f'Number of tasks running in parallel by column: {len(self._long_ali_loci_dict)}')



    def _alignment_line(
        self, 
        thread_indice
        
    ):
        """
        A multi-threaded function for extracting sequences and aligning them
        : param key_value   reference genome coordinates_map<sample, vector<chr:start-end> >   map<sample, vector<chr:start-end> >
        : param key1        reference genome coordinates, chr+"_"+start
        : param value1      map<sample, vector<chr:start-end> >
        :
        : return [key1, sample1_sample2_synpd, map<sample1, map<sample2, pd.DataFrame()> >]
        """
        # Store calculation results
        loci_sample1_sample2_synpd_list = []

        # Record iteration value
        index = -1
        for key1, value1 in thread_indice[2].items():
            # Iterator increment
            index += 1
            
            # Determines whether the index is in the thread's index
            if index < thread_indice[0]:
                continue
            if index > thread_indice[1]:
                # Returns the collinear coordinates of all sample combinations at this threads
                return loci_sample1_sample2_synpd_list

            # temporarily store collinear coordinates
            sample1_sample2_synpd = {}  # map<sample1, map<sample2, pd.DataFrame()> >

            # Store extracted sequence files
            file_name_map = {key2: os.path.abspath(f"{key1}_{key2}.fa") 
                        for key2 in value1.keys()}
            file_name_loc_map = {key2: {} for key2 in value1.keys()}
            sample_list = list(value1.keys())

            # chromosome
            chromosome = ""

            # ############### subseq ############### #
            for key2, value2 in value1.items():  # map<sample, vector<chr:start-end> >
                # Where the extracted sequences are saved
                output_file_name = file_name_map[key2]

                # submit task
                chromosome, file_name_loc_map_tmp = samtools(
                    self._path, 
                    self._ali_loci_filename,
                    key1, 
                    key2, 
                    value2, 
                    self._config_file_map, 
                    output_file_name
                )

                # Record the coordinate information of the sample
                for sample2_tmp, loci_map in file_name_loc_map_tmp.items():
                    file_name_loc_map[sample2_tmp] = loci_map


            # ############### minimap2 + syri ############### #
            # Record the sample collinear with sample1, which is used to judge whether the subsequent samples need to be compared with each other
            syn_sample_dict = {}  # map<sample1, map<sample2, 0/1> >   0->noSyn  1->Syn

            # sample1
            for idx1, sample1 in enumerate(sample_list):
                filename1 = file_name_map[sample1]

                # Record the sample collinear with sample1, which is used to judge whether the subsequent samples need to be compared with each other
                synSampleList = []
                # Record the non-collinear sample with sample1, which is used to judge whether the subsequent samples need to be compared with each other
                noSynSampleList = []

                # ref coordinates and length information
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
                
                # sample2
                for idx2 in range(idx1 + 1, len(sample_list)):
                    sample2 = sample_list[idx2]
                    filename2 = file_name_map[sample2]

                    qryLenSum = 0
                    for key3, value3 in file_name_loc_map[sample2].items():  # map<start, end>
                        qryLenTmp = abs(value3 - key3 + 1)
                        qryLenSum += qryLenTmp

                    # Determine whether to compare and calculate collinearity
                    if sample1 in syn_sample_dict:
                        if sample2 in syn_sample_dict[sample1]:
                            # 1 means self._args.synRatio is collinear, add directly
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
                                    # assignment
                                    sample1_sample2_synpd[sample1][sample2] = SynPDTmp
                                except KeyError:
                                    # initialization
                                    if sample1 not in sample1_sample2_synpd:
                                        sample1_sample2_synpd[sample1] = {}
                                    if sample2 not in sample1_sample2_synpd[sample1]:
                                        sample1_sample2_synpd[sample1][sample2] = pd.DataFrame()
                                    # assignment
                                    sample1_sample2_synpd[sample1][sample2] = SynPDTmp
                                continue
                            # 0 means non-collinearity, also no comparison, just skip the sample2
                            else:
                                continue
                    
                    # output filename
                    alignment_file_name = os.path.abspath(f'{key1}_{sample1}_{sample2}.sam')

                    # submit task
                    # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
                    synpd = minimap2_syn(self._args, self._path, filename1, filename2, alignment_file_name, ref_start_tmp_list, ref_end_tmp_list, ref_len_tmp_list)

                    # If it is empty, go directly to the next loop
                    if synpd.empty:
                        # Samples recorded as non-collinear
                        noSynSampleList.append(sample2)
                        continue

                    try:
                        # assignment
                        sample1_sample2_synpd[sample1][sample2] = synpd
                    except KeyError:
                        # initialization
                        if sample1 not in sample1_sample2_synpd:
                            sample1_sample2_synpd[sample1] = {}
                        if sample2 not in sample1_sample2_synpd[sample1]:
                            sample1_sample2_synpd[sample1][sample2] = pd.DataFrame()
                        # assignment
                        sample1_sample2_synpd[sample1][sample2] = synpd

                    # Add a list here to collect all the samples that are collinear, and skip these 
                    # samples for subsequent comparisons. The length of collinearity accounts for more than self._args.synRatio
                    qryTotal = synpd['bLen'].sum()
                    alignRatio = float(qryTotal)/max(qryLenSum, 1)
                    if alignRatio >= self._args.synRatio:
                        synSampleList.append(sample2)
                    elif alignRatio <= self._args.nosynRatio:
                        # Samples recorded as non-collinear
                        noSynSampleList.append(sample2)

                # Add collinear samples to the total dict for later comparison queries
                # sample1
                for idx2, sample2 in enumerate(synSampleList):
                    # initialization dict
                    if sample2 not in syn_sample_dict:
                        syn_sample_dict[sample2] = {}
                    # sample2
                    for idx3 in range(idx2 + 1, len(synSampleList)):
                        syn_sample_dict[sample2][synSampleList[idx3]] = 1
                # Add non-collinear samples to the total dict for later comparison queries
                # sample1
                for sample2 in synSampleList:
                    # initialization dict
                    if sample2 not in syn_sample_dict:
                        syn_sample_dict[sample2] = {}
                    # sample2
                    for sample3 in noSynSampleList:
                        syn_sample_dict[sample2][sample3] = 0

            # remove file (fasta)
            for key2, value2 in file_name_map.items():
                if os.path.exists(value2) and self._args.debug==False:
                    os.remove(value2)

            # The collinear coordinates of all sample combinations at the current position
            loci_sample1_sample2_synpd_list.append((key1, sample1_sample2_synpd))

        # Returns the collinear coordinates of all sample combinations at this threads
        return loci_sample1_sample2_synpd_list
    

    def _alignment_column(
        self, 
        key1, 
        value1, 
        pool
    ):
        """
        A multi-threaded function for extracting sequences and aligning them
        : param key1    reference genome coordinates, map<chr+"_"+start, map<sample, vector<chr:start-end> > >
        : param value1  map<sample, vector<chr:start-end> >
        : param pool    process pool
        :
        : return key1, sample1_sample2_synpd_tmp, map<sample1, map<sample2, pd.DataFrame()> >
        """
        # temporarily store collinear coordinates
        sample1_sample2_synpd_tmp = {}  # map<sample1, map<sample2, pd.DataFrame()> >

        # Store extracted sequence files
        file_name_map = {key2: os.path.abspath(f"{key1}_{key2}.fa") 
                    for key2 in value1.keys()}
        file_name_loc_map = {key2: {} for key2 in value1.keys()}
        sample_list = list(value1.keys())

        # ############### subseq ############### #
        results = []
        for key2, value2 in value1.items():  # map<sample, vector<chr:start-end> >
            # Where the extracted sequences are saved
            output_file_name = file_name_map[key2]

            # submit task
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

        # Get the return value of the function
        for result in results:
            chromosome, file_name_loc_map_tmp = result.get()

            # Record the coordinate information of the sample
            for sample2_tmp, loci_map in file_name_loc_map_tmp.items():
                file_name_loc_map[sample2_tmp] = loci_map


        # ############### minimap2 + syri ############### #
        # Store the final result of the comparison
        results = []
        
        # sample1
        for idx1, sample1 in enumerate(sample_list):
            filename1 = file_name_map[sample1]
            
            # ref coordinates and length information
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
            
            # sample2
            for idx2 in range(idx1 + 1, len(sample_list)):
                sample2 = sample_list[idx2]
                filename2 = file_name_map[sample2]
                
                # output filename
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

        # Get the return value of the function
        for result in results:
            # Get the output
            sample1, sample2, synpd = result
            synpd = synpd.get()
            
            # If it is empty, go directly to the next loop
            if synpd.empty:
                continue

            try:
                # assignment
                sample1_sample2_synpd_tmp[sample1][sample2] = synpd
            except KeyError:
                # initialization
                if sample1 not in sample1_sample2_synpd_tmp:
                    sample1_sample2_synpd_tmp[sample1] = {}
                if sample2 not in sample1_sample2_synpd_tmp[sample1]:
                    sample1_sample2_synpd_tmp[sample1][sample2] = pd.DataFrame()
                # assignment
                sample1_sample2_synpd_tmp[sample1][sample2] = synpd

        # remove file
        for key2, value2 in file_name_map.items():
            if os.path.exists(value2) and self._args.debug==False:
                os.remove(value2)

        # Returns the collinear coordinates of all sample combinations at the current position
        return key1, sample1_sample2_synpd_tmp
    

def samtools_index(
    path,
    config_file_map
):
    """
    :description:              build the genome index
    :param path:               environment variable
    :param: config_file_map    dict{sample, {genome: "path", aligns: "path", syri_out: "path"}} / dict{sample, genome_path}
    :
    :return:                   0
    """
    for key1, value1 in config_file_map.items():
        # the path of genome
        genomeFile = ""
        # Profiles submitted by SynDiv have ["genome"] key
        try:
            genomeFile = value1['genome']
        # Profiles submitted by SynDiv_p don't have ["genome"] key
        except TypeError:
            genomeFile = value1

        # the path of genome idx
        genomeIdxFile = f"{genomeFile}.fai"
        
        # Check if the index file exists
        if os.path.exists(genomeIdxFile):
            continue
        else:
            # build the index
            cmd = f"samtools faidx {genomeFile}"
            # submit task
            run_command(cmd, path)

    return 0


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
    :description:                 extract sequence
    :param path:                  environment variable
    :param ali_loci_filename:     SynDic_c no_syn, The output coordinates that need to be compared
    :param key1:                  chr+"_"+start
    :param key2:                  sample
    :param value2:                vector<chr:start-end>
    :param config_file_map:         input configuration file dict
    :param output_file_name:        output filename
    :return file_name_loc_map_tmp:    map<sample2, map<start, end> >
    """
    # ########## extract sequence ######### #
    # Record the coordinate information of sample2
    file_name_loc_map_tmp = {}
    file_name_loc_map_tmp[key2] = {}  # map<sample2, map<start, end> >

    for idx1, loci in enumerate(value2):
        # Check whether the chromosome number is in compliance with the regulations,  (chr1:257130-257225)
        if ":" not in loci or "-" not in loci:
            raise ValueError(f"Error: incorrect '{loci}' found in '{ali_loci_filename}'.")
        
        # Extract chromosome number
        chromosome = loci.split(":")[0]

        cmd = ""

        # If it is the first sequence, keep the chromosome number
        if idx1 == 0:
            # Profiles submitted by SynDiv have ["genome"] key
            try:
                cmd = f'echo \">{key1.split("_")[0]}\" > {output_file_name}; samtools faidx {config_file_map[key2]["genome"]} {loci} | grep -v \">\" >> {output_file_name}'
            # Profiles submitted by SynDiv_p don't have ["genome"] key
            except TypeError:
                cmd = f'echo \">{key1.split("_")[0]}\" > {output_file_name}; samtools faidx {config_file_map[key2]} {loci} | grep -v \">\" >> {output_file_name}'
        else:  # If it is a subsequent sequence, do not need the chromosome number
            # Profiles submitted by SynDiv have ["genome"] key
            try:
                cmd = f'samtools faidx {config_file_map[key2]["genome"]} {loci} | grep -v \">\" >> {output_file_name}'
            # Profiles submitted by SynDiv_p don't have ["genome"] key
            except TypeError:
                cmd = f'samtools faidx {config_file_map[key2]} {loci} | grep -v \">\" >> {output_file_name}'

        # record coordinates
        lociList = loci.split(":")[1].split("-")
        file_name_loc_map_tmp[key2][int(lociList[0])] = int(lociList[1])

        # submit task
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
    :description:                 Align and find collinear coordinates
    :param args:                  Input parameter information
    :param path:                  environment variable
    :param filename1:             the name of sample1
    :param filename2:             the name of sample2
    :param alignment_file_name:     the output of alignment
    :param ref_start_tmp_list:       List of starting coordinates for sample 1
    :param ref_end_tmp_list:         List of ending coordinates for sample 1
    :param ref_len_tmp_list:         List of lengths for sample 1
    """
    # minimap2
    cmd = f'minimap2 -ax asm5 --eqx {filename1} {filename2} -o {alignment_file_name}'
    # submit task
    run_command(cmd, path)

    # If the file does not exist, return directly
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
    :description:                              Parsing minimap2 results, looking for collinear coordinates
    :param args:                               Input parameter information
    :param ref_start_tmp_list:                 ref coordinates and length information
    :param ref_end_tmp_list:                   ref coordinates and length information
    :param ref_len_tmp_list:                   ref coordinates and length information
    :param alignment_file_name:                minimap2 Comparison result file
    :return: synpd                             synpd -> pd.DataFrame()
    """
    # temporarily store collinear coordinates
    # ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']
    synpd = pd.DataFrame()

    # syri, Look for Collinearity Information
    synDatas = pd.DataFrame()
    try:
        synDatas = syri_syn.main(args, alignment_file_name)
    except:
        return synpd

    # remove file
    if os.path.exists(alignment_file_name) and args.debug==False:
        os.remove(alignment_file_name)

    # If synDatas is empty, return a null value directly
    if synDatas.empty:
        return synpd
    
    # debug code
    if args.debug:
        logger.error(f"alignment_file_name:{alignment_file_name}\nsynDatas:{synDatas}")

    # Save the coordinates to table
    # Print the elements of each row by row
    # Used to determine how to calculate coordinates
    idxTmp = 0  # temporary index
    preLenTmp = 0  # The length (n) of the preceding sequence is used to determine whether it is an independent alignment segment.
    thisLenTmp = ref_len_tmp_list[0]  # The maximum value of the sum of the length of the preceding sequence and the next index (n+1) represents that there is no further sequence below.
    for index, row in synDatas.iterrows():
        try:
            # temporary DataFrame
            synDataTmp1 = synDatas.iloc[[index]].copy()

            if args.debug:
                logger.error(f"synDataTmp1:{synDataTmp1}")

            # If it is an empty set or a reverse comparison, the next cycle
            if synDataTmp1.empty or synDataTmp1.iloc[0, 7] == -1 or synDataTmp1.iloc[0, 8] == -1:
                continue

            # Judgment coordinate position
            while True and idxTmp < len(ref_len_tmp_list):
                # If the start is greater than the length of all current sequences, update the index and length
                while synDataTmp1.iloc[0, 0] > thisLenTmp and idxTmp < len(ref_len_tmp_list):
                    idxTmp += 1
                    preLenTmp = sum(ref_len_tmp_list[:idxTmp])

                    # If there are still coordinates below, update the next length
                    if idxTmp + 1 < len(ref_len_tmp_list):
                        thisLenTmp = sum(ref_len_tmp_list[:idxTmp + 1])
                    else:  # If not, assign the maximum value
                        thisLenTmp = sum(ref_len_tmp_list)

                # is an independent comparison
                if synDataTmp1.iloc[0, 0] >= preLenTmp and synDataTmp1.iloc[0, 1] <= thisLenTmp:
                    synDataTmp1.iloc[0, 0] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                    synDataTmp1.iloc[0, 1] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                    try:
                        # assignment
                        synpd = pd.concat([synpd, synDataTmp1], ignore_index=True)
                    except KeyError:
                        # assignment
                        synpd = pd.concat([synpd, synDataTmp1], ignore_index=True)

                    if args.debug:
                        logger.error(f"\n{synDataTmp1}\n")

                    break
                # spans two aligned sequences
                else:
                    # Extract and save the alignment coordinates of the current sequence
                    synDataTmp2 = synDataTmp1.copy()
                    synDataTmp2.iloc[0, 0] += ref_start_tmp_list[idxTmp] - 1 - preLenTmp
                    synDataTmp2.iloc[0, 1] = ref_end_tmp_list[idxTmp]
                    try:
                        # assignment
                        synpd = pd.concat([synpd, synDataTmp2], ignore_index=True)
                    except KeyError:
                        # assignment
                        synpd = pd.concat([synpd, synDataTmp2], ignore_index=True)
                    if args.debug:
                        logger.error(f"\n{synDataTmp2}\n")
                    # Update the coordinates of synDataTmp1 and calculate the remaining coordinates
                    synDataTmp1.iloc[0, 0] = thisLenTmp + 1
        except KeyError as e:  # handle exception
            if args.debug:
                logger.error(f"KeyError occurred: {e}")
            continue
        except IndexError as e:  # handle exception
            if args.debug:
                logger.error(f"IndexError occurred: {e}")
                logger.error(f'Error ref_start_tmp_list: {ref_start_tmp_list}')
                logger.error(f"Error ref_end_tmp_list: {ref_end_tmp_list}")
                logger.error(f"Error ref_len_tmp_list: {ref_len_tmp_list}")
                logger.error(f"Error index:{index}", index)
                logger.error(f"Error alignment_file_name: {alignment_file_name}")
            continue
    
    return synpd


# save the result to a file
def save_result(sample1, sample2, synpd):
    """
    :param sample1:     the name of sample1
    :param sample2:     the name of sample2
    :param synpd:      共线性的坐标   pd.DataFrame()
    :
    :return: 0
    """
    # save the output string
    outTxt = ""

    # write the result to the output string
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

    # output to file
    outputFilePath = os.path.abspath("{}_{}.syn.out".format(sample1, sample2))
    with open(outputFilePath, 'a', encoding="utf-8") as outputFile:
        outputFile.write(outTxt)

    return 0


# save the result to a file
def sort_result(filePath):
    """
    :param filePath:    The path of the file to be sorted
    :
    :return: 0
    """
    # open the file and read each line
    with open(filePath, 'r') as f:
        lines = f.readlines()
    
    # Sort and output by chromosome and start position
    sorted_lines = sorted(lines, key=lambda x: (x.split()[0], int(x.split()[1])))
    
    # write the sorted lines to a new file
    with open(filePath, 'w') as f:
        f.writelines(sorted_lines)

    return 0


# Calculates the left and right indexes of the thread pool
def calculate_indices(my_dict):
    # Total number of tasks
    dict_length = len(my_dict)

    # Number of tasks per process
    if dict_length > 10000:
        chunk_size = 100
    elif dict_length > 5000:
        chunk_size = 50
    else:
        chunk_size = 15
    
    import math
    jobs_num = math.ceil(dict_length / chunk_size)

    thread_indices = []

    for i in range(jobs_num):
        left_index = i * chunk_size
        right_index = min((i + 1) * chunk_size - 1, dict_length - 1)
        thread_indices.append((left_index, right_index, my_dict))

    return thread_indices


def main(args, config_file_map, no_synPath, workDir, code_path):
    """
    :param args:              the output of getParser
    :param config_file_map:   getParser outputs paths to all required samples for each sample
    :param no_synPath:        The non-collinear coordinates output by SynDiv_c no_syn that need to be compared
    :param workDir:           work path
    :param code_path:         environment variable
    :
    :return: coorPath
    """
    # mkdir
    makedir(workDir)
    workDir = workDir + os.sep

    # cd
    os.chdir(workDir)

    
    # Constructor
    ali_loci_class = AlignmentLociClass(no_synPath, config_file_map, args, code_path)


    # #################################### build index #################################### #
    logger.error(f'Building index.')

    ali_loci_class._index()

    samtools_index(code_path, config_file_map)


    # #################################### Alignment #################################### #
    logger.error(f'Alignment started.')

    # Create a process pool, specify the maximum number of processes as args.jobs
    pool = multiprocessing.Pool(processes=args.jobs)

    # task counter initialization
    task_count = 0
    task_total = len(ali_loci_class._long_ali_loci_dict) + len(ali_loci_class._short_ali_loci_dict)

    # Traverse ali_loci_class._short_ali_loci_dict dictionary, multi-process
    if len(ali_loci_class._short_ali_loci_dict) > 0:
        thread_indices = calculate_indices(ali_loci_class._short_ali_loci_dict)
        results = pool.imap_unordered(ali_loci_class._alignment_line, thread_indices)

        # #################################### Save #################################### #
        for loci_sample1_sample2_synpd_list in results:
            for loci_sample1_sample2_synpd in loci_sample1_sample2_synpd_list:
                loci, result = loci_sample1_sample2_synpd
                # After each task is completed, increment the task counter by 1
                task_count += 1
                # print progress bar
                progress = task_count / task_total * 100
                logger.error(f'Alignment Progress: {progress:.2f}% | {task_count}/{task_total} | {loci}')

                # save result
                for key1, value1 in result.items():
                    for key2, value2 in value1.items():
                        save_result(key1, key2, value2)

    # Traverse ali_loci_class._long_ali_loci_dict dictionary, multi-process
    for key1, value1 in ali_loci_class._long_ali_loci_dict.items():  # map<chr+"_"+start, map<sample, vector<chr:start-end> > >
        loci, result = ali_loci_class._alignment_column(key1, value1, pool)

        # After each task is completed, increment the task counter by 1
        task_count += 1
        # print progress bar
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

    # Close the process pool and wait for all tasks to complete
    pool.close()
    pool.join()


    # #################################### Output profile information #################################### #
    no_syn_alignmentPath = os.path.abspath(f"{args.prefix}no_syn_alignment.out")
    with open(no_syn_alignmentPath, 'w', encoding="utf-8") as f:
        for synOutFilePath in synOutFilePathList:
            fileNameList = os.path.basename(synOutFilePath).replace(".syn.out", "").split("_")
            f.write('\t'.join([fileNameList[0], fileNameList[1], os.path.abspath(synOutFilePath)]) + "\n")

    return no_syn_alignmentPath
