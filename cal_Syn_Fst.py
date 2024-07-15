#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2024/06/18"
__version__ = "1.0.2"
__author__ = "jbhe"
__email__ = "hejiabao@webmail.hzau.edu.cn or hejiabao2001@gmail.com"


# Importing libraries
import argparse
import os
import subprocess


# environment variables
env_path = {'PATH': os.environ.get('PATH')}

# Define a function to parse command-line arguments
def get_parser():
    # Create an ArgumentParser object for parsing command-line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    # Set description for the script
    parser.description = (
        "## To specify the populations configuration, use the following format in your file:\n"
        "# sub1 - Subpopulation 1\n"
        "# sub2 - Subpopulation 2\n"
        "# pop - Total population\n"
        "    sub1\tsub1_count\tsub1.cal.out\tsub1.cal.out.ins\tsubpopulation\n"
        "    sub2\tsub2_count\tsub2.cal.out\tsub2.cal.out.ins\tsubpopulation\n"
        "    pop\tpop_count\tpop.cal.out\tpop.cal.out.ins\tpopulation\n"
    )


    # Input parameters
    required_input = parser.add_argument_group("Input arguments")
    required_input.add_argument(
        "-c", "--config",
        dest="config",
        help="configuration file (format: name count cal_out_path cal_out_ins_path subpopulation/population)",
        required=True
    )

    # Output parameters
    required_out = parser.add_argument_group("Input arguments")
    required_out.add_argument(
        "-p", "--prefix",
        dest="prefix",
        help="output file prefix",
        required=True
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Get parsed input and output filenames
    config = args.config
    prefix = args.prefix

    return config, prefix

# Calling the command line
def run(command, envPath):
    command = subprocess.run(command, shell=True, check=True, text=True, capture_output=True, env=envPath)
    return command.stdout

# Read synteny diversity of populations and subpopulations
def read_sub(config):
    # Open the subs configuration file
    subs_configuration_dict = {}
    subs_ins_configuration_dict = {}
    pop_configuration_list = []
    pop_ins_configuration_list = []
    try:
        with open(config, 'r') as f:
            samplesets_list = f.readlines()
    except IOError as e:
        raise ValueError(f'Error occurred while reading the file {samplesets_list}: {e}')

    if not samplesets_list:
        raise ValueError(f'Empty file: {samplesets_list}')
    for sampleset in samplesets_list:
        if not sampleset.strip() or "#" in sampleset:
            continue

        sampleset_list = sampleset.split()
        if len(sampleset_list) == 5:
            sampleset_name = str(sampleset_list[0])
            sampleset_count = int(sampleset_list[1])
            sampleset_cal = os.path.abspath(sampleset_list[2])
            sampleset_cal_ins = os.path.abspath(sampleset_list[3])
            sampleset_type = str(sampleset_list[4])

            try:
                # sub
                if sampleset_type == "subpopulation":
                    subs_configuration_dict[sampleset_name] = [sampleset_count, sampleset_cal]
                    subs_ins_configuration_dict[sampleset_name] = [sampleset_count, sampleset_cal_ins]
                # pop
                elif sampleset_type == "population":
                    pop_configuration_list = [sampleset_name, sampleset_count, sampleset_cal]
                    pop_ins_configuration_list = [sampleset_name, sampleset_count, sampleset_cal_ins]
            except Exception as e:
                raise Exception(f'Error: The fourth column of the configuration file is wrong in {config}: {e}')

            # Check if files exist
            if not os.path.exists(sampleset_cal) or os.path.getsize(sampleset_cal) == 0:
                raise ValueError(f"'{sampleset_cal}': No such file or is empty.")

            if not os.path.exists(sampleset_cal_ins) or os.path.getsize(sampleset_cal_ins) == 0:
                raise ValueError(f"'{sampleset_cal_ins}': No such file or is empty.")

    return subs_configuration_dict, subs_ins_configuration_dict, pop_configuration_list, pop_ins_configuration_list

def ins_merge(subs_configuration_dict, subs_ins_configuration_dict, pop_configuration_list, pop_ins_configuration_list):
    # Read the synteny diversity of all subpopulations
    for sub_name, sub_value in subs_configuration_dict.items():
        sub_cal = os.path.abspath(sub_value[-1])
        sub_ins_cal = os.path.abspath(subs_ins_configuration_dict[sub_name][-1])
        sub_ins_tmp_file = f'{sub_name}.cal.out.ins.tmp'
        sub_command = (f'awk \'BEGIN {{ FS = "\t"; OFS = "\t" }} '
                       f'NR == FNR && !/^#/ {{ b[$1"\t"$2] = $3; next }} '
                       f'!/^#/ {{ key = $1"\t"$2; if (key in b) {{ $5 = b[key]; print $5,"INS"; }} else '
                       f'{{ print $5, "ALL"; }} }}\' {sub_ins_cal} {sub_cal} | sed \'1 i '
                       f'{sub_name}_INS_Syntenic_Diversity\t{sub_name}_INS_Syntenic_Diversity\n\' > {sub_ins_tmp_file}'
                       )
        run(sub_command, env_path)
    # Open the pop file
    pop_name = pop_ins_configuration_list[0]
    pop_cal = os.path.abspath(pop_configuration_list[-1])
    pop_ins_cal = os.path.abspath(pop_ins_configuration_list[-1])
    pop_tmp_file = f'{pop_name}.cal.out.ins.tmp'
    pop_command = (f'awk \'BEGIN {{ FS = "\t"; OFS = "\t" }} '
                   f'NR == FNR && !/^#/ {{ b[$1"\t"$2] = $3; next }} '
                   f'!/^#/ {{ key = $1"\t"$2; if (key in b) {{ $5 = b[key]; print $1,$2,$5,"INS"; }} else '
                   f'{{ print $1,$2,$5,"ALL"; }} }}\' {pop_ins_cal} {pop_cal} | sed \'1 i '
                   f'#CHROM\tPOS\t{pop_name}_INS_Syntenic_Diversity\n\' > {pop_tmp_file}'
                   )
    run(pop_command, env_path)
    pop_ins_on_tmp_file = f'{pop_name}.cal.out.ins.on.tmp'
    pop_ins_on_command = f'awk -F "\t" \'{{OFS="\t"; print $1, $2}}\' {pop_tmp_file} > {pop_ins_on_tmp_file}'
    run(pop_ins_on_command, env_path)
    pop_ins_down_tmp_file = f'{pop_name}.cal.out.ins.down.tmp'
    pop_ins_down_command = f'awk -F "\t" \'{{OFS="\t"}}; ' \
                       f'{{print $3,$4}}\' {pop_tmp_file} > {pop_ins_down_tmp_file} && rm {pop_tmp_file}'
    run(pop_ins_down_command, env_path)

    # Merge
    file_list = [pop_ins_on_tmp_file] + [f"{sub_name}.cal.out.ins.tmp" for sub_name in subs_configuration_dict.keys()] \
                + [pop_ins_down_tmp_file]
    all_cal_tmp_file = f'{pop_name}.cal.out.ins.all.tmp'
    paste_command = f'paste -d "\t" {" ".join(file_list)} > {all_cal_tmp_file} && rm  {" ".join(file_list)}'
    run(paste_command, env_path)
    modify_command = f'awk -F "\t" \'BEGIN {{ OFS = "\t" }} /INS/ {{ for (i = 4; i <= NF; i += 2) $i = ""; sub(/\t$/' \
                     f', ""); print }}\' {all_cal_tmp_file} > {all_cal_tmp_file}.modify && rm {all_cal_tmp_file}'
    run(modify_command, env_path)
    remove_command = (f'awk -F "\\t" \'{{ OFS = "\\t"; for (i=1; i<=NF; i++) if ($i != "") printf "%s%s", $i, '
                      f'(i==NF ? "\\n" : OFS) }}\' {all_cal_tmp_file}.modify > {all_cal_tmp_file} && sed -i \'1 i '
                      f'#Insertion\' {all_cal_tmp_file} && rm {all_cal_tmp_file}.modify'
    )
    run(remove_command, env_path)
    all_cal_ins_tmp_file_path = os.path.abspath(f'{all_cal_tmp_file}')

    return all_cal_ins_tmp_file_path

def all_merge(subs_configuration_dict, pop_configuration_list):
    # Read the synteny diversity of all subpopulations
    for sub_name, sub_value in subs_configuration_dict.items():
        sub_cal = os.path.abspath(sub_value[-1])
        sub_tmp_file = f'{sub_name}.cal.out.tmp'
        sub_command = f'awk -F "\t" \'{{OFS="\t"}} NR == 1 {{print "{sub_name}_Syntenic_Diversity"}} ' \
                      f'NR > 1 {{print $5}}\' {sub_cal} > {sub_tmp_file}'
        run(sub_command, env_path)
    # Open the pop file
    pop_name = pop_configuration_list[0]
    pop_cal = os.path.abspath(pop_configuration_list[-1])
    pop_on_tmp_file = f'{pop_name}.cal.out.on.tmp'
    pop_on_command = f'awk -F "\t" \'{{OFS="\t"}} {{print $1, $2}}\' {pop_cal} > {pop_on_tmp_file}'
    run(pop_on_command, env_path)
    pop_down_tmp_file = f'{pop_name}.cal.out.down.tmp'
    pop_down_command = f'awk -F "\t" \'{{OFS="\t"}} NR == 1 {{print "{pop_name}_Syntenic_Diversity"}}; ' \
                       f'NR > 1 {{print $5}}\' {pop_cal} > {pop_down_tmp_file}'
    run(pop_down_command, env_path)
    # Merge
    file_list = [pop_on_tmp_file] + [f"{sub_name}.cal.out.tmp" for sub_name in subs_configuration_dict.keys()] \
                + [pop_down_tmp_file]
    all_cal_tmp_file = f'{pop_name}.cal.out.all.tmp'
    paste_command = f'paste -d "\t" {" ".join(file_list)} > {all_cal_tmp_file} && rm  {" ".join(file_list)}'
    run(paste_command, env_path)
    all_cal_tmp_file_path = os.path.abspath(f'{all_cal_tmp_file}')

    return all_cal_tmp_file_path

def output_name(cal_tmp_file_path):
    # Output the filtered results
    with open(cal_tmp_file_path, 'r') as file:
        first_line = file.readline()
        if '#Insertion' in first_line:
            return "INS"
    return "ALL"

def syn_fst(subs_configuration_dict, cal_tmp_file_path, prefix):
    # Output the filtered results
    output_type = output_name(cal_tmp_file_path)
    if output_type == "ALL":
        all_syn_fst_file = f'{prefix}.syn_Fst.out'
        out_syn_file = open(all_syn_fst_file, 'wt')
    elif output_type == "INS":
        all_syn_fst_file = f'{prefix}.syn_Fst.out.ins'
        out_syn_file = open(all_syn_fst_file, 'wt')
    else:
        raise ValueError(f"Unsupported output type: {output_type}")
    # Start calculating
    with open(cal_tmp_file_path, 'rt', encoding='utf-8') as f:
        # Read the file content line by line
        for line in f:
            # Split each line into a list using tabs as separators
            line_list = line.strip().split("\t")
            if "#" in line:
                if "Insertion" in line:
                    out_syn_file.write(line)
                elif "CHROM" in line:
                    # If the line contains the "#" character, it is a header line; add identifiers as header information
                    temp_list = line_list + ["Syn_Fst"]
                    temp_txt = '\t'.join(temp_list) + "\n"
                    out_syn_file.write(temp_txt)
            else:
                # Otherwise, calculate the fac value based on certain rules and check if it meets the conditions
                sub_count_list = [subs_configuration_dict[sub_name][0] for sub_name in subs_configuration_dict]
                pop_count = sum(sub_count_list)
                sub_syn_list = [float(line_list[i + 2]) for i in range(len(subs_configuration_dict))]
                pop_syn = float(line_list[-1])
                # Calculate
                pop_value = 2 * pop_syn * (1 - pop_syn)

                if pop_syn != 0 and pop_syn != 1:
                    first_numerator_first = pop_value
                    first_numerator_second = sum(2 * sub_syn_list[i] * (1 - sub_syn_list[i]) * sub_count_list[i] for i in range(len(sub_count_list))) / pop_count
                    first_numerator = first_numerator_first - first_numerator_second
                    first_denominator = pop_value
                    first_item = first_numerator / first_denominator
                else:
                    first_item = 0

                second_item = pop_value - sum(2 * sub_syn_list[i] * (1 - sub_syn_list[i]) for i in range(len(sub_count_list)))

                syn_fst = first_item + second_item

                temp_list = line_list + [str(syn_fst)]
                temp_txt = '\t'.join(temp_list) + "\n"
                out_syn_file.write(temp_txt)

    # Close the file
    out_syn_file.close()
    rm_command = f'rm {cal_tmp_file_path} &'
    run(rm_command, env_path)

# Main function
def main():
    config, prefix = get_parser()
    subs_configuration_dict, subs_ins_configuration_dict, pop_configuration_list, pop_ins_configuration_list \
        = read_sub(config)
    all_cal_ins_tmp_file_path = ins_merge(subs_configuration_dict, subs_ins_configuration_dict,
                                          pop_configuration_list, pop_ins_configuration_list)
    syn_fst(subs_configuration_dict, all_cal_ins_tmp_file_path, prefix)
    all_cal_tmp_file_path = all_merge(subs_configuration_dict, pop_configuration_list)
    syn_fst(subs_configuration_dict, all_cal_tmp_file_path, prefix)

if __name__ == '__main__':
    main()