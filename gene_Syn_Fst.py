#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2024/06/19"
__version__ = "1.0.2"
__author__ = "jbhe"
__email__ = "hejiabao@webmail.hzau.edu.cn or hejiabao2001@gmail.com"


# Import libraries
import argparse  # Import argparse module for parsing command-line arguments
import re  # Import regular expression matching module
import subprocess
import os


# environment variables
env_path = {'PATH': os.environ.get('PATH')}

# Define a function to parse command-line arguments
def get_parser():
    # Create an ArgumentParser object for parsing command-line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Input parameters
    required_input = parser.add_argument_group("Input arguments")
    required_input.add_argument(
        "-g", "--gff",
        dest="gff",
        help="Refer to the gff annotation file of reference genome",
        required=True
    )
    required_input.add_argument(
        "-syn", "--syn_fst",
        dest="syn_fst",
        help="Syn-Fst value (syn_Fst.out)",
        required=True
    )
    required_input.add_argument(
        "-syn_ins", "--syn_fst_ins",
        dest="syn_fst_ins",
        help="Syn-Fst_ins value (syn_Fst.out.ins)",
    )
    # Output parameters
    required_out = parser.add_argument_group("Output arguments")
    required_out.add_argument(
        "-p", "--prefix",
        dest="prefix",
        help="output file prefix",
        required=True
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Get parsed input and output filenames
    gff = args.gff
    syn_fst = args.syn_fst
    syn_fst_ins = args.syn_fst_ins
    prefix = args.prefix

    return gff, syn_fst, syn_fst_ins, prefix

# Calling the command line
def run(command, envPath):
    command = subprocess.run(command, shell=True, check=True, text=True, capture_output=True, env=envPath)
    return command.stdout

# Build gene index
def indexer(gff, window):
    # Create an empty dictionary
    gff_dict = {}

    # Open the input file for reading
    with open(gff, 'r', encoding='utf-8') as f:
        # Read the file content line by line
        for line in f:
            # Split each line into a list using tabs as separators
            line_list = line.strip().split("\t")
            if "#" not in line_list[0][0] and "Mt" not in line_list[0][0] and \
                    "Pt" not in line_list[0][0]:  # Skip comment lines

                # Index the line containing gene information
                if line_list[2] in ["gene", "mRNA"]:
                    genename = re.search("ID=(?:gene:)?(\\w+)", line_list[8]).group(1)  # genename
                    gene_chr = line_list[0]  # Chromosome where the gene is located
                    # Determine the direction of the gene
                    if line_list[6] == "+":  # Gene direction is forward
                        gene_start = int(line_list[3])  # Gene start position
                        gene_end = int(line_list[4])  # Gene end position

                        gene_upstream = gene_start - window  # Gene upstream
                        gene_downstream = gene_end + window  # Gene downstream

                        if gene_upstream < 1:
                            gene_upstream = 1
                    else:  # Gene direction is reverse
                        gene_start = int(line_list[4])  # Gene start position
                        gene_end = int(line_list[3])  # Gene end position

                        gene_upstream = gene_start + window  # Gene upstream
                        gene_downstream = gene_end - window  # Gene downstream

                        if gene_downstream < 1:
                            gene_downstream = 1

                    gff_dict[genename] = {}  # Create key-value pairs
                    gff_dict[genename]['location'] = [gene_chr, gene_start, gene_end]
                    gff_dict[genename]['updownstream'] = [gene_upstream, gene_downstream]
                elif line_list[2] == "exon":
                    if line_list[6] == "+":
                        exon_start = int(line_list[3])  # Exon start position
                        exon_end = int(line_list[4])  # Exon end position
                    else:
                        exon_start = int(line_list[4])  # Exon start position
                        exon_end = int(line_list[3])  # Exon end position

                    exon_list = [exon_start, exon_end]  # Temporary list
                    gff_dict[genename].setdefault('exon', []).append(exon_list)  # Add exon information to the dictionary

    # Return gff_dict, key: genename, value: dict, dict_key: location,updownstream,exon
    return gff_dict

# syn_fst file tmp
def syn_fst_file(syn_fst):

    # Preprocess the syn_Fst.out.ins file
    syn_fst_path = os.path.abspath(syn_fst)
    syn_fst_tmp_file = os.path.abspath(f'{syn_fst}.tmp')
    awk_command = f'grep -v "#" {syn_fst_path} | awk -F "\t" \'{{OFS="\t"}} {{print $0, "ALL"}}\' > {syn_fst_tmp_file}'
    run(awk_command, env_path)
    syn_fst_tmp_path = os.path.abspath(syn_fst_tmp_file)
    syn_fst_type = "ALL"
    return syn_fst_tmp_path, syn_fst_type

# syn_fst_ins file tmp
def syn_fst_ins_file(syn_fst, syn_fst_ins):

    # Preprocess the syn_Fst.out.ins file
    syn_fst_path = os.path.abspath(syn_fst)
    syn_fst_ins_path = os.path.abspath(syn_fst_ins)
    syn_fst_ins_tmp_file = f'{syn_fst_ins}.tmp'
    awk_command = (
        f'awk \'BEGIN {{ FS = "\\t"; OFS = "\\t" }} NR == FNR && !/^#/ {{ b[$1"\\t"$2] = $3"\\t"$4"\\t"$5"\\t"$6; '
        f'next }} !/^#/ {{ key = $1"\\t"$2; if (key in b) {{ print $1,$2,b[key],"INS"; }} else {{ print $0, "ALL"; }} '
        f'}}\' {syn_fst_ins_path} {syn_fst_path} > {syn_fst_ins_tmp_file}'
    )
    run(awk_command, env_path)
    syn_fst_ins_tmp_path = os.path.abspath(syn_fst_ins_tmp_file)
    syn_fst_ins_type = "INS"
    return syn_fst_ins_tmp_path, syn_fst_ins_type

# Build index
def build_syn_fst_index(tmp_file):
    # Open
    with open(tmp_file, 'rt', encoding='utf-8') as f:
        # Define a dictionary to store file content
        syn_fst_dic = {}
        # Read the file content line by line
        for line in f:
            # Split each line into a list using tabs as separators
            line_list = line.strip().split("\t")
            if "#" in line:
                # Skip comment lines
                continue
            else:
                chrname = line_list[0]
                syn_fst = float(line_list[5])
                syn_fst_type = str(line_list[6])
                syn_fst_tuple = (syn_fst, syn_fst_type)
                # Dictionary list
                syn_fst_dic.setdefault(chrname, []).append(syn_fst_tuple)  # chrname:[tuple1, tuple2, ...]

    return syn_fst_dic

# Extract target genes
def target_gene(gff_dict, syn_fst_temp_dic, outfile, window):

    # Output the filtered results
    out = open(outfile, 'wt')
    # Write column names
    colume_name = ['#Genename', 'Chromosome', 'Upstream'+'(-'+str(window)+')', 'Downstream'+'(+'+str(window)+')',
                   'Up_score', 'Exon_score', 'Down_score', 'Average_score']
    colume_name_txt = '\t'.join(colume_name) + "\n"
    out.write(colume_name_txt)

    # Loop through each gene
    for i in gff_dict.keys():
        # Filter out mitochondrial and chloroplast genes
        if gff_dict[i]['location'][0] != "ChrPt" and gff_dict[i]['location'][0] != "ChrMt":
            if 'exon' not in gff_dict[i]:
                continue

            gene_chr = gff_dict[i]['location'][0]       # Chromosome where the gene is located
            gene_up = gff_dict[i]['updownstream'][0]    # Gene upstream
            gene_down = gff_dict[i]['updownstream'][1]  # Gene downstream
            gene_tss = gff_dict[i]['location'][1]       # Gene start point
            gene_tts = gff_dict[i]['location'][2]       # Gene end point
            gene_up_region = [gene_up, gene_tss-1]      # Gene upstream region
            gene_exon_region = gff_dict[i]['exon']      # Gene exon region
            gene_down_region = [gene_tts+1, gene_down]  # Gene downstream region
            all_region = [gene_up, gene_down]           # Gene upstream and downstream 2k region

            ## Upstream score (up_score)

            # Consider positive and negative strands
            syn_fst_up_start = min(gene_up_region[0], gene_up_region[1])
            syn_fst_up_end = max(gene_up_region[0], gene_up_region[1])

            up_region_syn_fst_tuple = syn_fst_temp_dic[gene_chr][syn_fst_up_start - 1:syn_fst_up_end]
            up_region_syn_fst = [t[0] for t in up_region_syn_fst_tuple]  # syn_fst on the region
            up_region_syn_fst_type = [t[1] for t in up_region_syn_fst_tuple]
            up_score = sum(up_region_syn_fst) / len(up_region_syn_fst)  # Average syn_fst in the region
            if "INS" in up_region_syn_fst_type:
                up_score_col = f"{up_score}(INS)"
            else:
                up_score_col = f"{up_score}(ALL)"


            ## Downstream score (down_score)

            # Consider positive and negative strands
            syn_fst_down_start = min(gene_down_region[0], gene_down_region[1])
            syn_fst_down_end = max(gene_down_region[0], gene_down_region[1])

            down_region_syn_fst_tuple = syn_fst_temp_dic[gene_chr][syn_fst_down_start - 1:syn_fst_down_end]
            down_region_syn_fst = [t[0] for t in down_region_syn_fst_tuple]  # syn_fst on the region
            down_region_syn_fst_type = [t[1] for t in down_region_syn_fst_tuple]
            down_score = sum(down_region_syn_fst) / len(down_region_syn_fst)  # Average syn_fst in the region
            if "INS" in down_region_syn_fst_type:
                down_score_col = f"{down_score}(INS)"
            else:
                down_score_col = f"{down_score}(ALL)"


            ## Exon score (exon_score)

            # Initialize exon score
            exon_score = 0
            exon_region_syn_fst_type = []
            # Loop through each exon region
            for index in range(len(gene_exon_region)):
                every_gene_exon_region = gene_exon_region[index]
                # Consider positive and negative strands
                syn_fst_exon_start = min(every_gene_exon_region[0], every_gene_exon_region[1])
                syn_fst_exon_end = max(every_gene_exon_region[0], every_gene_exon_region[1])

                sub_exon_region_syn_fst_tuple = syn_fst_temp_dic[gene_chr][syn_fst_exon_start - 1:
                                                            syn_fst_exon_end]  # syn_fst on each exon region
                sub_exon_region_syn_fst = [t[0] for t in sub_exon_region_syn_fst_tuple]  # syn_fst on the region
                sub_exon_region_syn_fst_type = [t[1] for t in sub_exon_region_syn_fst_tuple]
                sub_exon_score = sum(sub_exon_region_syn_fst) / len(sub_exon_region_syn_fst) \
                                 / len(gene_exon_region)  # Average syn_fst in each exon region
                exon_score += sub_exon_score  # Subtotal score for each exon region
                exon_region_syn_fst_type.extend(sub_exon_region_syn_fst_type)
            if "INS" in exon_region_syn_fst_type:
                exon_score_col = f"{exon_score}(INS)"
            else:
                exon_score_col = f"{exon_score}(ALL)"


            ## Average score for the entire region (average_score)

            # Consider positive and negative strands
            syn_fst_all_start = min(all_region[0], all_region[1])
            syn_fst_all_end = max(all_region[0], all_region[1])

            # syn_fst at the upstream and downstream 2k of the gene
            all_region_syn_fst_tuple = syn_fst_temp_dic[gene_chr][syn_fst_all_start - 1:syn_fst_all_end]
            all_region_syn_fst = [t[0] for t in all_region_syn_fst_tuple]  # syn_fst on the region
            all_region_syn_fst_type = [t[1] for t in all_region_syn_fst_tuple]
            average_score = sum(all_region_syn_fst) / len(all_region_syn_fst)  # Average syn_fst
            if "INS" in all_region_syn_fst_type:
                average_score_col = f"{average_score}(INS)"
            else:
                average_score_col = f"{average_score}(ALL)"

            temp_list = [str(i), str(gene_chr), str(gene_up), str(gene_down), str(up_score_col),
                         str(exon_score_col), str(down_score_col), str(average_score_col)]  # List to store the content of each line
            temp_txt = '\t'.join(temp_list) + "\n"
            out.write(temp_txt)  # Write to file

    # Close the file
    out.close()

# Main function
def main():
    # parse
    gff, syn_fst, syn_fst_ins, prefix = get_parser()
    # INS
    syn_fst_ins_tmp_path, syn_fst_ins_type = syn_fst_ins_file(syn_fst, syn_fst_ins)
    syn_fst_ins_type = "INS"
    syn_fst_ins_dict = build_syn_fst_index(syn_fst_ins_tmp_path)
    rm_commmand = f'rm {syn_fst_ins_tmp_path}'
    run(rm_commmand, env_path)
    # result
    if syn_fst_ins_type == "INS":
        window = 5000
        gff_dict = indexer(gff, window)
        out_syn_fst_gene_file = f'{prefix}.gene_syn_Fst.out.ins'
        target_gene(gff_dict, syn_fst_ins_dict, out_syn_fst_gene_file, window)
        sed_commmand = f'sed \'s/(ALL)//g\' {out_syn_fst_gene_file} | awk \'NR==1 || /INS/\' | sed \'1 i #Insertion\' ' \
                       f'> {out_syn_fst_gene_file}.proces && mv {out_syn_fst_gene_file}.proces {out_syn_fst_gene_file}'
        run(sed_commmand, env_path)
    # ALL
    syn_fst_tmp_path, syn_fst_type = syn_fst_file(syn_fst)
    syn_fst_type = "ALL"
    syn_fst_dict = build_syn_fst_index(syn_fst_tmp_path)
    rm_commmand = f'rm {syn_fst_tmp_path}'
    run(rm_commmand, env_path)
    if syn_fst_type == "ALL":
        window = 2000
        gff_dict = indexer(gff, window)
        out_syn_fst_gene_file = f'{prefix}.gene_syn_Fst.out'
        target_gene(gff_dict, syn_fst_dict, out_syn_fst_gene_file, window)
        sed_commmand = f'sed -i \'s/(ALL)//g\' {out_syn_fst_gene_file}'
        run(sed_commmand, env_path)

if __name__ == '__main__':
    main()