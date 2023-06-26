# SynDiv

[![GitHub Downloads](https://img.shields.io/github/downloads/JiaoLab2021/SynDiv/total.svg?style=social&logo=github&label=Download)](https://github.com/JiaoLab2021/SynDiv/releases)
[![BioConda Install](https://img.shields.io/conda/dn/duzezhen/syndiv.svg?style=flag&label=BioConda%20install)](https://anaconda.org/DuZeZhen/syndiv)
[![GitHub last commit](https://img.shields.io/github/last-commit/JiaoLab2021/syndiv.svg?label=Last%20commit&logo=github&style=flat)](https://github.com/JiaoLab2021/SynDiv/releases)
[![Build Status](https://github.com/JiaoLab2021/SynDiv/actions/workflows/ci.yaml/badge.svg)](https://github.com/JiaoLab2021/SynDiv/actions)

## Introduction

A tool for quick and accurate calculation of syntenic diversity.

## Requirements

[VG_url]: https://github.com/vgteam/vg
[GraphAligner_url]: https://github.com/maickrau/GraphAligner
[Paragraph_url]: https://github.com/Illumina/paragraph
[BayesTyper_url]: https://github.com/bioinformatics-centre/BayesTyper
[GraphTyper2_url]: https://github.com/DecodeGenetics/graphtyper
[PanGenie_url]: https://github.com/eblerjana/pangenie

Please note the following requirements before building and running the software:

* `Linux` operating system
* cmake version `3.12` or higher
* Python version `3.8` or higher
* C++ compiler that supports `C++14` or higher, and the `zlib` library installed (we recommend using GCC version `"4.9"` or newer) for building `graphvcf` and `fastAQ`
* The following dependencies must also be installed: [VG][VG_url], [GraphAligner][GraphAligner_url], [Paragraph][Paragraph_url], [BayesTyper][BayesTyper_url], [GraphTyper2][GraphTyper2_url], [PanGenie][PanGenie_url]

## Installation

**Building on Linux**

Use the following script to build the software:

1. First, obtain the source code.

```shell
git clone https://github.com/JiaoLab2021/EVG.git
cd EVG
```

2. Next, compile the software and add the current directory to your system's `PATH` environment variable. Please make sure that `EVG`, `graphvcf`, and `fastAQ` are all in the same folder, as `EVG` will call these two programs from its own directory.

```shell
cmake ./
make
chmod +x EVG.py
ln -sf EVG.py EVG
echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
source ~/.bashrc
```

3. Assuming that you have installed all the required software dependencies, please make sure they have been added to your environment path or activated in the corresponding `code` environment. If you haven't installed them yet, you can use the following `code` to install all the dependencies:

```shell
# To create a conda environment named graph (you can replace it with any other name), 
# make sure to replace all occurrences of graph in the following code with the name you 
# have chosen.
conda create -n graph
conda activate graph
# Install software using conda
conda install vg graphaligner paragraph bayestyper graphtyper2 pangenie
# To install PanGenie using conda under the same environment, replace pangenie in 
# the environment.yml file with the name you have chosen.
git clone https://github.com/eblerjana/pangenie.git
cd pangenie
sed  -i 's/pangenie/graph/' environment.yml
conda env update --file environment.yml
echo 'export PATH="$PATH:'$(pwd)/src'"' >> ~/.bashrc
# Install kmc using conda
conda install kmc
```

4. To verify that the software has been installed correctly, perform a test run using the following steps:

```shell
EVG -h
graphvcf -h
fastAQ -h
vg -h
GraphAligner -h
paragraph -h
bayesTyper -h
graphtyper -h
PanGenie -h
# test
cd test
EVG -r test.fa -v test.vcf.gz -s sample.txt --software VG-MAP VG-Giraffe GraphAligner Paragraph BayesTyper GraphTyper2 PanGenie &>log.txt &
```

## Usage

**Input Files**

* Reference Genome
* VCF File of Population Variants
* Sample File:

```shell
# Sample File
sample1 path_to_sample1_read1 path_to_sample1_read2
sample2 path_to_sample2_read1 path_to_sample2_read2
...
sampleN path_to_sampleN_read1 path_to_sampleN_read2
```

Please note that the Sample file must be formatted exactly as shown above, where each sample is listed with its corresponding read files.

**Running**

For convenience, let's assume the following file names for the input:

* `refgenome.fa`
* `input.vcf.gz`
* `sample.txt`

`EVG` automatically selects suitable software based on the genome, mutation and sequencing data. If desired, users can also use the `"--software"` command to specify their preferred software. The default running command is as follows:

```shell
EVG -r refgenome.fa -v input.vcf.gz -s sample.txt
```

The results are stored in the `merge/` folder, and each file is named after the corresponding sample listed in `sample.txt`: `sample1.vcf.gz`, `sample2.vcf.gz`, ..., `sampleN.vcf.gz`.

```shell
$ tree merge/
merge/
├── test1.vcf.gz
└── test2.vcf.gz

0 directories, 2 files
```

**Parameter**

* `--depth`: This parameter specifies the maximum sequencing data depth allowed for downstream analysis. If this value is exceeded, EVG will randomly downsample reads to the specified level in order to speed up the run. The default downsampling level is set at 15×, but it can be adjusted to meet specific requirements.
* `--mode`: This parameter determines the operating mode of `EVG`. In fast mode, only certain software is utilized to genotype SNPs and indels, while precise mode employs all software to genotype all variants.
* `--force`: If there are pre-existing files in the running directory of `EVG`, this parameter can be used to forcibly empty the folder. Otherwise, the software will encounter an error and exit.
* `--restart`: This parameter allows the software to resume from where it left off if it unexpectedly stops, enabling a breakpoint restart. Note that software completion is determined by file existence. It's recommended to manually check for incomplete or empty files before using this parameter and delete them.

## Citation

Please cite:

[article][article_url]

## License

MIT