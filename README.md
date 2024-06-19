# SynDiv

<!-- [![GitHub Downloads](https://img.shields.io/github/downloads/JiaoLab2021/SynDiv/total.svg?style=social&logo=github&label=Download)](https://github.com/JiaoLab2021/SynDiv/releases) -->
[![BioConda Install](https://img.shields.io/conda/dn/duzezhen/syndiv.svg?style=flag&label=BioConda%20install)](https://anaconda.org/DuZeZhen/syndiv)
[![GitHub last commit](https://img.shields.io/github/last-commit/JiaoLab2021/syndiv.svg?label=Last%20commit&logo=github&style=flat)](https://github.com/JiaoLab2021/SynDiv/releases)
[![Build Status](https://github.com/JiaoLab2021/SynDiv/actions/workflows/ci.yaml/badge.svg)](https://github.com/JiaoLab2021/SynDiv/actions)

## Introduction

A tool for quick and accurate calculation of syntenic diversity.

## Requirements

Please note the following requirements before building and running the software:

* `Linux` operating system
* cmake version `3.12` or higher
* Python version `3.6` or higher and the following packages: numpy-1.21.2, pandas-1.2.4, pysam-0.16.0.1 and psutil-4.4.2
* C++ compiler that supports `C++17` or higher, and the `zlib` library installed (we recommend using GCC version `"7.3.0"` or newer) for building `SynDic_c`
* `samtools ≥ 1.15`, `minimap2 ≥ 2.20`

## Installation

**Install via conda**

```shell
conda create -n syndiv
conda activate syndiv
# Install SynDiv with all dependencies
conda install -c bioconda -c conda-forge -c duzezhen syndiv
```

**Building on Linux**

Use the following script to build the software:

1. First, obtain the source code.

```shell
git clone https://github.com/JiaoLab2021/SynDiv.git
cd SynDiv
```

2. Next, compile the software and add the current directory to your system's `PATH` environment variable.

```shell
cmake ./
make
chmod +x SynDiv.py SynDiv_p.py genome2SynDiv_config.py
ln -sf SynDiv.py SynDiv
ln -sf SynDiv_p.py SynDiv_p
ln -sf genome2SynDiv_config.py genome2SynDiv_config
echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
source ~/.bashrc
```

3. Assuming that you have installed all the required software dependencies, please make sure they have been added to your environment path or activated in the corresponding `code` environment. If you haven't installed them yet, you can use the following `code` to install all the dependencies:

```shell
conda create -n syndiv
conda activate syndiv
# Install software using conda
conda install samtools minimap2
```

4. To verify that the software has been installed correctly, perform a test run using the following steps:

```shell
SynDiv -h
SynDiv_p -h
SynDiv_c -h
samtools
minimap2
# test
cd test
nohup /usr/bin/time -v SynDiv -r genome/refgenome.fa -c configuration.txt &>log.txt &
```

## Usage

**Input Files**

To quickly get started, you will need two input files: `aligns` files and `syri.out` files. Once you have obtained these files, make sure to prepare the Reference genome and configuration file.

* Reference Genome
* configuration file

Please note that the chromosome names in the query genome must match those in the reference genome.

```shell
# configuration file
sample1 sample1.fa sample1.aligns sample1.syri.out
sample2 sample2.fa sample2.aligns sample2.syri.out
...
sampleN sampleN.fa sampleN.aligns sampleN.syri.out
```

[configuration_article]: https://github.com/JiaoLab2021/SynDiv/wiki/Configuration-file

File should be separated by tabs. The code examples for generating `aligns` and `syri.out` files can be found on the [wiki][configuration_article].

**Running**

Before running the software, it is recommended to set the maximum number of open files using the `ulimit -n <number>` command. The maximum number of open files can be calculated based on the number of genomes (`n`) and the number of threads (`t`) using the following formula:

```shell
number = 10 + t*(2n - t - 1)
```

For convenience, let's assume the following file names for the input:

* `refgenome.fa`
* `configuration.txt`

```shell
ulimit -n 50000
SynDiv -r refgenome.fa -c configuration.txt &
```

## Citation

Please cite:

[article][article_url]

## License

MIT