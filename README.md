[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/olivar/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/olivar/badges/downloads.svg)](https://anaconda.org/bioconda/olivar)

# Olivar multiplex PCR tiling design

## Description
Olivar is a Python3 software for multiplex PCR tiling design. Olivar implements a novel algorithm to reduce non-specific amplifications in PCR, while avoiding primer design at SNPs and other undesired regions at the same time. Olivar also optimize for primer dimers with the [SADDLE](https://www.nature.com/articles/s41467-022-29500-4) algorithm. Olivar is also published as a preprint on [bioRxiv](https://doi.org/10.1101/2023.02.11.528155). 
![](Figures/Fig1.png)

## Web Interface

A web interface is available at [olivar.rice.edu](https://olivar.rice.edu/), although it does not support all available functions at the moment. 

## Install with conda (Linux or macOS Intel)

#### 1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if not installed already
e.g., download latest miniconda installer for macOS Intel x86, and in terminal
```
bash Miniconda3-latest-MacOSX-x86_64.sh
```

#### 2. Create a new conda environment and install olivar via [Bioconda](https://bioconda.github.io/)
```
conda create -n olivar olivar --channel conda-forge --channel bioconda --channel defaults --strict-channel-priority
```
Setting channel priority is important for [blast](https://bioconda.github.io/recipes/blast/README.html) to function properly. The Bioconda version of blast does not support Windows or Apple silicon at the moment. 

## Dependencies
```
python >=3.8
blast >=2.12.0
biopython
numpy
pandas
plotly >=5.13.0
tqdm
```
Use BLAST 2.12.0 or 2.13.0 to reproduce the results in example_output. 


## Usage

### Input files

 - (Required) Reference sequence in fasta format ([example](example_input/EPI_ISL_402124.fasta)). **Ambiguous bases are currently not supported and will raise errors.**

 - (Optional) List of sequence variations to be avoided, in csv format ([example](example_input/delta_omicron_loc.csv)). Column "START" and "STOP" are required, "FREQ" is considered as 1.0 if empty. Other columns are not required. Coordinates are 1-based. 

 - (Optional) A BLAST database of non-specific sequences ([example](example_input/Human)). To make your own BLAST database, check out the [NCBI BLAST User Manual](https://www.ncbi.nlm.nih.gov/books/NBK569841/). NCBI BLAST Command Line Applications are already installed along with Olivar.\
The example BLAST database is created with 23 Chromosomes and MT of human genome assembly [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/), with the command below
```
makeblastdb -in GRCh38_primary.fasta -dbtype nucl -title GRCh38_primary -parse_seqids -hash_index -out GRCh38_primary -max_file_sz 4GB -logfile makeblastdb.out -taxid 9606
```
`makeblastdb` is installed along with Olivar. 

### Command-line interface

The Olivar CLI tool comprises of four sub-commands: `build`, `tiling`, `save` and `validate`. Descriptions of command-line arguments can be found in [Command-line parameters](#command-line-parameters). 

#### 1. Build Olivar reference
A fasta reference sequence is required, coordinates of sequence variations and BLAST database are optional. 
```
olivar build example_input/EPI_ISL_402124.fasta -v example_input/delta_omicron_loc.csv -d example_input/Human/GRCh38_primary -o example_output -p 1
```
An Olivar reference file ([olivar-ref.olvr](example_output/olivar-ref.olvr)) will be generated. Use multiple CPU cores (`-p`) to accelerate this process. 

In this step, the input reference sequence is chopped into kmers, and GC content, sequence complexity and BLAST hits are calculated for each kmer. Sequence variations are also labeled if coordinates are provided. A risk score is assigned to each nucleotide of the reference sequence, guiding the placement of primer design regions. 

#### 2. Design tiled amplicons
An Olivar reference file generated in step 1 is required. Set random seed (`--seed`) to make the results reproducible. Use multiple CPU cores (`-p`) to accelerate this process. Output files are listed below (coordinates are 1-based). 
```
olivar tiling example_output/olivar-ref.olvr -o example_output --max-amp-len 420 --min-amp-len 252 --check-var --seed 10 -p 1
```
| Default name &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Description|
| :-------  | :-------- | 
| olivar-design.olvd| Olivar design file, keeping all intermediate results during the design. |
| olivar-design.csv| Sequences, coordinates (1-based) and pool assignment of primers, inserts and amplicons. |
| olivar-design.json| Design configurations. |
| olivar-design.fasta| Reference sequence. |
| olivar-design.html| An interactive plot to view primers and the risk array. |
| olivar-design_Loss.html| Loss of PDR optimization and primer dimer optimization. |
| olivar-design_risk.csv| Risk scores of each risk component. |
| olivar-design.scheme.bed| Primer sequences and coordinates in [ARTIC/PrimalScheme](https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019) format. |

In this step, the placement of primer design regions (PDRs) is optimized based on the risk array ([Fig.1d](Figures/Fig1.png)), and primer candidates are generated by SADDLE for each PDR in the optimized PDR set. SADDLE also minimizes primer dimer by exploring different combinations of primer candidates. 

#### (Optional) Load from a previous Olivar design and save output files
Output files in step 2 can be generated repeatedly as long as the Olivar deisng file (.olvd) is provided. 
```
olivar save example_output/olivar-design.olvd -o example_output
```

#### (Optional) Validate existing primer pools
Input should be a csv file, with four required columns: "amplicon_id" (amplicon name), "fP" (sequence of forward primer), "rP" (sequence of reverse primer) and "pool" (primer pool number, e.g., 1). This could be an Olivar designed primer pool generated in step 2, or primer pools that are not designed by Olivar. Output files are listed below (coordinates are 1-based). Use multiple CPU cores (`-p`) to accelerate this process. 
```
olivar validate example_output/olivar-design.csv --pool 1 -d example_input/Human/GRCh38_primary -o example_output -p 1
```
| Default name &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Description|
| :-------  | :-------- | 
| olivar-val_pool-1.csv| Basic information of each single primer, including dG, dimer score, BLAST hits, etc. |
| olivar-val_pool-1_ns-amp.csv| Predicted non-specific amplicons. |
| olivar-val_pool-1_ns-pair.csv| Predicted non-specific primer pairs. |

### Import Olivar as a Python package

Olivar can also be imported as a Python package, comprising of four functions with the same names and parameters as the four sub-commands in the CLI. 
```
from olivar import build, tiling, save, validate
```
Refer to [example.py](example.py) for more details.


## Command-line parameters

#### sub-command: `build`
```
olivar build fasta-file [--var <string>] [--db <string>] [--output <string>] 
[--title <string>] [--threads <int>]
```
| Argument &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Default &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Description|
| :-------  | :----- | :-------- | 
| fasta-file| | Positional argument. Path to the fasta reference sequence.|
|--var, -v| **None**| Optional, path to the csv file of SNP coordinates and frequencies. Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.|
|--db, -d| **None**| Optional, path to the BLAST database. Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").|
|--output, -o| ./| Output directory (output to current directory by default).|
|--title, -t| olivar-ref| Name of the Olivar reference file.|
|--threads, -p| 1| Number of threads.|

#### sub-command: `tiling`
```
olivar tiling olvr-file [--output <string>] [--title <string>] [--max-amp-len <int>] 
[--min-amp-len <int>] [--w-egc <float>] [--w-lc <float>] [--w-ns <float>] [--w-var <float>] 
[--temperature <float>] [--salinity <float>] [--dg-max <float>] [--min-gc <float>] 
[--max-gc <float>] [--min-complexity <float>] [--max-len <int>] [--check-var] 
[--fp-prefix <DNA>] [--rp-prefix <DNA>] [--seed <int>] [--threads <int>]
```
| Argument &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Default &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Description|
| :-------  | :----- | :-------- | 
| olvr-file| | Positional argument. Path to the Olivar reference file (.olvr).|
|--output, -o| ./| Output path (output to current directory by default).|
|--title, -t| olivar-design| Name of design.|
|--max-amp-len| 420| Maximum amplicon length.|
|--min-amp-len| None| Minimum amplicon length. 0.9*{max-amp-len} if not provided.|
|--w-egc| 1.0| Weight for extreme GC content.|
|--w-lc| 1.0| Weight for low sequence complexity.|
|--w-ns| 1.0| Weight for non-specificity.|
|--w-var| 1.0| Weight for variations.|
|--temperature| 60.0| PCR annealing temperature.|
|--salinity| 0.18| Concentration of monovalent ions in units of molar.|
|--dg-max| -11.8| Maximum free energy change of a primer in kcal/mol.|
|--min-gc| 0.2| Minimum GC content of a primer.|
|--max-gc| 0.75| Maximum GC content of a primer.|
|--min-complexity| 0.4| Minimum sequence complexity of a primer.|
|--max-len| 36| Maximum length of a primer.|
|--check-var| False| Boolean flag. Filter out primer candidates with variations within 5nt of 3' end. NOT recommended when a lot of variations are provided, since this would significantly reduce the number of primer candidates. |
|--fp-prefix| None| Prefix of forward primer. Empty by default.|
|--rp-prefix| None| Prefix of reverse primer. Empty by default.|
|--seed| 10| Random seed for optimizing PDRs and SADDLE.|
|--threads, -p| 1| Number of threads.|

#### sub-command: `save`
```
olivar save olvd-file [--output <string>]
```
| Argument &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Default| Description|
| :-------  | :----- | :-------- | 
| olvd-file| | Positional argument. Path to the Olivar design file (.olvd)|
|--output, -o| ./| Output directory (output to current directory by default).|

#### sub-command: `validate`
```
olivar validate csv-file [--pool <int>] [--db <string>] [--output <string>] 
[--title <string>] [--max-amp-len <int>] [--temperature <float>] [--threads <int>]
```
| Argument &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Default &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Description|
| :-------  | :----- | :-------- | 
| csv-file| | Positional argument. Path to the csv file of a primer pool. Required columns: "amplicon_id" (amplicon name), "fP" (sequence of forward primer), "rP" (sequence of reverse primer), "pool" (pool number, e.g., 1).|
|--pool| 1| Primer pool number. |
|--db, -d| None| Optional, path to the BLAST database. Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").|
|--output, -o| ./| Output directory (output to current directory by default).|
|--title, -t| olivar-val| Name of validation.|
|--max-amp-len| 1500| Maximum length of predicted non-specific amplicon. Ignored is no BLAST database is provided.|
|--temperature| 60.0| PCR annealing temperature.|
|--threads, -p| 1| Number of threads.|
