# SuRoQ (Small RNA Quality) - a pipeline for quick and dirty QC of your small RNA (piRNA-oriented) sequencing data
SuRoQ requires only your reads in FASTQ or FASTA format (change extensions to .fq or .fa, respectively), genome assembly FASTA and TE consensus sequences FASTA. It produces three kinds of plots:
- Reads size distribution for those that mapped to genome and TEs (NOT mutually exclusive!).
- SeqLogo plots for sense and antisense TE-mapped reads, useful for validating the U-bias.
- Ping-pong signature with Z score for 10-nt overlap indicated in the title.
  
<img src=https://github.com/foriin/SuRoQ/blob/main/example/suroq_output.png width="600">

SuRoQ heavily borrows from [piPipes](https://github.com/bowhan/piPipes), namely a concepts of .insert and .BED2 files (for clarification refer to piPipes) and a couple of C++ functions that deal with those formats.

## Installation
SuRoQ works only on Linux x64 systems, it wasn't tested on Mac, but it's possible for it to work there in theory. For installation, clone this repo via
```
git clone https://github.com/foriin/SuRoQ.git
```
Then, use suroq.yml file to prepare a conda environment (I use [mamba](https://anaconda.org/conda-forge/mamba), because it is infinitely faster):
```
mamba env create -f suroq.yml
```
If you don't want to set a conda environment, here's the list of the software that you'll need:
- [bowtie](https://bowtie-bio.sourceforge.net/index.shtml)
- [samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [R](https://www.r-project.org/)
  - [ggplot2](https://ggplot2.tidyverse.org/)
  - [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
  - [ggseqlogo](https://omarwagih.github.io/ggseqlogo/)
 
## Running
Run SuRoQ with:
```
./SuRoQ.sh <your_reads> <genome.fasta> <TEs.fasta> [number_of_cores] [output_directory]
```
The last two parameters are optional but you have to specify both if you want to set only the output directory name. I will work on improving arguments handling pretty soon.
After completion, you will find the plot in the `plots` directory and all the files used for its generation in the `tables` directory.

If you have any issues with this tool, open an issue here or [email me](mailto:liartom2@gmail.com)
