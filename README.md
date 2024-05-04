# CS50_final_project_v2

On the command line interface, type in:
$ conda activate ./envs/snakemakeCS50

To run snakemakeCS50, specify:
$ snakemakeCS50 --cores=1 

download hg38 reference chromosome from: 
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

to unzip all .fa.gz files:
$ gunzip directory_path/*.fa.gz

To index all .fa files:
$ cd dir_path
$ for file in *.fa; do bwa index $file; done

To generate fa.fai file: 
$ cd dir_path
$ for file in *.fa; do samtools faidx $file; done