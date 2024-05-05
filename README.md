# CS50_final_project_v2
Welcome! The project contains a webpage that allows users to identify DNA mutations present in a specific chromosome (excluding sex chromosomes) of a given sample. 

Please note that the project requires conda, which can be installed from here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Download the following folders from this google drive link: https://drive.google.com/drive/folders/1Ca8MNlBiNXYzEubC6i33GsYCerNvKAHN?usp=drive_link

The link includes: 
(1) a "fastq" folder, which contains the fastq files for sample "A", "B", and "C". In the webpage, users can select any one of these samples and identify variants specific to the sample.
note: The samples are provided by snakemake tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-3-creating-an-environment-with-the-required-software
  
(2) a "hg38" folder, which includes the human reference (hg38) chromosome data. In the webpage, users can select any one of these reference chromosomes. 
Note: the hg38 human reference chromosomes are downloaded from: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/. I than used samtools faidx commands to generate the .fa.fai files and used bwa index to generate the .fa.amb, .fa.ann, .fa.bwt, .fa.pac, and .fa.sa files. 

Do the following in a command line interface (bash):
Clone CS50_final_project git repo:
$ git clone https://github.com/LucyZhao19/CS50_final_project.git your_directory

Set "CS50_final_project" as the current working directory:
$ cd your_path/CS50_final_project

Place the downloaded folder in your "CS50_final_project" directory. 

Then, unzip the "hg38" folder:
$ tar -xzf hg38.tar.gz

Create the environment "snakemakeCS50" through conda:
$ conda init
$ conda activate base
$ conda env create -p ./envs/snakemakeCS50 --file ./envs/environment.yaml

Activate snakemakeCS50:
$ conda activate ./envs/snakemakeCS50

Run flask:
$ flask run

snakemakeCS50 can also be ran on the command-line interface, as instructed below:
To run snakemakeCS50, specify:
$ snakemakeCS50 --cores=1

The human reference chromosomes can be downloaded from:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

Unzip all .fa.gz files:
$ gunzip directory_path/*.fa.gz

Index all .fa files:
$ cd dir_path
$ for file in *.fa; do bwa index $file; done

Generate fa.fai file from all .fa files: 
$ cd dir_path
$ for file in *.fa; do samtools faidx $file; done