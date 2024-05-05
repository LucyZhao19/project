# CS50_final_project_v2
Welcome! The project contains a webpage that allows users to identify DNA mutations present in a specific chromosome (excluding sex chromosomes) of a given sample. 
note: The genomic data (fastq files) of three provided samples are taken from snakemake tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-3-creating-an-environment-with-the-required-software

Please note that the project requires conda, which can be installed from here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Download the following folders from this google drive link: https://drive.google.com/drive/folders/1Ca8MNlBiNXYzEubC6i33GsYCerNvKAHN?usp=drive_link

The link includes a "hg38.zip" folder, which contains the human reference (hg38) chromosome data. 
Note: the hg38 human reference chromosomes are downloaded from: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/. I than used samtools faidx commands to generate the .fa.fai files and used bwa index to generate the .fa.amb, .fa.ann, .fa.bwt, .fa.pac, and .fa.sa files. 

Do the following in a command line interface (bash):
Clone CS50_final_project git repo:
$ git clone https://github.com/LucyZhao19/CS50_final_project.git your_directory

Set "CS50_final_project" as the current working directory:
$ cd your_path/CS50_final_project

Place the downloaded folder in your "CS50_final_project" directory. 

Then, unzip hg38.zip:

`$ unzip hg38.zip`

Check if the unzipped folder exists in the directory:

`$ ls`

A "hg38" folder should appear in the directory. If the unzip is successful, the user can remove hg38.zip:

`$ rm hg38.zip`

Create the environment "snakemakeCS50" through conda:

`$ conda init`

`$ conda activate base`

`$ conda env create -p ./envs/snakemakeCS50 --file ./envs/environment.yaml`

Activate snakemakeCS50:

`$ conda activate ./envs/snakemakeCS50`

Run flask:

`$ flask run`

To open the webpage, click on the http link (which should appear in the "*Running on ..." line after executing the flask run command)
On the webpage, users can select any one of the samples and select a reference chromosome. Hit 

snakemakeCS50 can also be ran on the command-line interface, as instructed below:
To run snakemakeCS50, specify:

`$ snakemake --cores=1`

The human reference chromosomes can be downloaded from:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

Unzip all .fa.gz files:

`$ gunzip directory_path/*.fa.gz`

Index all .fa files:

`$ cd dir_path`

`$ for file in *.fa; do bwa index $file; done`

Generate fa.fai file from all .fa files: 

`$ cd dir_path`

`$ for file in *.fa; do samtools faidx $file; done`


TODO:
- Maybe create a return button so the user can go back to index.html?
- finalize design.md!
- make the video!

SIDE NOTE:
To update conda environment and remove any unnecessary packages from the environment:
$ conda env update -p ./envs/snakemakeCS50 --file ./envs/environment.yaml --prune
