# CS50_final_project: DNA mutation identification (via Snakemake)
Welcome! The project contains a webpage that allows users to identify DNA mutations present in a specific chromosome (excluding sex chromosomes) of a given sample. Mutations are defined as nucleotide changes in the sample that differs from the expected nucleotide in the reference chromosome.

Youtube link to the final project's video: https://www.youtube.com/watch?v=_1EQCSyjra4

note: The genomic data (fastq files) of three provided samples are taken from snakemake tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-3-creating-an-environment-with-the-required-software

Please note that the project requires conda, which can be installed from here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Download Option 1: 
Unzip the folder submitted on Gradescope:

`unzip project.zip`

Check if the folder has been successfully unzipped (there should be a folder name called "project" without the .zip extension):

`ls`

If the unzip was successful, remove the zipped folder:

`rm project.zip`

Download Option 2: 
Download the following folders from this google drive link: https://drive.google.com/drive/folders/1Ca8MNlBiNXYzEubC6i33GsYCerNvKAHN?usp=drive_link

The link includes a "hg38.zip" folder, which contains the human reference (hg38) chromosome data. 
Note: the hg38 human reference chromosomes are downloaded from: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/. I than used samtools faidx commands to generate the .fa.fai files and used bwa index to generate the .fa.amb, .fa.ann, .fa.bwt, .fa.pac, and .fa.sa files. 

Do the following in a command line interface (bash):
Clone CS50_final_project git repo:

`git clone https://github.com/LucyZhao19/CS50_final_project.git your_directory`

Set "CS50_final_project" as the current working directory:

`cd your_path/CS50_final_project`

Place the downloaded folder in your "CS50_final_project" directory. 

Then, unzip hg38.zip:

`unzip hg38.zip`

Check if the unzipped folder exists in the directory:

`ls`

A "hg38" folder should appear in the directory. If the unzip is successful, the user can remove hg38.zip:

`rm hg38.zip`

Create the environment "snakemakeCS50" through conda:

`conda init`

`conda activate base`

`conda env create -p ./envs/snakemakeCS50 --file ./envs/environment.yaml`

Activate snakemakeCS50:

`conda activate ./envs/snakemakeCS50`

Run flask:

`flask run`

To open the webpage, click on the http link (which should appear in the "*Running on ..." line after executing the flask run command)
On the home page, users should click on the "Select fastq file" to choose a sample in the dropdown. The user should also click on the "Select human reference (hg38) chromosome" dropdown to select a reference chromosome. Please make sure to choose an option for provided in the sample and reference dropdown. Failure to do so will return an error message. Click "Submit" to submit selection. It may take a while for the content to load (underneath the hood, snakemake is working hard to identify mutations in the specified chromosome of the selected sample!).

Once the content is fully loaded, the user will see the Snakemake outputs, which include a VCF output (.txt file) and a report.html embedded in the webpage. 
- The VCF file contains a list of mutations identified in the specified chromosome of the selected sample. 
- The report.html contains details about the job execution. The "Workflow" tab displays steps of the mutation calling pipeline in sequential order. User can see the codes for each step by clicking on corresponding dot of the pipeline. The "Statistics" tab displays the runtime (in seconds) of each step and the date and time of the job execution. 

The user can download the VCF output and report.html by clicking on the "Download VCF" and "Download Report" button, respectively. 

To return to the home page, click on the "Home" button located on the top left corner of the current webpage. 

Backup plan if human reference chromosome cannot be downloaded from google drive link:

In the project folder, create a "hg38" folder:
`cd your_path/project`

`mkdir ./hg38`

`cd ./hg38`

The human reference chromosomes can be downloaded from:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

Only download the .fa.gz files with filename chr1.fa.gz, chr2.fa.gz, chr3.fa.gz, ..., chr22.fa.gz (skip the files with "GI", "KI", "chrUn", "chrX", "md5sum", and/or "chrY" in the filename)
Move the downloaded reference chromosome files and unzip all .fa.gz files in the hg38 folder:

`gunzip your_path_to_project/hg38/*.fa.gz`

Index all .fa files in the hg38 folder:

`for file in *.fa; do bwa index $file; done`

Generate fa.fai file from all .fa files in the hg38 folder: 

`for file in *.fa; do samtools faidx $file; done`
