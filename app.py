# import os: https://docs.python.org/3/library/os.html
import os
# load flask
from flask import Flask, render_template, redirect, request
# load functions needed for uploading files
# https://flask.palletsprojects.com/en/3.0.x/patterns/fileuploads/
from werkzeug.utils import secure_filename
# import subprocess for running shell commands
import subprocess
# import csv module
import csv

# application
app = Flask(__name__)

# # set folder path of users' uploaded files
# FASTQ_FOLDER = '/workspaces/133161714/final_project/fastq'

# define a list of extension for fastq files
EXTENSIONS = ["fq", "fastq"]

# define a list of possible chromosome reference (chr1 to chr22)
REFERENCE = ["chr" + str(i) for i in range(1, 23)]

# # configure path to folder containing the user uploaded files.
# app.config["FASTQ_FOLDER"] = FASTQ_FOLDER

# configure app route
@app.route("/", methods=['GET', 'POST'])
def index():
    # set vcf and report to None.
    vcf = "None"
    report = "None"

    # # set error message "error: invalid fastq and/or reference selection. Choose a valid option!"
    # error = "error: invalid fastq and/or reference selection. Choose a valid option!"
    
    # if the method is POST:
    if request.method == "POST":
        
        # try get user fastq file
        try: 
            fastq = request.form.get("fastq")
        # except when request.form.get gives a ValueError
        except ValueError:
            # display error.html
            return render_template("error.html")
    
        # try split fastq file name into sample name and extension
        try: 
            sample, extension = fastq.split(".")  
        # except when command gives a an AttributeError
        except AttributeError:
            # display error.html
            return render_template("error.html")
       
        # try get user selected human reference genome
        try: 
            reference = request.form.get("reference")
        # except when request.form.get gives a ValueError
        except ValueError or None:
            # display error.html
            return render_template("error.html")
        
        # if reference is empty, display error.html
        if reference is None:
            return render_template("error.html")
        
        # if file is a fastq file and if reference is in REFERENCE:
        if extension in EXTENSIONS and reference in REFERENCE:   
            # delete all previously generated snakemake output under the static folder
            # find output files of snakemake workflow output
            outputs = os.listdir('./static')
            # iterate over each file in outputs:
            for file in outputs:
                # if the file ends with ".html", remove the file
                if file.endswith(".html"):
                    os.remove("./static/" + file)
                # if the file ends with ".txt", remove the file
                if file.endswith(".txt"):
                    os.remove("./static/" + file)
            # delete all files in ./sorted_reads
            outputs = os.listdir("./sorted_reads")
            for file in outputs:
                os.remove("./sorted_reads/" + file)

            # write user input into a csv file (citation: https://www.geeksforgeeks.org/writing-csv-files-in-python/)
            # save user-selected sample and reference input as a dictionary
            data = [{'sample': sample, 'reference': reference}]

            # set col headers as sample and reference
            colnames = ['sample', 'reference']

            # set csv filename
            filename = "user_input.csv"

            # write to csv file
            with open(filename, 'w') as file:
                # create a dict write object
                writer = csv.DictWriter(file, fieldnames=colnames)

                #write headers
                writer.writeheader()

                #writing data rows
                writer.writerows(data)
            # perform variant calling pipeline using user input and snakemake.
            # force snakemake to re-run the entire pipeline
            command = ["snakemake", "--forceall", "--cores", "3"]
            # captures the outputs of all intermediate process
            subprocess.run(command, capture_output=True, text=True)
            # store name of output vcf file
            vcf = sample + "_" + reference + ".txt"
            # store name of expected report
            report = "report.html"
                
            # render snakemake's vcf output to result.html
            return render_template("result.html", snakemake_vcf=vcf, snakemake_report=report)
        
        # for all other unforeseen cases, display error.html
        return render_template("error.html")
    
    # otherwise, stay on index.html without showing any content under snakemake result
    return render_template("index.html")
