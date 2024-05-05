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
    # set plot and vcf to None.
    # if POST, plot will store snakemake image output
    # if POST, vcf will store snakemake vcf output
    plot = "None"
    vcf = "None"

    # if the method is POST:
    if request.method == "POST":
        # try get user selected human reference genome
        try: 
            reference = request.form.get("reference")
        # except when request.form.get gives a ValueError (then reference = None)
        except ValueError:
            reference = None
        print(reference)
        # try get user fastq file
        try: 
            fastq = request.form.get("fastq")
            # split fastq file name into sample name and extension
            sample, extension = fastq.split(".")  
        # except when request.form.get gives a ValueError (then fastq = None)
        except ValueError:
            fastq = None
        print(fastq)
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
            command = ["snakemake", "--forceall", "--cores", "1"]
            # captures the outputs of all intermediate process
            # snakemake_output = subprocess.run(command, capture_output=True, text=True)
            subprocess.run(command, capture_output=True, text=True)
            # store name of output vcf file
            expected_vcf = sample + "_" + reference + ".txt"
            print(expected_vcf)

                
            # render plot to result.html
            return render_template("result.html", snakemake_plot=plot, snakemake_vcf=expected_vcf)
        
        # otherwise, display "error: invalid fastq and/or reference selection. Choose a valid option!"
        error = "error: invalid fastq and/or reference selection. Choose a valid option!"

        # display error message in error.html
        return render_template("error.html", error=error)
    
    # otherwise, stay on index.html without showing any content under snakemake result
    return render_template("index.html")
