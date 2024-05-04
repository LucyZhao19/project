# citation: https://snakemake.readthedocs.io/en/stable/tutorial/basics.html

import pandas

# use pandas to refer to user input (stored in 'user_input.csv')
DATA = pandas.read_csv('user_input.csv')
# store user sample selection
SAMPLE = DATA['sample'][0]
# store user reference selection
REFERENCE = DATA['reference'][0]

# before executing the snakemake workflow, delete all previously generated snakemake files
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#onstart-onsuccess-and-onerror-handlers
onstart:
    shell("snakemake --delete-all-output")

# execute all following rules
rule all:
    input:
        expand("static/{sample}_{reference}.txt", sample=SAMPLE, reference=REFERENCE)

# comma prevents concatenation of the two lines. 
rule bwa_map:
    input:
        fa=expand("hg38/{reference}.fa", reference=REFERENCE),
        fastq=expand("fastq/{sample}.fastq", sample=SAMPLE)
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input.fa} {input.fastq} | samtools view -Sb - > {output}"

# If the command is written across multiple lines, 
# have one white space at the end of each line except the last line. 
# python automatically concatenates the lines together
rule samtools_sort:
    input: 
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

# index read alignments
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

# generate the vcf file (which contains a list of called mutations) via beftools mpileup and call
rule bcftools_call:
    input:
        fa="hg38/{reference}.fa",
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai"
    output:
        "static/{sample}_{reference}.txt"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

# automatically generate snakemake report if snakefile was successfully executed
# https://stackoverflow.com/questions/58200594/is-it-possible-in-snakemake-to-produce-reports-and-the-dag-images-automatically
onsuccess:
    shell("snakemake --report static/report.html")

# # give summary statistics
# rule plot_quals:
#     input:
#         "static/{sample}_{reference}.txt"
#     output:
#         "static/{sample}_{reference}_quals.png"
#     script:
#         "plot-quals.py"
