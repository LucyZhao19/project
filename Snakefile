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

# execute all rules required to produce the final output ({sample}_{reference}.txt)
rule all:
    input:
        expand("static/{sample}_{reference}.txt", sample=SAMPLE, reference=REFERENCE)

# reference refers to the user-selected reference chromosome (i.e. string stored in REFERNENCE)
# sample refers to the user-selected fastq filename (i.e. string stored in SAMPLE)
rule align_and_sort_reads:
    input:
        fa=expand("hg38/{reference}.fa", reference=REFERENCE),
        fastq=expand("fastq/{sample}.fastq", sample=SAMPLE)
    output:
        "sorted_reads/{sample}.bam"
    # bwa mem aligns read to the reference chromosome, 
    # samtools sort re-order the reads in order from the first to the last chromosomal position and generates a BAM output.
    shell:
        "bwa mem {input.fa} {input.fastq} | "
        "samtools sort -O bam > {output}"

# samtools index serves to index read for DNA mutation calling and generates a BAM.BAI output.
# the BAM.BAI output is required for subsequent rule/step.
rule index_bam:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

# generate the vcf file (which contains a list of DNA mutations) via beftools mpileup and call
rule call_mutation:
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
