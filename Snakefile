# citation: https://snakemake.readthedocs.io/en/stable/tutorial/basics.html
# citation: https://eriqande.github.io/eca-bioinf-handbook/managing-workflows-with-snakemake.html
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
# https://bio-bwa.sourceforge.net/bwa.shtml
# https://www.htslib.org/doc/samtools-sort.html
rule align_and_sort_reads:
    # specifies input reference chromosome (.fa) and sample (.fastq)
    input:
        fa=expand("hg38/{reference}.fa", reference=REFERENCE),
        fastq=expand("fastq/{sample}.fastq", sample=SAMPLE)
    # outputs bam file
    output:
        "bam_files/{sample}.bam"
    # bwa mem aligns read to the reference chromosome, 
    # samtools sort re-order the reads in order from the first to the last chromosomal position and generates a BAM output (-O bam).
    shell:
        "bwa mem {input.fa} {input.fastq} | "
        "samtools sort -O bam > {output}"

# samtools index serves to index read for DNA mutation calling and generates a BAM.BAI output.
# the BAM.BAI output is required for subsequent rule/step.
# https://www.htslib.org/doc/samtools-index.html
rule index_bam:
    # inputs bam that needs to be indexed
    input:
        "bam_files/{sample}.bam"
    # outputs indexed bam file (.bam.bai)
    output:
        "bam_files/{sample}.bam.bai"
    # index bam file via samtools index
    shell:
        "samtools index {input}"

# generate the vcf file (which contains a list of DNA mutations) via beftools mpileup and call
# https://samtools.github.io/bcftools/bcftools.html#mpileup
rule call_mutation:
    # input specified .fa, .bam, and bam.bai files
    input:
        fa="hg38/{reference}.fa",
        bam="bam_files/{sample}.bam",
        bai="bam_files/{sample}.bam.bai"
    # outputs VCF file (.txt)
    output:
        "static/{sample}_{reference}.txt"
    # in bcftools mpileup, specifies the reference (--fasta-ref) and input bam file. 
    # in bcftools call, use the multiallelic calling model (-m) to identify mutations in samples 
    # and outputs only mutations (-v) in the vcf file. 
    shell:
        "bcftools mpileup --fasta-ref {input.fa} {input.bam} | "
        "bcftools call -m -v > {output}"

# automatically generate snakemake report if snakefile was successfully executed
# https://stackoverflow.com/questions/58200594/is-it-possible-in-snakemake-to-produce-reports-and-the-dag-images-automatically
onsuccess:
    shell("snakemake --report static/report.html")
