configfile: "config.yaml"

rule all:
    input:
        #expand("data/calls/{sample}.vcf", sample=config["samples"])
        #expand("data/var_freq/{sample}.csv", sample=config["samples"])
        "data/diversity.csv"

#filter sequences
rule filter_sequences:
    input:
        R1 = "data/fastq_files/{sample}_R1.fastq.gz",
        R2 = "data/fastq_files/{sample}_R2.fastq.gz"
    params:
        tm="0",
        qs="25", # the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
        ps="20" #  how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%
    log:
        "logs/fastq/{sample}.html"
    output:
        R1 = "data/fastq_filtered/{sample}_R1.fastq.gz",
        R2 = "data/fastq_filtered/{sample}_R2.fastq.gz"
    shell:
        "fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} "
        " --trim_tail1 {params.tm} --trim_tail2 {params.tm} --qualified_quality_phred {params.qs} "
        " --unqualified_percent_limit {params.ps} --html {log}"

def get_reference(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "data/reference_index/" + reference + "-reference.fa",

def get_bed(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "data/bed_file/" + reference + ".bed",

#bowtie2 mappping of filtered reads
rule bowtie2_map:
    input:
        #ref = "data/reference_index/Rp0-reference.fa" + reference,
        ref = get_reference,
        R1 = "data/fastq_filtered/{sample}_R1.fastq.gz",
        R2 = "data/fastq_filtered/{sample}_R2.fastq.gz"
    output:
        "data/mapped_reads/{sample}.bam"
    log:
        "logs/bowtie2/{sample}.log"
    threads: 8
    shell:
        "(bowtie2 -p {threads} -x {input.ref} -1 {input.R1} -2 {input.R2} "
        "-S {output}) 2> {log}"

rule samtools_sort:
    input:
        "data/mapped_reads/{sample}.bam"
    output:
        "data/sorted_reads/{sample}.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_index:
    input:
        "data/sorted_reads/{sample}.bam"
    output:
        "data/sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

def get_reference2(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return str(reference) + "-reference"

rule lofreq:
    input:
        fa=get_reference,
        bam="data/sorted_reads/{sample}.bam",
        bai="data/sorted_reads/{sample}.bam.bai",
        bed=get_bed
    params:
        p="1"
    output:
        "data/calls/{sample}.vcf"
    shell:
        #"lofreq faidx data/reference_sequence/HIV_REJOC_ENV_REF.fa "
        "lofreq call-parallel --pp-threads {params.p} -l {input.bed} -f {input.fa} -o {output} {input.bam} "

rule call_freq:
    input:
        "data/calls/{sample}.vcf",
        reference=get_reference
    output:
        "data/var_freq/{sample}.csv"
    script:
        "scripts/vcf_to_df.py"

rule diversity:
    input:
        #"data/var_freq/{sample}.csv"
        expand("data/var_freq/{sample}.csv", sample=config["samples"])
    output:
        "data/diversity.csv",
        "data/matrix.csv"
    script:
        "scripts/diversity_metrics.py"
