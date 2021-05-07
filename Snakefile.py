configfile: "config_H.yaml"

rule all:
    input:
        expand("{species}/Counts/{species}_counts.txt", species=config["species"]),

rule flexbar:
    input:
        R1="{species}/raw/{sample}_1.fq.gz",
        R2="{species}/raw/{sample}_2.fq.gz",
        adapters="illumina_multiplex.fa"
    output:
        T1="{species}/trim/{species}_{sample}_trimmed_1.fastq.gz",
        T2="{species}/trim/{species}_{sample}_trimmed_2.fastq.gz"
    conda: "envs/TRIM.yaml"
    threads: 4
    shell:
        """flexbar -r {input.R1} -p {input.R2} \
        -a {input.adapters} --adapter-trim-end ANY \
        -q TAIL -qf i1.5 -n {threads} -m 25 \
        -t {wildcards.species}/trim/{wildcards.species}_{wildcards.sample}_trimmed \
        -z GZ"""

rule map:
    input:
        T1="{species}/trim/{species}_{sample}_trimmed_1.fastq.gz",
        T2="{species}/trim/{species}_{sample}_trimmed_2.fastq.gz",
    params:
        Genome="UCSC/genome"
    output:
        BAM="{species}/bam/{species}_{sample}.bam"
    conda: "envs/MAPPING.yaml"
    threads: 8
    shell:
        """hisat2 -p {threads} -x {params.Genome} -1 {input.T1} -2 {input.T2}| \
        samtools view -b -F 0x4 > {output.BAM}"""

rule counts:
    input:
        BAM= expand("{species}/bam/{species}_{sample}.bam", species=config["species"], sample=config["sample"]),
        annotation="gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz",
        alias="GENCODE_chromosome_aliases.txt",
    conda: "envs/COUNTS.yaml"
    output:
        counts="{species}/Counts/{species}_counts.txt"
    shell:
        "featureCounts -p -t 'exon' -g 'gene_name' -a {input.annotation} -A {input.alias} -o {output.counts} {input.BAM}"
