in_dir = config['input']
out_dir = config['output']
samples = config['samples'].split()
R=['1', '2']

rule all:
    input:
        out_dir + "/multiqc_raw.html",
        expand(in_dir + "/{sample}_R1_trimmed.fastq.gz", sample=samples),
        expand(in_dir + "/{sample}_R2_trimmed.fastq.gz", sample=samples),
        out_dir + "/multiqc_trimmed.html",
        out_dir + "/reference_index",
        expand(out_dir + "/STAR/{sample}_aligned.bam", sample=samples),
        expand(out_dir + "/STAR/{sample}_aligned_sorted.bam", sample=samples),
        expand(out_dir + "/FEATURECOUNTS/{sample}_counts.txt", sample=samples)

rule fastqc:
    input:
        in_dir + "/{sample}.fastq.gz"
    output:
        out_dir + "/fastqc/{sample}_fastqc.zip",
        out_dir + "/fastqc/{sample}_fastqc.html"
    conda:
        "envs/environment.yaml"
    params:
        out_d = out_dir + "/fastqc/"
    shell:
        "fastqc {input} -o {params.out_d}"

rule multiqc:
    input:
        expand(out_dir + "/fastqc/{sample}_L001_R{R}_001_fastqc.zip", sample = samples, R=R)
    output:
        out_dir + "/multiqc_raw.html"
    conda:
        "envs/environment.yaml"
    shell:
        "multiqc {input} -o " + out_dir + " -n multiqc_raw"

rule multiqc_trimmed:
    input:
        expand(out_dir + "/fastqc/{sample}_R{R}_trimmed_fastqc.zip", sample = samples, R=R)
    output:
        out_dir + "/multiqc_trimmed.html"
    conda:
        "envs/environment.yaml"
    shell:
        "multiqc {input} -o " + out_dir + " -n multiqc_trimmed"

rule bbduk_trim:
    input:
        in_dir + "/{sample}_L001_R1_001.fastq.gz",
        in_dir + "/{sample}_L001_R2_001.fastq.gz"
    output:
        in_dir + "/{sample}_R1_trimmed.fastq.gz",
        in_dir + "/{sample}_R2_trimmed.fastq.gz"
    conda:
        "envs/environment.yaml"
    shell:
        "bbduk.sh ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 " +
        "in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]}"


rule star_index:
    input:
        reference_dna = "ref_annot/chr19_20Mb.fa",
        reference_gtf = "ref_annot/chr19_20Mb.gtf"
    output:
        out_dir + "/reference_index"
    shell:
        "mkdir {output} && " +
        "STAR --runMode genomeGenerate " +
        "--runThreadN 4 " +
        "--genomeDir {output} " +
        "--genomeFastaFiles {input.reference_dna} " +
        "--sjdbGTFfile {input.reference_gtf}"

rule run_star:
    input:
        in_dir + "/{sample}_R1_trimmed.fastq.gz",
        in_dir + "/{sample}_R2_trimmed.fastq.gz",
        genomeDir = out_dir + "/reference_index"
    output:
        out_dir + "/STAR/{sample}_aligned.bam"
    shell:
        "STAR --genomeDir {input.genomeDir} " +
        "--readFilesCommand zcat " +
        "--runThreadN 4 " +
        "--outFileNamePrefix {output} " +
        "--runMode alignReads " +
        "--readFilesIn {input[0]} {input[1]} " +
        "--outSAMtype BAM SortedByCoordinate && " +
        "mv {output}Aligned.sortedByCoord.out.bam {output}"

rule sort_bam:
    input:
        out_dir + "/STAR/{sample}_aligned.bam"
    output:
        out_dir + "/STAR/{sample}_aligned_sorted.bam"
    shell:
        "samtools sort -n {input} -o {output}"

rule index_bam:
    input:
        out_dir + "/STAR/{sample}_aligned_sorted.bam"
    output:
        out_dir + "/STAR/{sample}_aligned_sorted.bam.bai"
    shell:
        "samtools index {input}"

rule featurecounts:
    input:
        out_dir + "/STAR/{sample}_aligned.bam",
        reference_gtf = "ref_annot/chr19_20Mb.gtf"
    output:
        out_dir + "/FEATURECOUNTS/{sample}_counts.txt",
    shell:
        "featureCounts -p -t exon -g gene_id -a {input.reference_gtf} " +
        "-o {output} {input[0]} -s 1"
