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
