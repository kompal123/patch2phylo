import os
import pandas as pd

# Load sample table 
samples_df = pd.read_csv(config["samples_tsv"], sep=r"\s+", header=0)
fq1_dict = dict(zip(samples_df["sample"], samples_df["fq1"]))
fq2_dict = dict(zip(samples_df["sample"], samples_df["fq2"]))

rule fastp_trim:
    threads: 4
    input:
        fq1 = lambda wc: fq1_dict[wc.sample],
        fq2 = lambda wc: fq2_dict[wc.sample]
    output:
        trimmed1 = os.path.join(config["results"]["qc_trimmed"], "{sample}_1.trimmed.fastq.gz"),
        trimmed2 = os.path.join(config["results"]["qc_trimmed"], "{sample}_2.trimmed.fastq.gz"),
        html     = os.path.join(config["logs"]["fastp"], "{sample}_fastp.html"),
        json     = os.path.join(config["logs"]["fastp"], "{sample}_fastp.json")
    params:
        adapter = config["adapter_path"],
        quality = config["parameters"]["fastp"]["qualified_quality_phred"],
        unqual  = config["parameters"]["fastp"]["unqualified_percent_limit"],
        nlimit  = config["parameters"]["fastp"]["n_base_limit"],
        length  = config["parameters"]["fastp"]["length_required"],
        detect  = "--detect_adapter_for_pe" if config["parameters"]["fastp"]["detect_adapter_for_pe"] else ""
    conda:
        "../envs/fastp.yaml"
    log:
        os.path.join(config["logs"]["fastp"], "{sample}.log")
    shell:
        """
        fastp \
            --in1 {input.fq1} \
            --in2 {input.fq2} \
            --out1 {output.trimmed1} \
            --out2 {output.trimmed2} \
            --adapter_fasta {params.adapter} \
            --qualified_quality_phred {params.quality} \
            --unqualified_percent_limit {params.unqual} \
            --n_base_limit {params.nlimit} \
            --length_required {params.length} \
            {params.detect} \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            > {log} 2>&1
        """
