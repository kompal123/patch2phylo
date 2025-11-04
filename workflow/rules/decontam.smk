
# workflow/rules/decontam.smk
# ------------------------------------------------------------
# Optional read decontamination against known contaminants
# ------------------------------------------------------------

import pandas as pd, os

# Read samples
samples_list = pd.read_csv(config["samples_tsv"], sep=r"\s+", header=0)["sample"].tolist()

# Config settings
contam_fasta = config.get("contaminants_fasta", "")
trimmed_dir  = config["results"]["qc_trimmed"]
clean_dir    = config["results"]["qc_decont"]
log_dir      = config["logs"]["decontam"]
threads      = config.get("parameters", {}).get("bowtie2", {}).get("threads", 4)

if contam_fasta:
    rule build_contam_index:
        input: fasta=contam_fasta
        output: touch(os.path.join(log_dir, "contam_index.done"))
        params: prefix=os.path.splitext(contam_fasta)[0]
        conda: "../envs/bowtie2.yaml"
        log: os.path.join(log_dir, "bowtie2_build.log")
        shell: """
            mkdir -p {log_dir}
            bowtie2-build {input.fasta} {params.prefix} > {log} 2>&1
        """

    rule filter_contaminants:
        input:
            idx_done=rules.build_contam_index.output,
            r1=lambda wc: f"{trimmed_dir}/{wc.sample}_1.trimmed.fastq.gz",
            r2=lambda wc: f"{trimmed_dir}/{wc.sample}_2.trimmed.fastq.gz"
        output:
            clean1=os.path.join(clean_dir, "{sample}_1.clean.fastq.gz"),
            clean2=os.path.join(clean_dir, "{sample}_2.clean.fastq.gz")
        params: prefix=os.path.splitext(contam_fasta)[0]
        conda: "../envs/bowtie2.yaml"
        log: os.path.join(log_dir, "{sample}.log")
        shell: """
            mkdir -p {clean_dir}
            bowtie2 --very-sensitive --threads {threads} -x {params.prefix} \
                -1 {input.r1} -2 {input.r2} --un-conc-gz {clean_dir}/{wildcards.sample}.uncleaned.fq.gz \
                > {log} 2>&1
            mv {clean_dir}/{wildcards.sample}.uncleaned.fq.1.gz {output.clean1}
            mv {clean_dir}/{wildcards.sample}.uncleaned.fq.2.gz {output.clean2}
        """
else:
    rule pass_through:
        input:
            r1=lambda wc: f"{trimmed_dir}/{wc.sample}_1.trimmed.fastq.gz",
            r2=lambda wc: f"{trimmed_dir}/{wc.sample}_2.trimmed.fastq.gz"
        output:
            clean1=os.path.join(clean_dir, "{sample}_1.clean.fastq.gz"),
            clean2=os.path.join(clean_dir, "{sample}_2.clean.fastq.gz")
        params:
            clean_dir=clean_dir
        shell:
            r"""
            set -euo pipefail
            install -d {params.clean_dir}
            # use readlink -f (portable) and quote paths
            ln -sf "$(readlink -f {input.r1})" "{output.clean1}" || cp -a "{input.r1}" "{output.clean1}"
            ln -sf "$(readlink -f {input.r2})" "{output.clean2}" || cp -a "{input.r2}" "{output.clean2}"
        """

