# workflow/rules/mapping_and_consensus.smk
# ------------------------------------------------------------
# Indexing patched assemblies, mapping reads, and calling consensus
# ------------------------------------------------------------

import os

# Directories from config
READS_DIR = config["results"]["qc_decont"]
PATCH_DIR = config["results"]["assembly_dir"]     # patched assemblies under <sample>/patched
MAP_DIR   = config.get("results", {}).get("mapping_dir", "results/mapping")
CONS_DIR  = config.get("results", {}).get("consensus_dir", "results/consensus")

# Logs (with safe defaults if missing)
LOG_BWA       = config.get("logs", {}).get("bwa",       "logs/bwa")
LOG_BCFTOOLS  = config.get("logs", {}).get("bcftools",  "logs/bcftools")

# Threads (with safe defaults if missing)
BWA_THREADS      = config.get("parameters", {}).get("bwa",       {}).get("threads", 8)
BCFTOOLS_THREADS = config.get("parameters", {}).get("bcftools",  {}).get("threads", 4)

# ------------------------------------------------------------
rule bwa_index:
    """
    Index the patched reference once per sample
    """
    input:
        ref = f"{PATCH_DIR}/{{sample}}/patched/scaffolds.fasta"
    output:
        idx_done = touch(f"{PATCH_DIR}/{{sample}}/patched/.bwa_index.done")
    threads: 1
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        bwa index {input.ref}
        touch {output.idx_done}
        """

# ------------------------------------------------------------
rule map_reads:
    """
    Map cleaned reads to the patched assembly and produce sorted BAM + index
    """
    input:
        idx_done = rules.bwa_index.output.idx_done,
        ref      = f"{PATCH_DIR}/{{sample}}/patched/scaffolds.fasta",
        r1       = lambda wc: f"{READS_DIR}/{wc.sample}_1.clean.fastq.gz",
        r2       = lambda wc: f"{READS_DIR}/{wc.sample}_2.clean.fastq.gz"
    output:
        bam = f"{MAP_DIR}/{{sample}}.sorted.bam",
        bai = f"{MAP_DIR}/{{sample}}.sorted.bam.bai"
    threads: BWA_THREADS
    log:    f"{LOG_BWA}/{{sample}}.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        mkdir -p {MAP_DIR} "$(dirname {log})"
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} 2> {log} \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# ------------------------------------------------------------
rule consensus:
    """
    Call variants and generate consensus FASTA
    """
    input:
        bam = rules.map_reads.output.bam,
        ref = f"{PATCH_DIR}/{{sample}}/patched/scaffolds.fasta"
    output:
        fasta = f"{CONS_DIR}/{{sample}}.fasta"
    params:
        vcf = f"{MAP_DIR}/{{sample}}.calls.vcf.gz"
    threads: BCFTOOLS_THREADS
    log:    f"{LOG_BCFTOOLS}/{{sample}}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        mkdir -p {CONS_DIR} "$(dirname {log})"
        # 1) Generate VCF
        bcftools mpileup -f {input.ref} {input.bam} \
          | bcftools call -mv -Oz -o {params.vcf} 2>> {log}
        bcftools index {params.vcf} 2>> {log}
        # 2) Build consensus
        bcftools consensus -f {input.ref} {params.vcf} > {output.fasta} 2>> {log}
        # 3) Cleanup
        rm -f {params.vcf} {params.vcf}.csi
        """


