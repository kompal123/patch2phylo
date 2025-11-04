# workflow/rules/kraken.smk
# ------------------------------------------------------------
# Taxonomic screening with Kraken 2 (+ MultiQC summary).
# If config["kraken2_db"] is empty, the rules are skipped, but `screen` always exists.
# ------------------------------------------------------------

import os
import pandas as pd  # (consider replacing with csv.DictReader to avoid pandas parse-time ABI issues)

# Read samples
samples_list = pd.read_csv(config["samples_tsv"], sep=r"\s+", header=0)["sample"].tolist()

# Config parameters
kraken_db = config.get("kraken2_db", "")
out_dir   = config["results"]["kraken_dir"]
log_dir   = config["logs"]["kraken2"]
threads   = config["parameters"]["kraken2"]["threads"]

if kraken_db:

    rule kraken_classify:
        """
        Classify reads with Kraken2, generating a per-sample report.
        """
        input:
            r1=lambda wc: f"{config['results']['qc_trimmed']}/{wc.sample}_1.trimmed.fastq.gz",
            r2=lambda wc: f"{config['results']['qc_trimmed']}/{wc.sample}_2.trimmed.fastq.gz"
        output:
            report=os.path.join(out_dir, "{sample}.report"),
            kraken=os.path.join(out_dir, "{sample}.kraken")
        params:
            db=kraken_db
        threads: threads
        log:
            os.path.join(log_dir, "{sample}.log")
        conda:
            "../envs/kraken2.yaml"
        shell:
            r"""
            mkdir -p "$(dirname {output.report})" "$(dirname {log})"
            kraken2 --db {params.db} --threads {threads} --paired \
                --report {output.report} --output {output.kraken} \
                {input.r1} {input.r2} > {log} 2>&1
            """

    rule kraken_multiqc:
        """
        Aggregate all Kraken2 reports with MultiQC.
        """
        input:
            expand(os.path.join(out_dir, "{sample}.report"), sample=samples_list)
        output:
            html=os.path.join("results", "multiqc_kraken.html")
        log:
            os.path.join(log_dir, "multiqc.log")
        conda:
            "../envs/multiqc.yaml"
        shell:
            r"""
            mkdir -p "$(dirname {output.html})" "$(dirname {log})"
            multiqc {out_dir} -n "$(basename {output.html})" -o "$(dirname {output.html})" > {log} 2>&1
            """

    # Alias rule to run just screening
    rule screen:
        input:
            config["results"]["multiqc_report"]

else:

    rule screen:
        shell:
            r'echo "No Kraken2 database specified; skipping taxonomic screening."'

