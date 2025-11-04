# workflow/rules/multiqc.smk

import os
import pandas as pd

# Load sample list (to establish dependencies)
samples_df = pd.read_csv(config["samples_tsv"], sep=r"\s+", header=0)
SAMPLES = samples_df["sample"].tolist()

rule multiqc:
    input:
        # Ensure rule waits on all Fastp JSONs
        jsons = expand(
            os.path.join(config["logs"]["fastp"], "{sample}_fastp.json"),
            sample=SAMPLES
        )
    output:
        report = "results/multiqc_report.html"
    threads: 1
    conda:
        "../envs/multiqc.yaml"
    log:
        os.path.join(config["logs"]["fastp"], "multiqc.log")
    shell:
        """
        mkdir -p results
        multiqc {config[logs][fastp]} \
            -o results \
            -n multiqc_report.html \
            > {log} 2>&1
        """
