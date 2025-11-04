# workflow/rules/shannon_variability.smk

rule compute_variability:
    input:
        msa = config["results"]["msa_dir"] + "/aligned_scaffolds.fasta"
    output:
        txt      = config["results"]["variability_dir"] + "/variability.txt",            # <- match rule all
        windowed = config["results"]["variability_dir"] + "/variability_windowed.txt",
        plot     = config["results"]["variability_dir"] + "/variability_plot.png"
    threads: 1
    conda:
        "../envs/variability.yaml"
    log:
        config.get("logs", {}).get("variability", "logs/variability") + "/shannon.log"
    script:
        "../scripts/shannon_variability.py"   # <- one level up from rules/




