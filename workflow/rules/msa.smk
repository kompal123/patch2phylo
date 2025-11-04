# ------------------------------------------------------------
# Multi-sequence alignment of the patched scaffolds
# Uses seqkit to pick the longest scaffold per sample
# Aligns with MAFFT (faster & lighter); optional ClustalO try
# ------------------------------------------------------------

import os, csv

# If importing pandas in the Snakefile previously caused the NumPy ABI error,
# keep this pure-Python reader. If you *must* use pandas, see note below.
def read_samples(tsv_path):
    with open(tsv_path, newline="") as fh:
        # works with tab OR multiple spaces
        first = fh.readline().strip()
        sep = "\t" if "\t" in first else None  # None = split on whitespace
    samples = []
    with open(tsv_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t" if sep == "\t" else " ", skipinitialspace=True):
            s = row.get("sample") or row.get("Sample") or row.get("samples") or row.get("Samples")
            if s:
                samples.append(s.strip())
    if not samples:
        # Fallback: manual parse
        with open(tsv_path) as fh:
            hdr = fh.readline().strip().split()
            idx = hdr.index("sample")
            for line in fh:
                if line.strip():
                    samples.append(line.split()[idx])
    return samples

SAMPLES = read_samples(config["samples_tsv"])

ASSEMBLY_DIR = config.get("results", {}).get("assembly_dir", "results/assembly")
MSA_DIR      = config.get("results", {}).get("msa_dir", "results/msa")
LOG_DIR      = config.get("logs", {}).get("clustalo", "logs/clustalo")

# Per-sample scaffold path template
SCAFF = lambda s: f"{ASSEMBLY_DIR}/{s}/patched/scaffolds.fasta"

rule concat_scaffolds:
    input:
        samples=lambda wc: [f"results/assembly/{s}/patched/scaffolds.fasta" for s in SAMPLES]
    output:
        combined=f"{MSA_DIR}/scaffolds_all.fasta"
    conda:
        "../envs/msa.yaml"   # provides seqkit; awk is from the system
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{MSA_DIR}" "{LOG_DIR}"
        : > "{output.combined}"

        for f in {input.samples}; do
          sample=$(basename "$(dirname "$(dirname "$f")")")  # sample_XX

          if [[ ! -s "$f" ]]; then
            echo "[warn] $f is missing/empty — skipping" >&2
            continue
          fi
          if ! grep -q '^>' "$f"; then
            echo "[warn] $f has no FASTA headers — skipping" >&2
            continue
          fi

          # Fast path with seqkit: pick the longest record, rename header to sample
          if seqkit sort -l -r "$f" 2>> "{LOG_DIR}/concat.log" \
              | seqkit head -n 1 2>> "{LOG_DIR}/concat.log" \
              | seqkit replace -p '^[^ ]+' -r "$sample" 2>> "{LOG_DIR}/concat.log" \
              >> "{output.combined}"; then
            continue
          fi

          echo "[warn] seqkit failed on $f — using awk fallback" >&2
          if ! awk -v sample="$sample" -f scripts/pick_longest.awk "$f" >> "{output.combined}"; then
            echo "[error] awk fallback also failed for $f" >&2
            exit 1
          fi
        done

        # must have at least one record
        if ! grep -q '^>' "{output.combined}"; then
          echo "ERROR: empty combined FASTA" >&2
          exit 1
        fi
        '''

rule run_clustalo:
    threads: 4
    input:
        combined=f"{MSA_DIR}/scaffolds_all.fasta"
    output:
        aligned=f"{MSA_DIR}/aligned_scaffolds.fasta"
    log:
        f"{LOG_DIR}/clustalo.log"
    conda:
        "../envs/msa.yaml"   # make sure this includes mafft (and optionally clustalo)
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{MSA_DIR}" "{LOG_DIR}"
        echo "=== MSA start ($(date)) ===" > "{log}"
        echo "[info] threads: {threads}" >> "{log}"
        echo "[info] input: {input.combined}" >> "{log}"

        echo "[cmd] mafft --auto --thread {threads} '{input.combined}' > '{output.aligned}'" >> "{log}"
        if timeout 30m mafft --auto --thread {threads} "{input.combined}" > "{output.aligned}" 2>> "{log}"; then
            echo "[ok] MAFFT finished" >> "{log}"
        else
            status=$?
            echo "[warn] MAFFT failed (code $status). Trying ClustalO…" >> "{log}"
            rm -f "{output.aligned}"
            if command -v clustalo >/dev/null 2>&1; then
                echo "[cmd] clustalo -i '{input.combined}' -o '{output.aligned}' --threads {threads} --force" >> "{log}"
                if clustalo -i "{input.combined}" -o "{output.aligned}" --threads {threads} --force >> "{log}" 2>&1; then
                    echo "[ok] ClustalO finished" >> "{log}"
                else
                    echo "ERROR: ClustalO also failed" >&2
                    exit 1
                fi
            else
                echo "ERROR: ClustalO not available and MAFFT failed" >&2
                exit 1
            fi
        fi

        if [[ ! -s "{output.aligned}" ]]; then
          echo "ERROR: empty alignment" >&2
          exit 1
        fi

        echo "=== MSA end ($(date)) ===" >> "{log}"
        '''




