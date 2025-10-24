# patch2phylo
End-to-end Snakemake pipeline that patches assemblies to references, builds consensuses, aligns sequences, infers a phylogeny (IQ-TREE2), plots it (ETE3), and summarizes QC with MultiQC.
patch2phylo

From draft assemblies to a phylogenetic tree — end-to-end, reproducible, and powered by Snakemake + Conda.
This workflow patches sample contigs against the best reference, maps reads, builds per-sample consensus, aligns all consensuses, infers a tree, plots it, and summarizes QC with MultiQC. It also computes per-position and sliding-window Shannon entropy across the alignment.

<p align="left">
  <img alt="Python" src="https://img.shields.io/badge/Python-3.10+-blue.svg">
  <img alt="R" src="https://img.shields.io/badge/linux-%3E=4.0-276DC3.svg">
  <img alt="Jupyter" src="https://img.shields.io/badge/Snakemake-%F0%9F%93%9A-orange.svg">
</p>

## Pipeline overview

Per sample
(Optional upstream) QC / decontaminate reads (fastp, etc.)
Select best reference for each sample via BLAST.
Patch/scaffold contigs to the best reference with RagTag (minimap2 aligner).
Index & map reads to patched scaffolds with bwa + samtools.
Polish / consensus (e.g., racon) → per-sample scaffolds.fasta and consensus.

Across samples
6. Concatenate top scaffold per sample; MSA with MAFFT (fallbacks from ClustalO handled).
7. Phylogeny with IQ-TREE2 (ModelFinder + UFboot + SH-aLRT; safe DNA mode).
8. Plot tree with ETE3 (robust to IQ-TREE’s UFboot/SH-aLRT labels).
9. Variability: compute per-site entropy + sliding-window entropy and a PNG plot.
10. MultiQC summary report.

## Requirements

Linux/macOS with Bash
mamba or conda (miniforge/mambaforge recommended)

Snakemake ≥7
A POSIX shell (/usr/bin/bash)
Important: Configure strict channe l priority once:

conda config --set channel_priority strict

## Repository Layout
```
patch2phylo/
├─ README.md
├─ config/
│  └─ config.yaml
├─ workflow/
│  ├─ Snakefile
│  ├─ rules/
│  │  ├─ assembly.smk
│  │  ├─ mapping_and_consensus.smk
│  │  ├─ msa.smk
│  │  ├─ phylogeny.smk
│  │  └─ multiqc.smk
│  ├─ envs/
│  │  ├─ ragtag.yaml
│  │  ├─ bwa.yaml
│  │  ├─ msa.yaml
│  │  ├─ iqtree.yaml
│  │  └─ ete3.yaml
│  └─ scripts/
│     ├─ patch_with_ref.py
│     ├─ shannon_variability.py
│     └─ plot_ete3_tree.py
├─ resources/
│  └─ (reference FASTAs, BLAST DBs, etc.)
└─ logs/
   └─ (pipeline logs)

```

Quick start
##  clone
git clone https://github.com/<you>/patch2phylo.git
cd patch2phylo

##  Inspect and edit config
 - config/config.yaml
 - config/samples.tsv
 - resources/blastdb/ (candidate references as BLAST DB + FASTA)
 - data/references/candidate_refs.fasta

## Dry-run
snakemake -np

## Run (adjust cores and conda prefix as needed)
snakemake --use-conda --cores 8 \
  --conda-prefix /home/$USER/.snakemake/conda \
  --rerun-incomplete --latency-wait 60 --printshellcmds

Inputs
config/samples.tsv (tab or whitespace-separated)
...


| Component                                                | What it does                                                                                                                   |
| -------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| **config/config.yaml**                                   | Central knobs: directories (`results`, `logs`), threads, and parameters (RagTag, MSA/MAFFT, IQ-TREE, variability window/step). |
| **samples.tsv**                                          | Lists `sample` IDs; may include `r1`/`r2` columns for paired reads used by mapping/consensus rules.                            |
| **resources/blastdb/candidate_db.***                     | BLAST database of candidate references used to select the best reference per sample.                                           |
| **data/references/candidate_refs.fasta**                 | Same references as FASTA (required by RagTag).                                                                                 |
| **results/assembly/<sample>/patched/scaffolds.fasta**    | Patched scaffolds per sample produced by RagTag (`--aligner mm2`, safe minimap2 params).                                       |
| **results/mapping/<sample>.sorted.bam(.bai)**            | Read mapping (BWA-MEM → samtools sort/index) if reads are provided in `samples.tsv`.                                           |
| **results/consensus/<sample>.fa**                        | Per-sample consensus sequence (optional rule).                                                                                 |
| **results/msa/scaffolds_all.fasta**                      | Concatenated per-sample sequences (longest scaffold per sample) for alignment.                                                 |
| **results/msa/aligned_scaffolds.fasta**                  | Multiple sequence alignment (MAFFT `--auto`; falls back from ClustalO if needed).                                              |
| **results/msa/tree.nwk**                                 | Maximum-likelihood phylogeny (IQ-TREE2 with ModelFinder + UFboot 1000 + SH-aLRT 1000; DNA mode; `U/u→T` sanitize).             |
| **results/plots/tree.{svg,pdf,png}**                     | Tree rendered with ETE3; plotting script normalizes IQ-TREE support labels (e.g., `100/100:0.01`) for parsing.                 |
| **results/shannon_variability/variability.txt**          | Per-site Shannon entropy across the MSA.                                                                                       |
| **results/shannon_variability/variability_windowed.txt** | Sliding-window mean entropy (window/step from config).                                                                         |
| **results/shannon_variability/variability_plot.png**     | Plot of per-site (optional) and windowed entropy.                                                                              |
| **results/multiqc_report.html**                          | Combined QC report across logs/metrics (MultiQC).                                                                              |
| **scripts/plot_ete3_tree.py**                            | Reads `tree.nwk`, fixes IQ-TREE branch labels, renders SVG/PDF/PNG (ETE3, headless).                                           |
| **scripts/shannon_variability.py**                       | Computes entropy tables and plot (Biopython + matplotlib/Agg); robust to gaps.                                                 |
| **envs/ragtag.yaml**                                     | `ragtag`, `minimap2`, `seqkit`, `blast` (pins to stable combos).                                                               |
| **envs/bwa.yaml**                                        | `bwa 0.7.19`, `samtools 1.18` with compatible `libzlib`/`zlib` pins to avoid solver/ABI issues.                                |
| **envs/msa.yaml**                                        | `mafft` (primary), optional `clustalo`, plus `seqkit` for FASTA ops.                                                           |
| **envs/iqtree.yaml**                                     | `iqtree2` for ML tree inference with ModelFinder/UFboot.                                                                       |
| **envs/ete3.yaml**                                       | `ete3` (+ deps like `reportlab`, matplotlib) for tree plotting (headless).                                                     |
| **envs/qc.yaml**                                         | `multiqc` (and `fastp` if you aggregate read QC).                                                                              |

## Conda environments

Environment YAMLs live in envs/. We pin a few versions to avoid common ABI and zlib conflicts:
pandas/numpy versions compatible (prevents “numpy.dtype size changed” errors).
samtools and zlib/libzlib pins to satisfy solver compatibility.
ete3, mafft, iqtree, seqkit, blast, ragtag, minimap2, bwa, racon, multiqc.

If you see solving errors, try with mamba:

snakemake --use-conda --cores 8 --conda-frontend mamba
e workflow uses --latency-wait. You can increase it, e.g. --latency-wait 120.

## Reproducibility

All steps are versioned through Conda envs.
Snakemake tracks code, inputs, and conda env definitions; changing any triggers appropriate re-runs.

## Citation / credits

Please cite the tools you use:
Snakemake; minimap2; RagTag; BLAST+; bwa / samtools; racon; seqkit; MAFFT / Clustal Omega; IQ-TREE2; ETE3; MultiQC; Biopython; matplotlib.

## License

MIT (or your preferred license)

## Contact

Issues and questions: open a GitHub issue or email the maintainer.
