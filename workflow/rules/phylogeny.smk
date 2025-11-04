# workflow/rules/phylogeny.smk
# ------------------------------------------------------------
# Phylogenetic inference using IQ-TREE v2 and tree plotting with ETE3
# ------------------------------------------------------------

# Resolve paths & knobs from config with safe defaults
MSA_DIR    = config.get("results", {}).get("msa_dir", "results/msa")
PLOTS_DIR  = config.get("results", {}).get("plots_dir", "results/plots")
IQLOG_DIR  = config.get("logs", {}).get("iqtree", "logs/iqtree")
IQ_THREADS = config.get("parameters", {}).get("iqtree", {}).get("threads", 4)

# workflow/rules/phylogeny.smk

MSA_DIR    = config.get("results", {}).get("msa_dir", "results/msa")
PLOTS_DIR  = config.get("results", {}).get("plots_dir", "results/plots")
IQLOG_DIR  = config.get("logs", {}).get("iqtree", "logs/iqtree")
IQ_THREADS = config.get("parameters", {}).get("iqtree", {}).get("threads", 4)

rule run_iqtree:
    threads: IQ_THREADS
    input:
        aln = f"{MSA_DIR}/aligned_scaffolds.fasta"
    output:
        tree = f"{MSA_DIR}/tree.nwk"
    params:
        prefix  = f"{MSA_DIR}/iqtree",
        msa_dir = MSA_DIR,
        log_dir = IQLOG_DIR
    log:
        f"{IQLOG_DIR}/iqtree.log"
    conda:
        "../envs/iqtree.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.msa_dir}" "{params.log_dir}"

        # 0) Count sequences; require >=4 taxa
        nseq="$(grep -c '^>' "{input.aln}" 2>/dev/null || true)"
        [ -z "$nseq" ] && nseq=0
        if [ "$nseq" -lt 4 ]; then
          echo "ERROR: alignment has only $nseq sequences; IQ-TREE needs ≥4." >&2
          exit 2
        fi

        # 1) Sanitize: convert U/u -> T only on sequence lines
        tmp_aln="{params.msa_dir}/aligned_scaffolds.dna.fasta"
        sed '/^>/! s/[uU]/T/g' "{input.aln}" > "$tmp_aln"

        # 2) De-duplicate identical sequences (keeps first per sequence content)
        dedup_aln="{params.msa_dir}/aligned_scaffolds.dedup.fasta"
        seqkit rmdup -s "$tmp_aln" -o "$dedup_aln" 2> "{log}" || true

        nuniq="$(grep -c '^>' "$dedup_aln" 2>/dev/null || echo 0)"
        if [ "$nuniq" -lt 4 ]; then
          echo "ERROR: only $nuniq unique sequences remain after collapsing duplicates; need ≥4." >&2
          echo "Unique IDs kept:" >> "{log}"
          seqkit seq -n "$dedup_aln" >> "{log}" 2>/dev/null || true
          exit 2
        fi

        echo "=== IQ-TREE start ($(date)) ===" >> "{log}"
        echo "[info] nseq=$nseq, nuniq=$nuniq" >> "{log}"

        # 3) Try IQ-TREE with ModelFinder on deduplicated DNA, safe mode
        echo "[cmd] iqtree2 -s '$dedup_aln' -st DNA -m MFP -B 1000 --alrt 1000 -T {threads} -safe -pre '{params.prefix}' -redo" >> "{log}"
        if iqtree2 \
              -s "$dedup_aln" \
              -st DNA \
              -m MFP \
              -B 1000 \
              --alrt 1000 \
              -T {threads} \
              -safe \
              -pre "{params.prefix}" \
              -redo >> "{log}" 2>&1; then
          :
        else
          echo "[warn] MFP run failed; trying fixed model GTR+G…" >> "{log}"
          rm -f "{params.prefix}.treefile" || true

          # 3b) Fixed model fallback
          if iqtree2 \
                -s "$dedup_aln" \
                -st DNA \
                -m GTR+G \
                -B 1000 \
                --alrt 1000 \
                -T {threads} \
                -safe \
                -pre "{params.prefix}" \
                -redo >> "{log}" 2>&1; then
            :
          else
            echo "[warn] Fixed model also failed." >> "{log}"
            rm -f "{params.prefix}.treefile" || true

            # 3c) Optional gap trimming fallback (if trimal exists)
            if command -v trimal >/dev/null 2>&1; then
              trim_aln="{params.msa_dir}/aligned_scaffolds.trim.fasta"
              echo "[cmd] trimal -automated1 -in '$dedup_aln' -out '$trim_aln'" >> "{log}"
              trimal -automated1 -in "$dedup_aln" -out "$trim_aln" >> "{log}" 2>&1 || true
              ntrim="$(grep -c '^>' "$trim_aln" 2>/dev/null || echo 0)"
              if [ "$ntrim" -ge 4 ] && [ -s "$trim_aln" ]; then
                echo "[cmd] iqtree2 (trimmed) -m GTR+G" >> "{log}"
                iqtree2 \
                  -s "$trim_aln" \
                  -st DNA \
                  -m GTR+G \
                  -B 1000 \
                  --alrt 1000 \
                  -T {threads} \
                  -safe \
                  -pre "{params.prefix}" \
                  -redo >> "{log}" 2>&1 || true
              fi
            fi
          fi
        fi

        if [ ! -s "{params.prefix}.treefile" ]; then
          echo "ERROR: IQ-TREE finished but no treefile was produced." >&2
          exit 2
        fi

        mv "{params.prefix}.treefile" "{output.tree}"
        echo "=== IQ-TREE end ($(date)) ===" >> "{log}"
        """

# in rules/phylogeny.smk
rule clean_newick:
    input:  "results/msa/tree.nwk"
    output: "results/msa/tree.cleaned.nwk"
    shell:
        r'''
        # Keep only the first support value: ")X/Y(" or ")X/Y:" or ")X/Y," -> ")X("
        perl -pe 's/\)([0-9eE+\-\.]+)\/([0-9eE+\-\.]+)(?=[:\),])/)$$1/g' "{input}" > "{output}"
        '''

rule plot_ete3_tree:
    input:  tree="results/msa/tree.cleaned.nwk"
    output: svg="results/plots/tree.svg",
            pdf="results/plots/tree.pdf",
            png="results/plots/tree.png"
    conda:  "../envs/ete3.yaml"
    script: "../scripts/plot_ete3_tree.py"






