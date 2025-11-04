# workflow/rules/assembly.smk
# ------------------------------------------------------------
# De novo assembly, polishing, reference selection, patching, and final consensus
# ------------------------------------------------------------

import os, pandas as pd

# Sample list
SAMPLES = pd.read_csv(config["samples_tsv"], sep=r"\s+", header=0)["sample"].tolist()

# Directories and inputs
ASSEMBLY_DIR   = config["results"]["assembly_dir"]
READS_DIR      = config["results"]["qc_decont"]

# Absolute, space-safe critical paths
PROJECT_ROOT   = os.path.abspath(".")
CANDIDATE_REFS = os.path.abspath(config.get("candidate_refs") or config["reference"])
BLASTDB_DIR    = os.path.abspath("resources/blastdb")
DB_BASE        = os.path.join(BLASTDB_DIR, "candidate_db")

# ------------------------------------------------------------
# 1) Build BLAST database of candidate references
# ------------------------------------------------------------
rule build_blast_db:
    input:
        fasta=CANDIDATE_REFS
    output:
        nhr=os.path.join(BLASTDB_DIR, "candidate_db.nhr"),
        nin=os.path.join(BLASTDB_DIR, "candidate_db.nin"),
        nsq=os.path.join(BLASTDB_DIR, "candidate_db.nsq")
    conda: "../envs/blast.yaml"
    log:   "logs/blast/build_db.log"
    shell:
        r"""
        set -euo pipefail
        install -d "{BLASTDB_DIR}" "logs/blast"
        test -s "{input.fasta}" || {{ echo "ERROR: missing/empty {input.fasta}" >&2; exit 2; }}
        grep -q '^>' "{input.fasta}" || {{ echo "ERROR: no FASTA headers in {input.fasta}" >&2; exit 2; }}
        makeblastdb -in "{input.fasta}" -dbtype nucl -parse_seqids \
          -out "{DB_BASE}" -logfile "{log}"
        """
# ------------------------------------------------------------
# 2) De novo assembly with SPAdes
# ------------------------------------------------------------
rule spades_assemble:
    threads: config["parameters"]["spades"]["threads"]
    input:
        r1=lambda wc: f"{READS_DIR}/{wc.sample}_1.clean.fastq.gz",
        r2=lambda wc: f"{READS_DIR}/{wc.sample}_2.clean.fastq.gz"
    output:
        contigs=f"{ASSEMBLY_DIR}/{{sample}}/contigs/contigs.fasta"
    params:
        memory=config["parameters"]["spades"]["memory"]
    conda: "../envs/spades.yaml"
    log:   f"{config['logs']['spades']}/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        install -d "{ASSEMBLY_DIR}/{wildcards.sample}/contigs" "{config[logs][spades]}"
        spades.py --phred-offset 33 \
          -1 "{input.r1}" -2 "{input.r2}" \
          -o "{ASSEMBLY_DIR}/{wildcards.sample}/contigs" \
          -t {threads} -m {params.memory} \
          > "{log}" 2>&1
        """

# ------------------------------------------------------------
# 3) Polish assembly with Racon (two-step: overlaps then racon)
# ------------------------------------------------------------
# produce PAF (smaller than SAM) for racon
# 3a) Make SAM overlaps (short reads preset)
# 3a) Make PAF overlaps (short reads) from one uncompressed FASTQ containing R1+R2
# OVERLAPS: build a *clean* FASTQ with unique qnames, then map from that exact file
rule overlaps_minimap2:
    threads: config.get("parameters", {}).get("racon", {}).get("threads", 4)
    input:
        contigs=lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/contigs/contigs.fasta",
        r1     =lambda wc: f"{READS_DIR}/{wc.sample}_1.clean.fastq.gz",
        r2     =lambda wc: f"{READS_DIR}/{wc.sample}_2.clean.fastq.gz"
    output:
        paf=f"{ASSEMBLY_DIR}/{{sample}}/contigs/overlaps.paf",
        fq =f"{ASSEMBLY_DIR}/{{sample}}/contigs/reads_all.clean.fq"
    conda: "../envs/minimap2.yaml"
    shell: r"""
        set -euo pipefail
        install -d "{ASSEMBLY_DIR}/{wildcards.sample}/contigs"

        # 1) write R1 with '/1' suffix on qname (strip anything after first space)
        zcat "{input.r1}" \
          | awk 'NR%4==1{{split($0,a," "); sub(/^@/,"",a[1]); print "@"a[1]"/1"}} NR%4==2||NR%4==3||NR%4==0' \
          > "{output.fq}.tmp1"

        # 2) write R2 with '/2' suffix on qname
        zcat "{input.r2}" \
          | awk 'NR%4==1{{split($0,a," "); sub(/^@/,"",a[1]); print "@"a[1]"/2"}} NR%4==2||NR%4==3||NR%4==0' \
          > "{output.fq}.tmp2"

        # 3) combine (R1 then R2) to one uncompressed FASTQ
        cat "{output.fq}.tmp1" "{output.fq}.tmp2" > "{output.fq}"
        rm -f "{output.fq}.tmp1" "{output.fq}.tmp2"

        # 4) create PAF overlaps from that exact FASTQ
        minimap2 -t {threads} -x sr "{input.contigs}" "{output.fq}" > "{output.paf}" 2>/dev/null
    """

rule racon_polish:
    threads: config.get("parameters", {}).get("racon", {}).get("threads", 4)
    input:
        contigs = lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/contigs/contigs.fasta",
        overlaps= lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/contigs/overlaps.paf",
        reads   = lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/contigs/reads_all.clean.fq"
    output:
        polished=f"{ASSEMBLY_DIR}/{{sample}}/contigs/polished.fasta"
    conda: "../envs/racon.yaml"   # keep racon 1.4.x here
    shell: r"""
        set -euo pipefail
        racon --threads {threads} "{input.reads}" "{input.overlaps}" "{input.contigs}" > "{output.polished}"
    """

# ------------------------------------------------------------
# 4) Select best reference via BLAST
# ------------------------------------------------------------
rule blast_ref:
    threads: config.get("parameters", {}).get("blast", {}).get("threads", 1)
    input:
        polished=lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/contigs/polished.fasta",
        db_files=rules.build_blast_db.output
    output:
        best_ref=f"{ASSEMBLY_DIR}/{{sample}}/best_ref.txt"
    params:
        db_base=DB_BASE,
        blast_out=f"{ASSEMBLY_DIR}/{{sample}}/blast.out",
        script="scripts/find_best_reference.py"
    conda: "../envs/blast.yaml"
    log:   f"logs/blast/{{sample}}.blastn.log"
    shell:
        r"""
        set -euo pipefail
        install -d "logs/blast" "{ASSEMBLY_DIR}/{wildcards.sample}"
        blastn -query "{input.polished}" -db "{params.db_base}" \
          -outfmt 6 -max_target_seqs 1 -num_threads {threads} \
          > "{params.blast_out}" 2>> "{log}"
        python "{params.script}" "{params.blast_out}" "{params.db_base}" "{output.best_ref}" >> "{log}" 2>&1
        """

# ------------------------------------------------------------
# 5) Patch polished contigs to chosen reference (use samtools faidx)
# ------------------------------------------------------------
rule patch_with_ref:
    threads: config["parameters"]["ragtag"]["threads"]
    input:
        polished=lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/contigs/polished.fasta",
        best_ref=lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/best_ref.txt"
    output:
        patched=f"{ASSEMBLY_DIR}/{{sample}}/patched/scaffolds.fasta"
    conda: "../envs/ragtag.yaml"
    log:   f"{config['logs']['ragtag']}/{{sample}}.log"
    params:
        db_base="resources/blastdb/candidate_db",
        refs_fa=CANDIDATE_REFS
    shell: r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "{ASSEMBLY_DIR}/{wildcards.sample}/patched"
        python scripts/patch_with_ref.py \
          --polished "{input.polished}" \
          --best-ref "{input.best_ref}" \
          --db-base  "{params.db_base}" \
          --refs-fa  "{params.refs_fa}" \
          --outdir   "{ASSEMBLY_DIR}/{wildcards.sample}/patched" \
          --threads  {threads} \
          --log      "{log}"

        # RagTag writes ragtag.scaffold.fasta — expose it at the expected path:
        install -D "{ASSEMBLY_DIR}/{wildcards.sample}/patched/ragtag.scaffold.fasta" "{output.patched}"
        # or: ln -sf ragtag.scaffold.fasta "{output.patched}"
    """
# ------------------------------------------------------------
# 6) Generate final consensus with iVar (mpileup → ivar)
# ------------------------------------------------------------
rule final_consensus:
    threads: config.get("parameters", {}).get("ivar", {}).get("threads", 4)
    input:
        ref=lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/patched/scaffolds.fasta",
        fq1=lambda wc: f"{READS_DIR}/{wc.sample}_1.clean.fastq.gz",
        fq2=lambda wc: f"{READS_DIR}/{wc.sample}_2.clean.fastq.gz"
    output:
        consensus=f"{ASSEMBLY_DIR}/{{sample}}/final_consensus.fasta"
    conda: "../envs/ivar.yaml"  # needs minimap2, samtools, ivar
    log:   f"{config['logs']['ivar']}/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        install -d "{config[logs][ivar]}" "{ASSEMBLY_DIR}/{wildcards.sample}"
        # Map, sort, index
        minimap2 -ax sr "{input.ref}" "{input.fq1}" "{input.fq2}" 2>> "{log}" \
          | samtools sort -@ {threads} -O bam -o "{ASSEMBLY_DIR}/{wildcards.sample}/tmp.bam"
        samtools index "{ASSEMBLY_DIR}/{wildcards.sample}/tmp.bam" 2>> "{log}"
        # Mpileup -> iVar consensus
        samtools mpileup -aa -A -d 0 -Q 0 -f "{input.ref}" "{ASSEMBLY_DIR}/{wildcards.sample}/tmp.bam" 2>> "{log}" \
          | ivar consensus -p "{ASSEMBLY_DIR}/{wildcards.sample}/consensus" -q 20 -t 0.6 -m 10 2>> "{log}"
        mv "{ASSEMBLY_DIR}/{wildcards.sample}/consensus.fa" "{output.consensus}"
        rm -f "{ASSEMBLY_DIR}/{wildcards.sample}/tmp.bam" "{ASSEMBLY_DIR}/{wildcards.sample}/tmp.bam.bai"
        """


