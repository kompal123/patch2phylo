#!/usr/bin/env python3
import argparse, os, re, subprocess, sys, shutil

def run(cmd, logfh, check=True):
    logfh.write(f"[cmd] {' '.join(cmd)}\n")
    logfh.flush()
    p = subprocess.run(cmd, stdout=logfh, stderr=logfh)
    if check and p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd)
    return p.returncode

def first_token(s: str) -> str:
    s = s.strip()
    s = re.sub(r"^>+", "", s)              # drop leading '>'
    s = re.split(r"\s+", s, maxsplit=1)[0] # keep token before whitespace
    return s

def file_has_header(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    with open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                return True
    return False

def main():
    ap = argparse.ArgumentParser(description="Extract best reference and run RagTag scaffold.")
    ap.add_argument("--polished", required=True, help="polished contigs FASTA")
    ap.add_argument("--best-ref", required=True, help="file containing best reference ID (one line)")
    ap.add_argument("--db-base",  required=True, help="BLAST DB base path (without .nin/.nsq)")
    ap.add_argument("--refs-fa",  required=True, help="Candidate references FASTA")
    ap.add_argument("--outdir",   required=True, help="Output directory (RagTag will write here)")
    ap.add_argument("--threads",  type=int, default=1, help="Threads for RagTag/minimap2")
    ap.add_argument("--log",      required=True, help="Path to log file")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.dirname(args.log), exist_ok=True)
    with open(args.log, "a", buffering=1) as logfh:
        logfh.write("=== patch_with_ref.py start ===\n")

        # 1) Read & sanitize best ID
        if not os.path.exists(args.best_ref) or os.path.getsize(args.best_ref) == 0:
            logfh.write("ERROR: best_ref file missing or empty\n")
            sys.exit(2)
        with open(args.best_ref, "rt") as fh:
            best_id = first_token(fh.readline())
        logfh.write(f"[info] best_id={best_id}\n")

        ref_fa = os.path.join(args.outdir, "ref.fasta")

        # 2) Try BLAST DB extraction, then fall back to FASTA grep
        try:
            run(["blastdbcmd", "-db", args.db_base, "-entry", best_id, "-out", ref_fa], logfh)
        except subprocess.CalledProcessError:
            logfh.write("[warn] blastdbcmd failed; trying FASTA grep fallback\n")
            # Try exact token (start-of-line, after '>')
            patt_exact = rf"^{re.escape(best_id)}\b"
            try:
                run(["seqkit", "grep", "-nr", "-p", patt_exact, args.refs_fa], logfh, check=True)
                # seqkit prints to stdout by default; we want to capture to file
                # Re-run directing to file to keep logging simple:
                with open(ref_fa, "wb") as outfh:
                    subprocess.check_call(
                        ["seqkit", "grep", "-nr", "-p", patt_exact, args.refs_fa],
                        stdout=outfh, stderr=logfh
                    )
            except Exception:
                # Try versionless (strip .version)
                base = best_id.split(".", 1)[0]
                logfh.write(f"[info] trying versionless token: {base}\n")
                patt_base = rf"^{re.escape(base)}(\.|$)"
                try:
                    with open(ref_fa, "wb") as outfh:
                        subprocess.check_call(
                            ["seqkit", "grep", "-nr", "-p", patt_base, args.refs_fa],
                            stdout=outfh, stderr=logfh
                        )
                except Exception:
                    logfh.write(f"ERROR: could not find {best_id} (or {base}) in {args.refs_fa}\n")
                    sys.exit(3)

        # 3) Normalize header: keep only first token after '>'
        if not file_has_header(ref_fa):
            logfh.write("ERROR: ref.fasta missing/invalid FASTA header\n")
            sys.exit(4)
        tmp_clean = ref_fa + ".clean"
        with open(ref_fa, "rt") as fh, open(tmp_clean, "wt") as out:
            for line in fh:
                if line.startswith(">"):
                    tok = first_token(line)
                    out.write(">" + tok + "\n")
                else:
                    out.write(line)
        os.replace(tmp_clean, ref_fa)

        # 4) Sanity: polished exists & looks like FASTA
        if not file_has_header(args.polished):
            logfh.write("ERROR: polished FASTA missing/invalid headers\n")
            sys.exit(5)

        # 5) Run RagTag (v2.1.0): use minimap2 and -t for threads
        ragtag_out = os.path.join(args.outdir, "ragtag.scaffold.fasta")
        if os.path.exists(ragtag_out):
            os.remove(ragtag_out)

        cmd = [
            "ragtag.py", "scaffold",
            ref_fa,
            os.path.abspath(args.polished),
            "-o", args.outdir,
            "--aligner", "minimap2",
            "-t", str(args.threads)
        ]
        run(cmd, logfh, check=False)

        # 6) Verify output
        if not (os.path.exists(ragtag_out) and os.path.getsize(ragtag_out) > 0):
            logfh.write("ERROR: RagTag did not produce ragtag.scaffold.fasta.\n")
            sys.exit(6)

        logfh.write("=== patch_with_ref.py done ===\n")

if __name__ == "__main__":
    main()



