# Usage: awk -v sample="sample_XX" -f scripts/pick_longest.awk input.fasta
BEGIN { RS=">"; ORS=""; max=0; rec="" }
NR>1 {
  n = split($0, L, /\r?\n/);
  seq = "";
  for (i=2; i<=n; i++) seq = seq L[i];
  gsub(/[^A-Za-z]/, "", seq);
  if (length(seq) > max) { max = length(seq); rec = ">" sample "\n" seq "\n" }
}
END { if (max>0) print rec; else exit 1 }
