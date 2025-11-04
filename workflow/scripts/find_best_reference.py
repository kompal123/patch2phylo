import sys
import subprocess
from collections import defaultdict

def main(blast_file, db_path, output_file):
    # Parse BLAST results
    reference_scores = defaultdict(float)
    
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            ref_id = parts[1]  # Subject ID
            bitscore = float(parts[11])  # Bitscore
            reference_scores[ref_id] += bitscore
    
    # Find best reference
    if not reference_scores:
        print("Error: No BLAST matches found!")
        sys.exit(1)
    
    best_ref = max(reference_scores, key=reference_scores.get)
    print(f"Selected best reference: {best_ref} with score {reference_scores[best_ref]:.1f}")
    
    # Extract best reference sequence using blastdbcmd
    try:
        result = subprocess.run(
            ["blastdbcmd", "-db", db_path, "-entry", best_ref],
            capture_output=True,
            text=True,
            check=True
        )
        
        with open(output_file, 'w') as f_out:
            f_out.write(result.stdout)
            
        print(f"Successfully extracted {best_ref} to {output_file}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error extracting sequence: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: blastdbcmd not found. Make sure BLAST+ is installed.")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python find_best_reference.py <blast.out> <blast_db_path> <output.fasta>")
        print("Example: python find_best_reference.py sample1_blast.out blast_db/candidate_db references/sample1_best_ref.fasta")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
