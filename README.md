# Amino Acid Point Mutant Analysis

Mutate residues throughout a protein-ligand complex and compute scores like delta delta G.


## How to run

```
python ./aa-point-mutant-analysis/perform_score.py <pos>,<residue aa> <output-dir>

python ./aa-point-mutant-analysis/perform_score.py 256,A ./outputs/

## Using slurm job

check out slurm.sbatch in the root!
