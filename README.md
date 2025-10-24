# Amino Acid Point Mutant Analysis

Mutate residues throughout a protein-ligand complex and compute scores like delta delta G.


## How to run

```
python ./aa-point-mutant-analysis/perform_score.py <pos>,<residue aa> <output-dir>
```

```
python ./aa-point-mutant-analysis/perform_score.py 256,A ./outputs/
```

## Running with SLURM

The run.sbatch script in the project root launches 50 parallel tasks to evaluate point mutations. It requires an input text file listing mutation sites and target residues.

### Basic steps:

1. Create a input.txt file with lines of the following format:

```
256,A
257,A
258,A
259,A
```

2. sbatch run.batch
