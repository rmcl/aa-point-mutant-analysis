# Amino Acid Point Mutant Analysis

Mutate residues throughout a protein-ligand complex and compute scores like delta delta G.


## Installation

The main dependency is pyrosetta which you can get here: https://www.pyrosetta.org/downloads

## How to run

```
python ./aa-point-mutant-analysis/perform_score.py <pos>,<residue aa> <output-dir>
```

```
python ./aa-point-mutant-analysis/perform_score.py 256,A ./outputs/
```

### Output

The script outptus a jsonl (list of json dictionary) file with each line representing one mutation at a site. An example is in `examples/`.

```
{"mutated_aa": "A", "residue_position": 189, "total_score": -607.7893564554801, ...}
{"mutated_aa": "A", "residue_position": 190, "total_score": -607.7893564554801, ...}
...
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



## Questions?

Drop me a line... rmcl@uoregon.edu


