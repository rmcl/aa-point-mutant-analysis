import random
from pyrosetta.toolbox import mutate_residue


def random_amino_acid():
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # Standard 20 amino acids
    return random.choice(amino_acids)

def create_mutant_pose(original_pose, positions_to_mutate):
    mutant_pose = original_pose.clone()

    changes_made = []

    for pos in positions_to_mutate:
        new_aa = random_amino_acid()
        changes_made.append(f'{mutant_pose.residue(pos).name1()}{pos}{new_aa}')
        mutate_residue(mutant_pose, pos, new_aa)

    return mutant_pose, changes_made

def create_mutated_pose_by_residue_dict(original_pose, mutated_residues):
    """Create a mutated pose from an original pose and a dictionary of mutated residues.

    Args:
        original_pose (pyrosetta.Pose): The original pose.
        mutated_residues (dict): A dictionary of mutated residues, where keys are residue indices and values are new amino acids.
    """
    mutant_pose = original_pose.clone()

    for residue_idx, new_aa in mutated_residues.items():
        mutate_residue(mutant_pose, residue_idx, new_aa)

    return mutant_pose

def compare_poses(pose1, pose2):
    assert pose1.total_residue() == pose2.total_residue(), "Poses have different number of residues."
    for i in range(1, pose1.total_residue() + 1):
        res1 = pose1.residue(i)
        res2 = pose2.residue(i)

        if res1.name() != res2.name():
            print(f"Residue {i}: {res1.name()} (pose1) vs {res2.name()} (pose2)")
