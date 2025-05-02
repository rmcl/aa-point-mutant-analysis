from typing import Dict, List
import pandas as pd

def compute_mutation_deltas(
    df_results: pd.DataFrame,
    starting_aa: Dict[int, str],
    columns: List[str] = ["total_score", "interface_total_score"]
) -> pd.DataFrame:
    """
    Compute delta values for all mutations relative to the wild-type amino acid at each position.

    Args:
        df_results (pd.DataFrame): DataFrame with mutation results.
        starting_aa (Dict[int, str]): Mapping from residue position to wild-type amino acid.
        columns (List[str]): List of score columns to compute deltas for.

    Returns:
        pd.DataFrame: New DataFrame with delta columns added, sorted by position and amino acid.
    """
    aa_order = list("ACDEFGHIKLMNPQRSTVWY")
    all_results = []

    for position, wt_aa in starting_aa.items():
        df_pos = df_results[df_results['residue_position'] == position].copy()
        if df_pos.empty:
            continue  # skip missing positions

        ref_row = df_pos.loc[df_pos['mutated_aa'] == wt_aa]
        if ref_row.empty:
            continue  # skip if no wild-type match

        for col in columns:
            ref_value = ref_row[col].values[0]
            df_pos[f'delta_{col}'] = df_pos[col] - ref_value

        # Reindex to ensure canonical amino acid order
        df_pos = df_pos.set_index('mutated_aa').reindex(aa_order).reset_index()
        df_pos['residue_position'] = position  # re-add in case it was dropped
        all_results.append(df_pos)

    # Combine and return
    df_all = pd.concat(all_results, ignore_index=True)
    return df_all
