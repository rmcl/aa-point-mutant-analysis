from typing import Dict, List
import matplotlib.pyplot as plt
import pandas as pd


def analyze_mutation(
    df_results: pd.DataFrame,
    position: int,
    starting_aa: Dict[int, str],
    columns: List[str] = ["total_score", "interface_total_score"]
):
    """
    Analyze the mutation at a specific position and visualize the results.
    Args:
        df_results (pd.DataFrame): DataFrame containing mutation results.
        position (int): The residue position to analyze.
        starting_aa (Dict[int, str]): Dictionary mapping positions to their starting amino acids.
        columns (List[str]): List of columns to visualize.
    """
    df = df_results.copy()

    # Get wild-type amino acid at this position
    wt_aa = starting_aa.get(position)
    if wt_aa is None:
        raise ValueError(f"No starting amino acid defined for position {position}.")

    print('STARTING', wt_aa)

    df_pos = df[df['residue_position'] == position].copy()
    if df_pos.empty:
        raise ValueError(f"No results found for position {position}.")

    # Locate the reference wild-type row
    ref_row = df_pos.loc[df_pos['mutated_aa'] == wt_aa]
    if ref_row.empty:
        raise ValueError(f"No reference row found for {wt_aa}{position} in df_results.")

    # Compute deltas and colors
    for col in columns:
        ref_value = ref_row[col].values[0]
        df_pos[f'delta_{col}'] = df_pos[col] - ref_value

    # Define colors for bars
    base_colors = ['gray' if aa != wt_aa else 'black' for aa in df_pos['mutated_aa']]
    alt_colors = ['orange' if aa != wt_aa else 'darkorange' for aa in df_pos['mutated_aa']]

    num_cols = len(columns)
    fig, axs = plt.subplots(2, num_cols, figsize=(7 * num_cols, 10), squeeze=False)

    for i, col in enumerate(columns):
        # Original scores
        axs[0, i].bar(df_pos['mutated_aa'], df_pos[col], color=base_colors if 'interface' not in col else alt_colors)
        axs[0, i].set_title(f'{col.replace("_", " ").title()} at {wt_aa}{position}')
        axs[0, i].set_ylabel(col.replace("_", " ").title())
        axs[0, i].set_xlabel('Mutated Amino Acid')

        # Delta scores
        delta_col = f'delta_{col}'
        axs[1, i].bar(df_pos['mutated_aa'], df_pos[delta_col], color=base_colors if 'interface' not in col else alt_colors)
        axs[1, i].set_title(f'Δ {col.replace("_", " ").title()} (vs {wt_aa}{position})')
        axs[1, i].set_ylabel(f'Δ {col.replace("_", " ").title()}')
        axs[1, i].axhline(0, color='black', linestyle='--')
        axs[1, i].set_xlabel('Mutated Amino Acid')

        for row in (0, 1):
            axs[row, i].tick_params(axis='x', rotation=0)

    plt.tight_layout()
    plt.show()
