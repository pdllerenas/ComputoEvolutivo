import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

def plot_combined_results(directory):
    # 1. Load the landscape data (final penalties)
    landscape_dfs = []
    for i in range(5):
        filename = f"{directory}/kfls_results_rank_{i}.csv"
        if os.path.exists(filename):
            landscape_dfs.append(pd.read_csv(filename))
        else:
            print(f"Warning: {filename} not found.")

    if not landscape_dfs:
        print("No landscape CSV files found.")
        return

    landscape_df = pd.concat(landscape_dfs, ignore_index=True)
    mean_landscape_df = landscape_df.groupby('k')['penalty'].mean().reset_index()

    # Automatically identify the best k based on the lowest mean penalty
    best_k = mean_landscape_df.loc[mean_landscape_df['penalty'].idxmin(), 'k']
    best_k = int(best_k)
    print(f"Optimal parameter found: k = {best_k}")

    # 2. Load the convergence data strictly for the best_k
    conv_dfs = []
    for i in range(5):
        filename = f"{directory}/convergence_rank_{i}.csv"
        if os.path.exists(filename):
            df = pd.read_csv(filename)
            df_k = df[df['k'] == best_k].copy()
            if not df_k.empty:
                df_k['global_run_id'] = f"Rank_{i}_" + df_k['run_id'].astype(str)
                conv_dfs.append(df_k)
        else:
            print(f"Warning: {filename} not found.")

    if not conv_dfs:
        print(f"No convergence data found for k={best_k}.")
        return

    conv_df = pd.concat(conv_dfs, ignore_index=True)
    mean_conv_df = conv_df.groupby('generation')['best_penalty'].mean().reset_index()

    # 3. Set up the 1x2 subplot figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # --- LEFT PLOT: k-Parameter Landscape ---
    ax1.scatter(landscape_df['k'], landscape_df['penalty'], color='black', s=3, marker='o', zorder=1, alpha=0.8)
    ax1.scatter(mean_landscape_df['k'], mean_landscape_df['penalty'], color='tab:red', s=15, marker='o', zorder=2)
    
    ax1.set_xlabel("k", fontsize=11)
    ax1.set_ylabel("penalty", fontsize=11)
    ax1.set_title("k-Parameter Sweep Landscape", fontsize=12)

    # --- RIGHT PLOT: Convergence Rate of Best k ---
    for run_id, group_df in conv_df.groupby('global_run_id'):
        ax2.plot(group_df['generation'], group_df['best_penalty'], color='black', linewidth=1.0, alpha=0.15, zorder=1)
    
    ax2.plot(mean_conv_df['generation'], mean_conv_df['best_penalty'], color='tab:red', linewidth=2.5, zorder=2)
    
    ax2.set_xlabel("Generation", fontsize=11)
    ax2.set_ylabel("Best Penalty", fontsize=11)
    ax2.set_title(f"Convergence Rate (Optimal k = {best_k})", fontsize=12)

    # --- APPLY STYLING TO BOTH AXES ---
    for ax in (ax1, ax2):
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.0)
        ax.tick_params(direction='out', color='black', labelsize=10)

    # 4. Save and Show
    plt.tight_layout()
    output_filename = f"{directory}_kfls_full_analysis.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Plot saved successfully as '{output_filename}'.")
    plt.show()

def main():
    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <directory>")
        return 1
    plot_combined_results(sys.argv[1])

if __name__ == "__main__":
    main()
