import pandas as pd
import numpy as np
import argparse
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS, TSNE
from sklearn.decomposition import PCA
import warnings

# --- 1. Suppress warning messages ---
warnings.filterwarnings("ignore", category=UserWarning, module="umap")
warnings.filterwarnings("ignore", category=FutureWarning)

# Attempt to import umap-learn
try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

# --- Configure plotting parameters for Adobe Illustrator compatibility ---
import matplotlib
matplotlib.use('Agg')
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 7
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 6
plt.rcParams['axes.titlesize'] = 9

def parse_sample_name(name):
    """
    Parses the sample ID into metadata components.
    Expected format: Species_Haplotype_Chromosome
    """
    parts = name.split('_')
    if len(parts) < 3:
        return name, "Unknown", "Unknown", "Unknown"
    chromosome = parts[-1]
    haplotype = parts[-2]
    species = "_".join(parts[:-2])
    genome = species.split('_')[0]
    return species, genome, haplotype, chromosome

def load_config(config_path):
    """
    Loads color and naming configurations from a TSV file.
    Expected format: ID \t DisplayName \t HexColor
    """
    try:
        df = pd.read_csv(config_path, sep='\t', header=None, names=['ID', 'Name', 'Color'])
        df['ID'] = df['ID'].astype(str).str.strip()
        df['Name'] = df['Name'].astype(str).str.strip()
        df['Color'] = df['Color'].astype(str).str.strip()
        id_to_color = dict(zip(df['ID'], df['Color']))
        id_to_name = dict(zip(df['ID'], df['Name']))
        return id_to_color, id_to_name
    except Exception as e:
        print(f"[Error] Failed to load config file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Generate Main Figure (Dimensionality Reduction & Visualization)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("input_matrix", help="Path to the Mash distance matrix (TSV).")
    parser.add_argument("config_file", help="Path to the configuration file (ID, Name, Color).")
    parser.add_argument("output_dir", help="Directory to save output files.")
    
    parser.add_argument("--method", default="umap", 
                        choices=['umap', 'nmds', 'mds', 'pca', 'tsne'], 
                        help="Dimensionality reduction method (default: umap)")
    
    # --- 2. Adjustable UMAP parameters ---
    parser.add_argument("--umap_neighbors", type=int, default=15, 
                        help="UMAP: n_neighbors (default: 15). Higher = global structure, Lower = local structure.")
    parser.add_argument("--umap_mindist", type=float, default=0.1, 
                        help="UMAP: min_dist (default: 0.1). Higher = more spread out, Lower = tighter clusters.")
    
    parser.add_argument("--width", type=float, default=12.0, help="Figure width (inches).")
    parser.add_argument("--height", type=float, default=4.0, help="Figure height (inches).")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print(f"Loading config: {args.config_file} ...")
    color_map, name_map = load_config(args.config_file)

    print(f"Reading matrix: {args.input_matrix} ...")
    df = pd.read_csv(args.input_matrix, sep='\t', index_col=0)
    
    # Generate metadata from index
    meta_list = []
    for sample in df.index:
        spec, gen, hap, chr_name = parse_sample_name(sample)
        meta_list.append({
            'SampleID': sample, 'Species': spec, 'Genome': gen,
            'Haplotype': hap, 'Chromosome': chr_name
        })
    meta_df = pd.DataFrame(meta_list).set_index('SampleID')
    
    matrix_data = df.values
    print(f"Running {args.method.upper()} transformation...")

    coords = None
    eig_pct = None

    if args.method == 'umap':
        if not HAS_UMAP:
            print("[ERROR] umap-learn is not installed.")
            sys.exit(1)
        print(f"  -> Parameters: n_neighbors={args.umap_neighbors}, min_dist={args.umap_mindist}")
        reducer = umap.UMAP(
            metric='precomputed', 
            n_neighbors=args.umap_neighbors, 
            min_dist=args.umap_mindist, 
            random_state=42
        )
        coords = reducer.fit_transform(matrix_data)

    elif args.method == 'nmds':
        mds = MDS(n_components=2, dissimilarity='precomputed', metric=False, 
                  n_init=50, max_iter=1000, random_state=42, n_jobs=-1)
        coords = mds.fit_transform(matrix_data)
        print(f"  -> NMDS Stress: {mds.stress_:.4f}")

    elif args.method == 'mds':
        mds = MDS(n_components=2, dissimilarity='precomputed', metric=True, random_state=42)
        coords = mds.fit_transform(matrix_data)
        # Using PCA to estimate variance for axis labels in metric MDS
        pca_temp = PCA(n_components=2).fit(matrix_data)
        eig_pct = pca_temp.explained_variance_ratio_ * 100

    elif args.method == 'tsne':
        perp = min(30, matrix_data.shape[0] - 1)
        coords = TSNE(n_components=2, metric='precomputed', random_state=42, init='random', perplexity=perp).fit_transform(matrix_data)

    elif args.method == 'pca':
        pca = PCA(n_components=2)
        coords = pca.fit_transform(matrix_data)
        eig_pct = pca.explained_variance_ratio_ * 100

    # Prepare final result dataframe
    result_df = meta_df.copy()
    result_df['Dim1'] = coords[:, 0]
    result_df['Dim2'] = coords[:, 1]
    result_df['DisplayName'] = result_df['Species'].map(name_map).fillna(result_df['Species'])
    
    # Save coordinates
    coords_path = os.path.join(args.output_dir, f"{args.method}_coordinates.tsv")
    result_df.to_csv(coords_path, sep='\t')

    # --- Plotting ---
    print("Generating publication-quality figure...")
    aa_df = result_df[result_df['Genome'] == 'AA'].copy()
    
    fig, axes = plt.subplots(1, 3, figsize=(args.width, args.height), dpi=300)
    plt.subplots_adjust(wspace=0.3)

    # Determine axis labels based on method
    if eig_pct is not None:
        prefix = "PC" if args.method == 'pca' else "PCoA"
        xlab = f"{prefix} 1 ({eig_pct[0]:.1f}%)"
        ylab = f"{prefix} 2 ({eig_pct[1]:.1f}%)"
    else:
        if args.method == 'umap':
            xlab = f"UMAP 1 (n={args.umap_neighbors})"
            ylab = f"UMAP 2"
        else:
            xlab = f"{args.method.upper()} 1"
            ylab = f"{args.method.upper()} 2"

    def plot_panel(ax, data, hue_col, title):
        if data.empty: return
        sns.scatterplot(
            data=data, x='Dim1', y='Dim2', 
            hue=hue_col, palette=color_map, 
            s=25, alpha=0.9, ax=ax, 
            linewidth=0.2, edgecolor='black', legend='full'
        )
        ax.set_title(title, fontweight='bold', pad=10)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        
        # Format the legend with configured display names
        handles, labels = ax.get_legend_handles_labels()
        new_labels = [name_map.get(lbl, lbl) for lbl in labels]
        ax.legend(handles, new_labels, title="", frameon=False, 
                  bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0., fontsize=6)

    # Generate the three panels
    plot_panel(axes[0], result_df, 'Species', f"All Genomes ({args.method.upper()})")
    plot_panel(axes[1], aa_df, 'Chromosome', "AA Genome (By Chromosome)")
    plot_panel(axes[2], aa_df, 'Species', "AA Genome (By Species)")

    # Save as PDF for vector editing
    pdf_out = os.path.join(args.output_dir, f"{args.method}_main_figure.pdf")
    plt.savefig(pdf_out, format='pdf', bbox_inches='tight')
    print(f"Done! Plot saved to: {pdf_out}")

if __name__ == "__main__":
    main()