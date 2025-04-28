#######################################################################################################################################################
# Proteomic insights into Novel Food Insects: Homology-based proteome characterization and allergenicity considerations for EU-regulated insect species 
# Tobias Meisinger, Hannes Planatscher, Albert Braeuning, Eva-Maria Ladenburger, Dieter Stoll, Cristiano Garino, Hermann Broll, Oliver Poetz
# 2025

# Script for preparation of Figure 4: Phylogenetic tree, Jaccard index dendrogram and compareMS2 2.0 dendrogram. 
# See section "Hierarchical clustering analysis" of the manuscript.
# Instructions: 1. Install all required packages, as defined in the import section
#               2. Save a newick tree file containing the phylogenetic tree derived from NCBI Taxonomy, a newick tree file derived from compareMS2 2.0 
#                  analysis, a CSV file containing a list of all unique peptides from all six species, and this script in the same directory.
#               3. Run from CMD (Windows) or bash (Linux): 
#                   python figure_combined-tree.py --phylo <newick phylogenetic tree> --dendro <CSV of unique peptides> --comparems2 <newick>
#######################################################################################################################################################

# Imports
import os
import argparse
import pandas as pd
from Bio import Phylo
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Label replacements for compareMS2 tree
LABEL_REPLACEMENTS = {
    "H._illucens": "Hermetia illucens",
    "T._molitor": "Tenebrio molitor",
    "L._migratoria": "Locusta migratoria",
    "A._domesticus": "Acheta domesticus",
    "G._sigillatus": "Gryllodes sigillatus",
    "A._diaperinus": "Alphitobius diaperinus"
}

def rename_terminals(tree):
    for clade in tree.get_terminals():
        if clade.name in LABEL_REPLACEMENTS:
            clade.name = LABEL_REPLACEMENTS[clade.name]

# Subplot A: Phylogenetic tree
def plot_phylogenetic_tree(newick_file, ax):
    tree = Phylo.read(newick_file, "newick")
    Phylo.draw(
        tree,
        do_show=False,
        label_func=lambda clade: clade.name.replace("_", " ") if clade.name else "",
        axes=ax,
        branch_labels=None,
    )
    for text in ax.texts:
        text.set_fontsize(16)
    ax.set_title("Phylogenetic Tree of Edible Insect Species", fontsize=18, ha='center', position=(0.6, 1))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
    ax.set_ylabel('')
    ax.set_xlabel('')

# Subplot B: Jaccard index calculation and creation of a dendrogram.
def plot_rotated_dendrogram(csv_file, ax):
    df = pd.read_csv(csv_file)
    df['presence'] = 1
    matrix = df.pivot_table(index='sample_taxonomy', columns='pep_seq', values='presence', fill_value=0)
    distances = pdist(matrix, metric='jaccard')
    linkage_matrix = linkage(distances, method='weighted')
    dendrogram(
        linkage_matrix,
        labels=matrix.index,
        orientation='left',
        ax=ax,
        leaf_font_size=16
    )
    ax.set_title("Hierarchical Clustering Dendrogram", fontsize=18, ha='center', position=(0.7, 1))
    ax.set_xlabel("distance", fontsize=14)
    ax.set_ylabel("species", fontsize=14)

# Subplot C: compareMS2 2.0 dendrogram creation.
def plot_comparems2_tree(nwk_file, ax):
    tree = Phylo.read(nwk_file, "newick")
    rename_terminals(tree)
    Phylo.draw(
        tree,
        do_show=False,
        axes=ax,
        branch_labels=None
    )
    ax.set_title("compareMS2 2.0 Dendrogram", fontsize=18, ha='center', position=(0.7, 1))
    ax.set_ylabel("species", fontsize=14)
    ax.set_xlabel("branch length", fontsize=14)
    ax.set_yticks([])
    ax.set_xlim(-0.5, 5.0)
    for text in ax.texts:
        text.set_fontsize(16)

# Main
def main():
    parser = argparse.ArgumentParser(description="Generate a 3-panel figure: phylogenetic tree, clustering dendrogram, and compareMS2 dendrogram.")
    parser.add_argument("-p", "--phylo", required=True, help="Newick file for subplot A (phylogenetic tree)")
    parser.add_argument("-d", "--dendro", required=True, help="CSV file for subplot B (hierarchical clustering)")
    parser.add_argument("-c", "--comparems2", required=True, help="Newick file for subplot C (compareMS2 dendrogram)")
    args = parser.parse_args()

    fig, axes = plt.subplots(1, 3, figsize=(24, 10))
    plt.subplots_adjust(wspace=0.7)

    plot_phylogenetic_tree(args.phylo, axes[0])
    axes[0].text(-0.12, 1.035, 'A', transform=axes[0].transAxes, fontsize=18, weight='bold', va='top', ha='left')

    plot_rotated_dendrogram(args.dendro, axes[1])
    axes[1].text(-0.12, 1.035, 'B', transform=axes[1].transAxes, fontsize=18, weight='bold', va='top', ha='left')

    plot_comparems2_tree(args.comparems2, axes[2])
    axes[2].text(-0.12, 1.035, 'C', transform=axes[2].transAxes, fontsize=18, weight='bold', va='top', ha='left')

    output_filename = "Figure_4_Arthropoda_Dendrogram"
    plt.savefig(f"{output_filename}.pdf", dpi=1200, bbox_inches='tight')
    print(f"Combined plot saved as {output_filename}.png and {output_filename}.pdf")

if __name__ == "__main__":
    main()
