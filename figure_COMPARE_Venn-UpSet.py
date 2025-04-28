#######################################################################################################################################################
# Proteomic insights into Novel Food Insects: Homology-based proteome characterization and allergenicity considerations for EU-regulated insect species 
# Tobias Meisinger, Hannes Planatscher, Albert Braeuning, Eva-Maria Ladenburger, Dieter Stoll, Cristiano Garino, Hermann Broll, Oliver Poetz
# 2025

# Script for preparation of Figure 10, UpSet plot analysis, Venn analysis of COMPARE MASCOT results.
# Instructions: 1. Install all required packages, as defined in the import section
#               2. Save cleaned species-specific MASCOT result CSV files and this script in the same directory.
#               3. Run from CMD (Windows) or bash (Linux): 
#                   python figure_COMPARE_Venn-UpSet.py
#######################################################################################################################################################

# Imports
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import venn
import fitz
from upsetplot import from_memberships, plot
import tempfile
import shutil
from itertools import combinations

# Helper function to read and process species-specific CSV files from MASCOT
def process_csv_files():
    csv_files = [f for f in os.listdir() if f.endswith('.csv')]
    combined_df = pd.DataFrame()
    for file in csv_files:
        df = pd.read_csv(file)
        df['file_name'] = file
        df['origin_species'] = file.split('_')[0]
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    return combined_df

# Helper function to create a species dictionary
def create_species_dict(df, species_list):
    filtered_df = df[df['origin_species'].isin(species_list)]
    species_dict = {
        species: set(filtered_df[filtered_df['origin_species'] == species]['prot_acc'])
        for species in species_list
    }
    return species_dict

# Define color palettes for Venn diagrams
bright_colors_abc = ["#4477AA", "#228833", "#AA3377"]
bright_colors_def = ["#66CCEE", "#EE6677", "#CCBB44"]

# Set up the figure for Venn diagrams (with two subplots)
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Process CSV files
combined_df = process_csv_files()

# Helper to compute intersections manually for Venn diagrams
def compute_intersections(species_dict):
    sets = list(species_dict.values())
    labels = list(species_dict.keys())
    intersection_dict = {}

    n = len(sets)
    for i in range(1, n+1):
        for combo in combinations(range(n), i):
            # Elements present in all sets of the combo
            intersect = sets[combo[0]].copy()
            for idx in combo[1:]:
                intersect &= sets[idx]
            # Exclude elements that are also in other sets not in the combo
            other_indices = set(range(n)) - set(combo)
            for idx in other_indices:
                intersect -= sets[idx]
            if intersect:
                key = '&'.join([labels[idx] for idx in combo])
                intersection_dict[key] = list(intersect)
    return intersection_dict

# Subplot A: Venn for ACHDO, GRYSI, LOCMI
species_list_abc = ['ACHDO', 'GRYSI', 'LOCMI']
species_dict_abc = create_species_dict(combined_df, species_list_abc)
venn.venn(species_dict_abc, cmap=ListedColormap(bright_colors_abc), ax=axes[0], legend_loc="lower right")

# Export CSV for Venn A
intersections_abc = compute_intersections(species_dict_abc)
max_len_abc = max(len(v) for v in intersections_abc.values())
export_df_abc = pd.DataFrame({k: v + [''] * (max_len_abc - len(v)) for k, v in intersections_abc.items()})
export_df_abc.to_csv("SuppT10_Venn_COMPARE_Orthopteran_intersections.csv", index=False)

# Subplot B: Venn for ALPDA, HERIL, TENMO
species_list_def = ['ALPDA', 'HERIL', 'TENMO']
species_dict_def = create_species_dict(combined_df, species_list_def)
venn.venn(species_dict_def, cmap=ListedColormap(bright_colors_def), ax=axes[1], legend_loc="lower right")

# Export CSV for Venn B
intersections_def = compute_intersections(species_dict_def)
max_len_def = max(len(v) for v in intersections_def.values())
export_df_def = pd.DataFrame({k: v + [''] * (max_len_def - len(v)) for k, v in intersections_def.items()})
export_df_def.to_csv("SuppT11_Venn_COMPARE_Tenebrioninae_intersections.csv", index=False)


# Adjust font size for better readability
for ax in axes:
    for text in ax.texts:
        text.set_fontsize(15)
       
       
axes[0].text(0.05, 0.9, 'B', transform=axes[0].transAxes, fontsize=20, fontweight='bold')      
axes[1].text(0.1, 0.9, 'C', transform=axes[1].transAxes, fontsize=20, fontweight='bold')

# Adjust layout
plt.subplots_adjust(wspace=-0.6)

# Create temporary file for Venn diagram and save it
with tempfile.NamedTemporaryFile(delete=False, suffix='.pdf') as tmp_venn_file:
    plt.tight_layout()
    tmp_venn_file_path = tmp_venn_file.name
    plt.savefig(tmp_venn_file_path, dpi=1200, bbox_inches='tight')

# Create and save the UpSet plot
prot_to_species = combined_df.groupby('prot_acc')['origin_species'].apply(set).reset_index()
memberships = prot_to_species['origin_species'].tolist()
data = from_memberships(memberships)

# Set global font size for UpSet plot
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 14
})

# Create and customize UpSet plot
fig_upset = plt.figure(figsize=(12, 6))
axes_upset = plot(data,
                  subset_size='count',
                  show_counts='%d',
                  sort_by='degree',
                  orientation='horizontal',
                  intersection_plot_elements=10,
                  totals_plot_elements=5,
                  fig=fig_upset,
                  element_size=15)

# Add bold "A" in the top-left corner using fig.text
fig_upset.text(0.09, 0.98, "A", fontsize=20, fontweight='bold', color='black')

# Rotate and adjust positioning of bar-top labels in the intersection bar plot
if 'intersections' in axes_upset:
    ax_upset = axes_upset['intersections']
    for text in ax_upset.texts:
        text.set_rotation(90)
        x, y = text.get_position()
        text.set_position((x, y + 2))

if 'totals' in axes_upset:
    ax_upset = axes_upset['totals']
    for text in ax_upset.texts:
        x, y = text.get_position()
        text.set_position((x + 2, y))

# Export CSV of entries for each intersection group
grouped = prot_to_species.groupby(prot_to_species['origin_species'].apply(lambda x: '&'.join(sorted(x))))
grouped_dict = grouped['prot_acc'].apply(list).to_dict()

# Normalize lengths for CSV export
max_len = max(len(v) for v in grouped_dict.values())
export_df = pd.DataFrame({k: v + [''] * (max_len - len(v)) for k, v in grouped_dict.items()})

export_df.to_csv("SuppT9_UpSet_COMPARE_intersections.csv", index=False)  

# Create temporary file for UpSet plot and save it
with tempfile.NamedTemporaryFile(delete=False, suffix='.pdf') as tmp_upset_file:
    tmp_upset_file_path = tmp_upset_file.name
    plt.savefig(tmp_upset_file_path, dpi=1200, bbox_inches='tight')

# Open the temporary PDFs using fitz
doc1 = fitz.open(tmp_venn_file_path)
doc2 = fitz.open(tmp_upset_file_path)

page1 = doc1.load_page(0)
page2 = doc2.load_page(0)

# Get the dimensions of both pages
rect1 = page1.rect
rect2 = page2.rect

# New page dimensions: width is the max of both, height is the sum of both heights
new_width = max(rect1.width, rect2.width)
new_height = rect1.height + rect2.height

# Create a new PDF to combine both pages
combined = fitz.open()
new_page = combined.new_page(width=new_width, height=new_height)

# Center the first plot (Venn diagram)
x_offset_venn = (new_width - rect1.width) / 2
top_rect = fitz.Rect(x_offset_venn, new_height - rect1.height, x_offset_venn + rect1.width, new_height)
new_page.show_pdf_page(top_rect, doc1, 0)

# Place the second plot (UpSet plot) offcenter (e.g., aligned to left)
x_offset_upset = 0
bottom_rect = fitz.Rect(x_offset_upset, 0, x_offset_upset + rect2.width, rect2.height)
new_page.show_pdf_page(bottom_rect, doc2, 0)

# Save the combined PDF
combined.save("Figure_10_COMPARE_UpSet-Venn-plots.pdf")

# Clean up the temporary files
os.remove(tmp_venn_file_path)
os.remove(tmp_upset_file_path)
