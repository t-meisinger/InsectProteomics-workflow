#######################################################################################################################################################
# Proteomic insights into Novel Food Insects: Homology-based proteome characterization and allergenicity considerations for EU-regulated insect species 
# Tobias Meisinger, Hannes Planatscher, Albert Braeuning, Eva-Maria Ladenburger, Dieter Stoll, Cristiano Garino, Hermann Broll, Oliver Poetz
# 2025

# Script for validation of putative allergens derived from MASCOT searches against the COMPARE reference database. 
# See section "Allergen candidate identification" of the manuscript.
# Instructions: 1. Install all required packages, as defined in the import section
#               2. Run the MASCOT-csv_clean.py script to generate cleaned CSV files from MASCOT export CSV files.
#               3. Concatenate all cleaned CSV files and add a species column to track sample origin.
#               4. Run from CMD (Windows) or bash (Linux): 
#                   python allergen_80AA-validation.py <concatenated-csv> <output.csv>
#######################################################################################################################################################

# Imports
import pandas as pd
import argparse

def filter_peptides(input_csv, output_csv):
    # Load the CSV file
    data = pd.read_csv(input_csv)
    
    # Rename columns for easier handling
    data.columns = data.columns.str.strip()  # Remove leading/trailing spaces
    data = data.rename(columns={
        data.columns[0]: "Species",
        data.columns[4]: "ProteinAccession",
        data.columns[32]: "PeptideSequence",
        data.columns[24]: "StartPosition",
        data.columns[25]: "EndPosition"
    })
    
    # Group by Species and ProteinAccession
    grouped = data.groupby(["Species", "ProteinAccession"])
    
    # Function to find clusters within an 80 AA window and calculate identity
    def find_clusters(group):
        group = group.sort_values(by="StartPosition")
        cluster = []
        clusters = []
        
        # Iterate through peptides
        for _, row in group.iterrows():
            if not cluster:
                cluster = [row]
            else:
                # Check if current peptide fits in the 80 AA window (start + 80 limit)
                if row["StartPosition"] <= cluster[0]["StartPosition"] + 79:
                    cluster.append(row)
                else:
                    # Finalize the current cluster
                    clusters.append(cluster)
                    cluster = [row]  # Start a new cluster
        
        # Check the last cluster
        if cluster:
            clusters.append(cluster)
        
        # Calculate cluster stats and prepare output
        output = []
        for cluster in clusters:
            if len(cluster) > 1:  # Only keep clusters with multiple peptides
                df_cluster = pd.DataFrame(cluster)
                cluster_start = df_cluster["StartPosition"].min()
                cluster_end = cluster_start + 79  # Fixed 80 AA window
                
                # Adjust peptides to only count overlap within the 80 AA window
                peptides_within_window = df_cluster.copy()
                peptides_within_window["AdjustedStart"] = peptides_within_window["StartPosition"].clip(lower=cluster_start)
                peptides_within_window["AdjustedEnd"] = peptides_within_window["EndPosition"].clip(upper=cluster_end)
                
                peptides_total_size = (
                    peptides_within_window["AdjustedEnd"] - peptides_within_window["AdjustedStart"] + 1
                ).sum()
                
                cluster_total_size = 80  # Fixed size
                identity = (peptides_total_size / cluster_total_size) * 100
                
                df_cluster["ClusterTotalSize"] = cluster_total_size
                df_cluster["PeptidesTotalSize"] = peptides_total_size
                df_cluster["Identity"] = identity
                df_cluster["Putative Allergen"] = df_cluster["Identity"].apply(
                    lambda x: "putative allergen" if x >= 35 else "n/a"
                )
                output.append(df_cluster)
        
        return pd.concat(output, ignore_index=True) if output else pd.DataFrame()
    
    # Apply the clustering function to each group
    filtered_data = grouped.apply(find_clusters).reset_index(drop=True)
    
    # Save the filtered data to a new CSV file
    filtered_data.to_csv(output_csv, index=False)

def main():
    parser = argparse.ArgumentParser(description="Filter peptides within an 80 AA window, calculate cluster stats, and identify putative allergens.")
    parser.add_argument("input_csv", help="Path to the input CSV file")
    parser.add_argument("output_csv", help="Path to the output filtered CSV file")
    args = parser.parse_args()
    
    filter_peptides(args.input_csv, args.output_csv)
    print(f"Filtered peptides with cluster stats and allergen info saved to {args.output_csv}")

if __name__ == "__main__":
    main()
