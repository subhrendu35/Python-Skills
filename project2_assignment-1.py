# Project 2
# Completed by Subhrendu Gupta
# University ID: 14359440

# This project loads a CSV, 
# Applies QC filters, prints the summaries,
# Saves a TSV file
# Also optionally generates figures


import argparse
import pandas as pd
import matplotlib.pyplot as plt

# Using Parse user arguments

def parse_arguments():
    """ Collect all command-line arguments."""
    parser = argparse.ArgumentParser(
        description="RNA-disease associated QC pipeline."
    )
    
    
    parser.add_argument("--input_CSV", required=True,
                        help="Path to the input CSV file.")
    parser.add_argument("--output_CSV", required=True,
                        help="Path to save the filtered TSV file.")
    parser.add_argument("--rna_thr", type=float, required=True,
                        help="RNA-locus confidence threshold (strict>).")
    parser.add_argument("--gene_thr", type=float, required=True,
                        help="Gene-disease confidence threshold (strict>).")
    parser.add_argument("--plots", action="store_true",
                        help="Include this flag to generate figures.")
                        
    
    return parser.parse_args()
                        
                        
# Applying QC Filters
                        
def quality_control(df, thr_rna, thr_gene):
    """Filter rows that are meeting both the RNA and gene confidence thresholds."""
   
    df_no_na = df.dropna(subset=["rna2locus_conf_score", "gene2disease_conf_score"])

    filtered = df_no_na[
        (df_no_na["rna2locus_conf_score"] > thr_rna) &
        (df_no_na["gene2disease_conf_score"] > thr_gene)
    ]
    
    return filtered.copy()
                        
                        
# Displaying summary 
                        
def show_summary(df, label):
    """Print statistics for a dataframe."""
    print(f"\n---- Summary ({label}) -----")
    # Using 'count' along with the specified metrics for a more complete picture
    print(df[["rna2locus_conf_score", "gene2disease_conf_score", "NofPmids", "NofSnps"]].describe())
                        
                        
# Generating plots 
                        
def generate_plots(before_df, after_df):
    """Generate all the important figures."""
                        
    # 1. Top 5 organisms after QC
    plt.figure(figsize=(7,5))
    after_df["Organism"].value_counts().head(5).plot(kind="bar")
    plt.title("Top 5 Organisms After QC")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("Organisms_after_qc.png")
    plt.close()
    print("Saved 'Organisms_after_qc.png'")
                        
    # Preparing the metric list
    metrics = ["mean", "median", "min", "max"]
                        
                        
    # 2. RNA-locus score comparison
    rna_stats = pd.DataFrame({
        "Before QC": before_df["rna2locus_conf_score"].agg(metrics),
        "After QC": after_df["rna2locus_conf_score"].agg(metrics),
    })
                        
    rna_stats.plot(kind="bar", rot=0)
    plt.title("RNA-Locus Confidence Score (Before vs After QC)")
    plt.ylabel("Score")
    plt.tight_layout()
    plt.savefig("rna2locus_comparison.png")
    plt.close()
    print("Saved 'rna2locus_comparison.png'")
                        
                        
                        
    # 3. Gene-disease score comparison
    gene_stats = pd.DataFrame({
        "Before QC": before_df["gene2disease_conf_score"].agg(metrics),
        "After QC": after_df["gene2disease_conf_score"].agg(metrics),    
    })
            
                        
    gene_stats.plot(kind="bar", rot=0)
    plt.title("Gene-Disease Confidence Score (Before vs After QC)")
    plt.ylabel("Score")
    plt.tight_layout()
    plt.savefig("gene2disease_comparison.png")
    plt.close()
    print("Saved 'gene2disease_comparison.png'")
                        
                        
                        # Main workflow

def run():
    """Main workflow controller for the QC script."""    
                        
    # Collect arguments
    args = parse_arguments()
    
    # Print settings (Corrected to use consistent attribute names from args object)
    print("\n--- Settings ---")
    print("Input file:", args.input_CSV)
    print("Output file:", args.output_CSV)
    print("RNA threshold:", args.rna_thr)
    print("Gene threshold:", args.gene_thr)
    print("Generate plots:", args.plots)
                        
    # Loading the  input CSV            
    df_raw = pd.read_csv(args.input_CSV)
    print("\nLoaded dataset with", len(df_raw), "rows.")
                        
    # Showing the pre-QC summary
    show_summary(df_raw, "Before QC")
                        
    # Applying the  QC filtering
    df_clean = quality_control(df_raw, args.rna_thr, args.gene_thr)
    print("\nThe Filtered dataset contains", len(df_clean), "rows.")
                        
    # Showing the  post-QC summary
    show_summary(df_clean, "After QC")
                        
    # Top 3 organisms after QC
    print("\n--- Top 3 Organisms After QC---")
    print(df_clean["Organism"].value_counts().head(3))
                        
    # Save TSV output (Corrected to use consistent attribute name output_CSV)
    df_clean.to_csv(args.output_CSV, sep="\t", index=False)
    print(f"\nFiltered data saved to: {args.output_CSV}")
                        
    # Generating plots
    if args.plots:
        print("\nGenerating figures...")
        generate_plots(df_raw, df_clean)
        print("All the photos are generated")


# Calling main workflow
if __name__ == "__main__":
    run()