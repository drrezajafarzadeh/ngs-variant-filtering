import pandas as pd
import argparse
import os

def filter_variants(input_path, output_path, max_af=0.01, min_cadd=20.0):
    """
    Filters NGS variants based on allele frequency, pathogenicity predictions,
    and ClinVar classification to prioritize candidates for clinical review.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"[*] Loading variants from {input_path}...")
    df = pd.read_csv(input_path, sep='\t')

    initial_count = len(df)
    print(f"[*] Total variants loaded: {initial_count}")

    # 1. Filter by Allele Frequency (rare variants)
    # Handle cases where AF is missing (NaN) - often true for novel pathogenic variants
    rare_mask = (df['gnomAD_AF'].isna()) | (df['gnomAD_AF'] <= max_af)

    # 2. Filter by pathogenicity/impact
    # Prioritize Pathogenic/Likely Pathogenic in ClinVar OR high computational prediction (CADD)
    pathogenic_clinvar = df['ClinVar_Significance'].str.contains('Pathogenic|Likely Pathogenic', case=False, na=False)
    high_cadd = df['CADD_Phred'] >= min_cadd
    deleterious_consequence = df['Consequence'].str.contains('stop_gained|frameshift|missense|splice_acceptor|splice_donor', case=False, na=False)

    # Combined criteria: Must be rare AND (ClinVar Pathogenic OR (High CADD AND deleterious consequence))
    filtered_df = df[rare_mask & (pathogenic_clinvar | (high_cadd & deleterious_consequence))]

    final_count = len(filtered_df)
    print(f"[+] Filtering complete. Prioritized {final_count} / {initial_count} variants.")

    # Sort by clinical priority: Pathogenic first, then by CADD score descending
    filtered_df = filtered_df.sort_values(
        by=['ClinVar_Significance', 'CADD_Phred'],
        ascending=[False, False]
    )

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    filtered_df.to_csv(output_path, sep='\t', index=False)
    print(f"[+] Prioritized variants saved to {output_path}")

if __name__ == "__main__":
    # Default paths for demonstration
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_input = os.path.join(script_dir, '../data/example_variants.tsv')
    default_output = os.path.join(script_dir, '../data/prioritized_variants.tsv')

    filter_variants(default_input, default_output)
