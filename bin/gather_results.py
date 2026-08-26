#!/usr/bin/env python3
"""
Gather all merged_summary.csv files from CrispRVariants analysis results.
Handles files with different column sets by creating a unified schema.
Adds BLAST_AAV counts from external text files in 'merged_reads' directory.
"""

import pandas as pd
import os
import sys
from pathlib import Path

def load_sample_order(csv_path):
    """
    Extract sample names in order from a sample sheet CSV.
    """
    if not csv_path or not os.path.exists(csv_path):
        return None
    try:
        df_samples = pd.read_csv(csv_path)
        sample_cols = ['Run', 'sample_name', 'Sample_Name', 'sample', 'Sample', 'sample_id', 'Sample_ID', 'id', 'ID', 'name', 'Name']
        target_col = None
        for col in sample_cols:
            if col in df_samples.columns:
                target_col = col
                break
        if target_col is None and len(df_samples.columns) > 0:
            target_col = df_samples.columns[0]
            
        if target_col:
            order = []
            for val in df_samples[target_col].dropna():
                s = str(val).strip()
                if s and s not in order:
                    order.append(s)
            return order
    except Exception as e:
        print(f"    WARNING: Could not parse sample order from CSV {csv_path}: {e}")
    return None

def gather_results(samples_csv_path=None):
    """
    Gather all *_merged_summary.csv files and attach AAV counts from the current Nextflow working directory.
    Reorders rows according to samples_csv_path if provided.
    """
    data_frames = []
    
    # Find all sample-specific summary files in the current working directory
    summary_files = list(Path('.').glob('*_merged_summary.csv'))
    
    for result_file in summary_files:
        # Extract sample name from filename (e.g., "Sample1_merged_summary.csv" -> "Sample1")
        sample_name = result_file.name.replace('_merged_summary.csv', '')
        
        try:
            df = pd.read_csv(result_file)
            df.insert(0, 'sample_name', sample_name)
            
            # Look for the matching AAV count file in the same directory
            aav_file = Path('.') / f"{sample_name}_aavcount.txt"
            aav_count = 0
            if aav_file.exists():
                with open(aav_file, 'r') as f:
                    content = f.read().strip()
                    if content:
                        aav_count = int(content.split()[0])
            
            # Look for raw read count file
            raw_count_file = Path('.') / f"{sample_name}_raw_read_count.txt"
            raw_reads = None
            if raw_count_file.exists():
                with open(raw_count_file, 'r') as f:
                    content = f.read().strip()
                    if content:
                        raw_reads = int(content.split()[0])
            
            df.insert(1, 'Raw_Reads', raw_reads)
            df['BLAST_AAV'] = aav_count
            data_frames.append(df)
            
        except Exception as e:
            print(f"  ERROR reading {result_file}: {e}")

    if not data_frames:
        print("ERROR: No valid merged_summary.csv files were found!")
        sys.exit(1)
        
    combined_df = pd.concat(data_frames, ignore_index=True, sort=False)
    
    # Sort combined_df according to sample_csv if provided
    sample_order = load_sample_order(samples_csv_path)
    if sample_order and 'sample_name' in combined_df.columns:
        combined_df['sample_name'] = combined_df['sample_name'].astype(str)
        order_map = {name: idx for idx, name in enumerate(sample_order)}
        max_idx = len(sample_order) + 1000
        combined_df['sort_key'] = combined_df['sample_name'].map(lambda x: order_map.get(x, max_idx))
        combined_df = combined_df.sort_values('sort_key').drop(columns=['sort_key']).reset_index(drop=True)
        print(f"Reordered results matching samples CSV order ({len(sample_order)} samples defined).")

    return combined_df

def main():
    samples_csv = None
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if os.path.exists(arg) and not arg.endswith('_merged_summary.csv') and not arg.endswith('all_results_merged_summary.csv'):
                if arg.endswith('.csv') or os.path.isfile(arg):
                    samples_csv = arg
                    break
    
    print("=" * 60)
    print("Gathering all merged_summary.csv files")
    if samples_csv:
        print(f"Using sample sheet order from: {samples_csv}")
    print()
    
    # Gather all results
    combined_df = gather_results(samples_csv)
    
    # Save combined results (write to current directory for Nextflow)
    output_file = "all_results_merged_summary.csv"
    combined_df.to_csv(output_file, index=False)
    
    print()
    print("=" * 60)
    print("Results saved")
    print("=" * 60)
    print(f"Output file: {output_file}")
    print(f"Total rows: {len(combined_df)}")
    
    # Check if BLAST_AAV was populated
    if 'BLAST_AAV' in combined_df.columns:
        filled = combined_df['BLAST_AAV'].notna().sum()
        print(f"BLAST_AAV column populated for {filled}/{len(combined_df)} rows")
        if filled == 0:
             print("WARNING: Still 0 matches. Please check that filenames in 'merged_reads' match folder names.")
    
    print()
    print("=" * 60)
    print("DONE!")
    print("=" * 60)

if __name__ == "__main__":
    main()