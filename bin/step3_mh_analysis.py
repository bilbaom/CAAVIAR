#!/usr/bin/env python3
"""
enhanced Microhomology Analysis Script (original from Julen Torrens)
Author: bilbaom
Purpose: Calculate microhomology patterns in CRISPR-induced deletions

Usage: python mh_analysis_simple.py <all_events_del.tsv> <amplicon_sequence> <cut_site_position> <working_directory>
"""

import pandas as pd
import numpy as np
import sys
import os


def calculate_microhomology(amplicon, cut_site, deletion_pos, deletion_size):
    """
    Calculate microhomology for a single deletion event.
    
    Args:
        amplicon (str): Reference amplicon sequence
        cut_site (int): Position of the cut site (1-based)
        deletion_pos (int): Position where deletion starts (1-based)
        deletion_size (int): Size of the deletion
    
    Returns:
        dict: Dictionary containing microhomology metrics
    """
    
    # Normalize to 0-based indexing
    del_pos_0 = deletion_pos - 1
    actual_del_pos = cut_site + del_pos_0
    # crisprvariants skips position 0 so we need to add 1 to negative positions
    if deletion_pos < 0:
        actual_del_pos = actual_del_pos + 1
    
    # Empty result dictionary
    empty = {
        'mh_sequence': '',
        'mh_length': 0,
        'mh_gc_content': 0,
        'mh_location': 'none',
        # offset-aware metrics
        'mh_length_off': 0,
        'mh_sequence_off': '',
        'mh_gc_content_off': 0,
        'mh_location_off': 'none',
        'deletion_size': 0,
        'deleted_sequence': ''
    }

    # Sanity check
    if actual_del_pos + deletion_size > len(amplicon):
        return empty

    # Extract sequences
    left_seq   = amplicon[:actual_del_pos]
    deleted_seq = amplicon[actual_del_pos:actual_del_pos + deletion_size]
    right_seq  = amplicon[actual_del_pos + deletion_size:]

    # Edge case: very short flanks
    if not left_seq or not right_seq or not deleted_seq:
        return empty

    # Helper: find homology length
    def find_mh(seq1, seq2, from_end=False):
        length = min(len(seq1), len(seq2))
        for i in range(length):
            a = seq1[-(i+1)] if from_end else seq1[i]
            b = seq2[-(i+1)] if from_end else seq2[i]
            if a != b:
                return i
        return length

    left_mh  = find_mh(left_seq, deleted_seq, from_end=True)
    right_mh = find_mh(right_seq, deleted_seq, from_end=False)

    # Determine microhomology
    if left_mh > right_mh:
        mh_location, mh_len, mh_sequence = 'left', left_mh, left_seq[-left_mh:]
    elif right_mh > left_mh:
        mh_location, mh_len, mh_sequence = 'right', right_mh, right_seq[:right_mh]
    elif left_mh > 0:  # equal tie
        mh_location, mh_len, mh_sequence = 'both', left_mh, left_seq[-left_mh:]
    else:
        mh_location, mh_len, mh_sequence = 'none', 0, ''

    # GC content
    gc_content = (sum(base in "GC" for base in mh_sequence) / mh_len * 100) if mh_len else 0
    
    # --- Offset microhomology search (allow ±1 bp shift) ---
    min_off_size = 3
    off_candidates = []

    # Left offset: drop 1 bp from end of left_seq
    if len(left_seq) > 1:
        left_off = find_mh(left_seq[:-1], deleted_seq, from_end=True)
        if left_off >= min_off_size:
            seq = left_seq[-1-left_off:-1]
            gc2 = (sum(b in "GC" for b in seq) / left_off) * 100
            off_candidates.append((left_off, seq, gc2, 'left_offset'))

    # Right offset: drop 1 bp from start of right_seq
    if len(right_seq) > 1:
        right_off = find_mh(right_seq[1:], deleted_seq, from_end=False)
        if right_off >= min_off_size:
            seq = right_seq[1:1+right_off]
            gc2 = (sum(b in "GC" for b in seq) / right_off) * 100
            off_candidates.append((right_off, seq, gc2, 'right_offset'))

    # Default: strict metrics also serve as offset metrics
    mh_len_off, mh_seq_off, gc_content_off, mh_loc_off = mh_len, mh_sequence, gc_content, mh_location

    # If any offset MH is longer, replace
    if off_candidates:
        best_off = max(off_candidates, key=lambda x: x[0])
        if best_off[0] > mh_len_off:
            mh_len_off, mh_seq_off, gc_content_off, mh_loc_off = best_off
    
    # Do not allow MH for deletions of 1bp
    if deletion_size == 1:
        mh_sequence = ''
        mh_len = 0
        gc_content = 0
        mh_location = None
        mh_len_off = 0
        mh_seq_off = '',
        gc_content_off = 0,
        mh_location_off = None

    return {
        'mh_sequence': mh_sequence,
        'mh_length': mh_len,
        'mh_gc_content': gc_content,
        'mh_location': mh_location,
        # offset-aware metrics
        'mh_length_off': mh_len_off,
        'mh_sequence_off': mh_seq_off,
        'mh_gc_content_off': gc_content_off,
        'mh_location_off': mh_loc_off,
        'deletion_size': deletion_size,
        'deleted_sequence': deleted_seq
    }



def calculate_sample_statistics(df, sample_name):
    """
    Calculate various microhomology statistics for a sample.
    
    Args:
        df (pd.DataFrame): DataFrame with deletion data including MH calculations
        sample_name (str): Name of the sample to analyze
    
    Returns:
        dict: Dictionary with calculated statistics
    """
    sample_data = df[df['sample'] == sample_name].copy()
    
    if len(sample_data) == 0:
        return {}
    
    total_deletions = len(sample_data)
    
    # Basic MH counts
    mh_1 = len(sample_data[sample_data['mh_length'] == 1])
    mh_1_plus = len(sample_data[sample_data['mh_length'] >= 1])
    mh_2_plus = len(sample_data[sample_data['mh_length'] >= 2])
    mh_1_plus_off = len(sample_data[sample_data['mh_length_off'] >= 1])
    mh_2_plus_off = len(sample_data[sample_data['mh_length_off'] >= 2])
    mh_4_plus = len(sample_data[sample_data['mh_length'] >= 4])
    no_mh = len(sample_data[sample_data['mh_length'] < 2])
    no_mh_strict = len(sample_data[sample_data['mh_length'] == 0])
    
    # Percentages
    mh_2_plus_pct = (mh_2_plus / total_deletions) * 100 if total_deletions > 0 else 0
    mh_1_plus_pct = (mh_1_plus / total_deletions) * 100 if total_deletions > 0 else 0
    mh_2_plus_off_pct = (mh_2_plus_off / total_deletions) * 100 if total_deletions > 0 else 0
    mh_1_plus_off_pct = (mh_1_plus_off / total_deletions) * 100 if total_deletions > 0 else 0
    
    # Deletion size + MH combinations
    D1 = len(sample_data[sample_data['size'] == 1])
    D24 = len(sample_data[(sample_data['size'] > 1) & (sample_data['size'] < 5) & (sample_data['mh_length'] < 2)])
    D24MH = len(sample_data[(sample_data['size'] > 1) & (sample_data['size'] < 5) & (sample_data['mh_length'] >= 2)])
    D5plus = len(sample_data[(sample_data['size'] > 4) & (sample_data['mh_length'] < 2)])
    D5plusMH = len(sample_data[(sample_data['size'] > 4) & (sample_data['mh_length'] >= 2)])
    long_deletions = sample_data[sample_data['size'] > 9]
    long_del_9plus = len(long_deletions)
    
    # Calculate del means (normal arithmetic mean)
    mean_all_del = sample_data['size'].mean() if len(sample_data) > 0 else 0

    mean_mh_del = sample_data.loc[sample_data['mh_length'] >= 2, 'size'].mean() \
                     if any(sample_data['mh_length'] >= 2) else 0

    mean_no_mh_del = sample_data.loc[sample_data['mh_length'] < 2, 'size'].mean() \
                        if any(sample_data['mh_length'] < 2) else 0
    
    mean_mh_del_strict = sample_data.loc[sample_data['mh_length'] >= 1, 'size'].mean() \
                     if any(sample_data['mh_length'] >= 1) else 0

    mean_no_mh_del_strict = sample_data.loc[sample_data['mh_length'] < 1, 'size'].mean() \
                        if any(sample_data['mh_length'] < 1) else 0
    
    # MH sequence length and GC content means
    mh_seq_data = sample_data[sample_data['mh_length'] >= 2]
    mean_mh_length = mh_seq_data['mh_length'].mean() if len(mh_seq_data) > 0 else 0
    mean_mh_gc = mh_seq_data['mh_gc_content'].mean() if len(mh_seq_data) > 0 else 0
    
    # Diversity metric - number of different deletion types needed to represent 75% of events
    deletion_counts = sample_data.groupby(['pos', 'size']).size().sort_values(ascending=False)
    cumulative_pct = deletion_counts.cumsum() / deletion_counts.sum()
    diversity_75 = len(cumulative_pct[cumulative_pct <= 0.75]) + 1
    
    return {
        'sample': sample_name,
        'total_deletions': total_deletions,
        'mh_1': mh_1,
        'mh_1_plus': mh_1_plus,  # Considering MH len 1 as MH
        'mh_2_plus': mh_2_plus,
        'mh_1_plus_off': mh_1_plus_off,  # Considering MH len 1 as MH
        'mh_2_plus_off': mh_2_plus_off,
        'mh_4_plus': mh_4_plus,
        'no_mh': no_mh,
        'no_mh_strict': no_mh_strict,  # Considering MH len 1 as MH
        'mh_1_plus_pct': mh_1_plus_pct, # Considering MH len 1 as MH
        'mh_2_plus_pct': mh_2_plus_pct,
        'mh_1p_off_pct': mh_1_plus_off_pct, # Considering MH len 1 as MH
        'mh_2p_off_pct': mh_2_plus_off_pct,
        'D1': D1,
        'D24': D24,
        'D24MH': D24MH,
        'D5plus': D5plus,
        'D5plusMH': D5plusMH,
        'long_del_9plus': long_del_9plus,
        'mean_all_del': mean_all_del,
        'mean_mh_del': mean_mh_del,
        'mean_no_mh_del': mean_no_mh_del,
        'mean_mh_del_strict': mean_mh_del_strict,
        'mean_no_mh_del_strict': mean_no_mh_del_strict,
        'mean_mh_length': mean_mh_length,
        'mean_mh_gc': mean_mh_gc,
        'diversity_75': diversity_75
    }


def sigmoid_normalize(values, mean_val, std_val, steepness=1.0):
    """Apply sigmoid normalization with adjustable steepness."""
    std_val = std_val if std_val != 0 else 1e-8
    z_scores = (values - mean_val) / std_val
    return 1 / (1 + np.exp(-steepness * z_scores))


def calculate_mh_scores(df, sample_name):
    """
    Calculate microhomology-based scores for a given sample.
    
    Args:
        df (pd.DataFrame): Must contain at least columns:
            ['sample','mh_length','mh_gc_content','size']
        sample_name (str): Sample identifier
    
    Returns:
        dict: Dictionary of scores
    """
    sample_data = df[df['sample'] == sample_name].copy()
    if len(sample_data) == 0:
        return {}

    # --- Universal params ---
    universal = {
        "len": {"mean": 7.5, "std": 17},
        "mh":  {"mean": 1.32, "std": 0.66}
    }

    # --- Sample-specific params ---
    sample = {
        "len": {"mean": sample_data['size'].mean(), "std": sample_data['size'].std()},
        "mh":  {"mean": sample_data['mh_length'].mean(), "std": sample_data['mh_length'].std()},
        "gc":  {"mean": sample_data['mh_gc_content'].mean(), "std": sample_data['mh_gc_content'].std()},
    }

    # --- Normalized columns ---
    sample_data['len_norm_sample'] = sigmoid_normalize(
        sample_data['size'], sample['len']['mean'], sample['len']['std'] or 1, steepness=2.0
    )
    sample_data['len_norm_universal'] = sigmoid_normalize(
        sample_data['size'], universal['len']['mean'], universal['len']['std'], steepness=2.0
    )

    # Ensure columns exist (fill with NaN by default)
    sample_data['mh_norm_sample'] = np.nan
    sample_data['mh_norm_universal'] = np.nan
    sample_data['gc_norm_sample'] = np.nan

    mh_mask = sample_data['mh_length'] >= 1
    if mh_mask.any():
        sample_data.loc[mh_mask, 'mh_norm_sample'] = sigmoid_normalize(
            sample_data.loc[mh_mask, 'mh_length'], sample['mh']['mean'], sample['mh']['std'] or 1
        )
        sample_data.loc[mh_mask, 'mh_norm_universal'] = sigmoid_normalize(
            sample_data.loc[mh_mask, 'mh_length'], universal['mh']['mean'], universal['mh']['std']
        )

        if 'mh_gc_content' in sample_data.columns:
            gc_std = sample['gc']['std'] or 1
            sample_data.loc[mh_mask, 'gc_norm_sample'] = sigmoid_normalize(
                sample_data.loc[mh_mask, 'mh_gc_content'], sample['gc']['mean'], gc_std
            )

    # --- Scoring formulas ---
    formulas = {
        # ------------------ threshold = 1 ------------------
        "f1_sample": lambda df: df.loc[df['mh_length'] >= 1, 'len_norm_sample'].sum(),
        "f1_universal": lambda df: df.loc[df['mh_length'] >= 1, 'len_norm_universal'].sum(),

        "f2_sample": lambda df: (df.loc[df['mh_length'] >= 1, 'len_norm_sample'] *
                                 df.loc[df['mh_length'] >= 1, 'mh_norm_sample']).sum(),
        "f2_universal": lambda df: (df.loc[df['mh_length'] >= 1, 'len_norm_universal'] *
                                    df.loc[df['mh_length'] >= 1, 'mh_norm_universal']).sum(),

        "f3_sample": lambda df: (df.loc[df['mh_length'] >= 1, 'len_norm_sample'] *
                                 df.loc[df['mh_length'] >= 1, 'mh_norm_sample'] *
                                 df.loc[df['mh_length'] >= 1, 'gc_norm_sample']).sum()
                              if 'gc_norm_sample' in df else 0,

        "f4_sample": lambda df: df.loc[df['mh_length'] == 0, 'len_norm_sample'].sum(),
        "f4_universal": lambda df: df.loc[df['mh_length'] == 0, 'len_norm_universal'].sum(),

        "f5_sample": lambda df: df.loc[(df['mh_length'] == 0) & (df['size'] > 9), 'len_norm_sample'].sum(),
        "f5_universal": lambda df: df.loc[(df['mh_length'] == 0) & (df['size'] > 9), 'len_norm_universal'].sum(),

        # ------------------ threshold = 2 ------------------
        "f1_thr2_sample": lambda df: df.loc[df['mh_length'] >= 2, 'len_norm_sample'].sum(),
        "f1_thr2_universal": lambda df: df.loc[df['mh_length'] >= 2, 'len_norm_universal'].sum(),

        "f2_thr2_sample": lambda df: (df.loc[df['mh_length'] >= 2, 'len_norm_sample'] *
                                      df.loc[df['mh_length'] >= 2, 'mh_norm_sample']).sum(),
        "f2_thr2_universal": lambda df: (df.loc[df['mh_length'] >= 2, 'len_norm_universal'] *
                                         df.loc[df['mh_length'] >= 2, 'mh_norm_universal']).sum(),

        "f3_thr2_sample": lambda df: (df.loc[df['mh_length'] >= 2, 'len_norm_sample'] *
                                      df.loc[df['mh_length'] >= 2, 'mh_norm_sample'] *
                                      df.loc[df['mh_length'] >= 2, 'gc_norm_sample']).sum()
                                   if 'gc_norm_sample' in df else 0,

        "f4_thr2_sample": lambda df: df.loc[df['mh_length'] <= 1, 'len_norm_sample'].sum(),
        "f4_thr2_universal": lambda df: df.loc[df['mh_length'] <= 1, 'len_norm_universal'].sum(),

        "f5_thr2_sample": lambda df: df.loc[(df['mh_length'] <= 1) & (df['size'] > 9), 'len_norm_sample'].sum(),
        "f5_thr2_universal": lambda df: df.loc[(df['mh_length'] <= 1) & (df['size'] > 9), 'len_norm_universal'].sum(),
    }

    # --- Calculate scores ---
    scores = {name: func(sample_data) for name, func in formulas.items()}

    return scores



def calculate_combined_scores(df):
    """
    Calculate combined scores for ALL variants:
    - Both thresholds: MH ≥1 and MH ≥2  
    - Both normalization types: sample and universal
    
    This creates 36 new combined score columns total:
    - 9 combinations × 2 thresholds × 2 normalization types = 36 columns
    
    Args:
        df (pd.DataFrame): DataFrame containing microhomology scores and total_deletions
    
    Returns:
        pd.DataFrame: DataFrame with all combined score columns added
    """
    # Make a copy to avoid modifying the original DataFrame
    df_result = df.copy()
    
    # Handle division by zero
    total_deletions_safe = df['total_deletions'].replace(0, float('nan'))
    
    # Define all combinations to calculate
    normalization_types = ['sample', 'universal']
    thresholds = ['', '_thr2']  # '' for threshold≥1, '_thr2' for threshold≥2
    
    for norm_type in normalization_types:
        for threshold in thresholds:
            # Create suffix for column names
            suffix = f"{threshold}_{norm_type}"
            
            # Define score column names for this combination
            f1_col = f'f1{suffix}'
            f2_col = f'f2{suffix}'
            f3_col = f'f3{suffix}'
            f4_col = f'f4{suffix}'
            f5_col = f'f5{suffix}'
            
            # Check if required columns exist for this combination
            required_cols = [f1_col, f2_col, f4_col, f5_col]
            if not all(col in df.columns for col in required_cols):
                print(f"Warning: Missing columns for {norm_type} normalization, threshold{threshold}. Skipping.")
                continue
            
            # Create combined score column names
            combo_suffix = f"{threshold}_{norm_type}" if threshold else f"_{norm_type}"
            
            # Calculate combined scores for this variant
            df_result[f'c1{combo_suffix}'] = df[f1_col] / total_deletions_safe
            df_result[f'c2{combo_suffix}'] = df[f2_col] / total_deletions_safe
            
            # c3 only if f3 column exists
            if f3_col in df.columns:
                df_result[f'c3{combo_suffix}'] = df[f3_col] / total_deletions_safe
            else:
                df_result[f'c3{combo_suffix}'] = float('nan')
            
            # Combined formulas
            df_result[f'c14{combo_suffix}'] = (df[f1_col] + df[f4_col]) / total_deletions_safe
            df_result[f'c24{combo_suffix}'] = (df[f2_col] + df[f4_col]) / total_deletions_safe
            
            if f3_col in df.columns:
                df_result[f'c34{combo_suffix}'] = (df[f3_col] + df[f4_col]) / total_deletions_safe
            else:
                df_result[f'c34{combo_suffix}'] = float('nan')
            
            df_result[f'c15{combo_suffix}'] = (df[f1_col] + df[f5_col]) / total_deletions_safe
            df_result[f'c25{combo_suffix}'] = (df[f2_col] + df[f5_col]) / total_deletions_safe
            
            if f3_col in df.columns:
                df_result[f'c35{combo_suffix}'] = (df[f3_col] + df[f5_col]) / total_deletions_safe
            else:
                df_result[f'c35{combo_suffix}'] = float('nan')
    
    return df_result


COLUMNS_TO_DROP = [
    # insertion/deletion size percentages
    'ins_gt4_pct', 'ins_gt9_pct', 'ins_gt14_pct',
    'del_gt4_pct', 'del_gt9_pct', 'del_gt14_pct',
    # MH percentages
    'nhej_pct', 'mh_1_plus_pct', 'mh_2_plus_pct', 'mh_1p_off_pct', 'mh_2p_off_pct',
    # long deletion count
    'long_del_9plus',
    # f-scores (sample & universal, both thresholds)
    'f1_sample', 'f1_universal', 'f2_sample', 'f2_universal',
    'f3_sample', 'f4_sample', 'f4_universal', 'f5_sample', 'f5_universal',
    'f1_thr2_sample', 'f1_thr2_universal', 'f2_thr2_sample', 'f2_thr2_universal',
    'f3_thr2_sample', 'f4_thr2_sample', 'f4_thr2_universal', 'f5_thr2_sample', 'f5_thr2_universal',
    # c-scores (sample)
    'c1_sample', 'c2_sample', 'c3_sample',
    'c14_sample', 'c24_sample', 'c34_sample',
    'c15_sample', 'c25_sample', 'c35_sample',
    'c1_thr2_sample', 'c2_thr2_sample', 'c3_thr2_sample',
    'c14_thr2_sample', 'c24_thr2_sample', 'c34_thr2_sample',
    'c15_thr2_sample', 'c25_thr2_sample', 'c35_thr2_sample',
    # c-scores (universal, subset)
    'c3_universal', 'c34_universal', 'c35_universal',
    'c3_thr2_universal', 'c34_thr2_universal', 'c35_thr2_universal',
]


def drop_unwanted_columns(df):
    """Drop columns listed in COLUMNS_TO_DROP if they exist in the DataFrame."""
    cols_present = [c for c in COLUMNS_TO_DROP if c in df.columns]
    return df.drop(columns=cols_present)


def main():
    """Main function to run the microhomology analysis."""
    
    # Check command line arguments
    if len(sys.argv) != 5:
        print("Usage: python mh_analysis_simple.py <all_events_del.tsv> <amplicon_sequence> <cut_site_position> <working_directory>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    amplicon = sys.argv[2].upper()
    cut_site = int(sys.argv[3])
    work_dir = sys.argv[4]
    
    # Change to working directory
    os.chdir(work_dir)
    
    print("Loading deletion data...")
    # Load the input data
    try:
        deletions_df = pd.read_csv(input_file, sep='\t')
        print(f"Loaded {len(deletions_df)} deletion events")
    except Exception as e:
        print(f"Error loading input file: {e}")
        sys.exit(1)
    
    # Validate required columns
    required_cols = ['pos', 'size', 'type', 'sample']
    missing_cols = [col for col in required_cols if col not in deletions_df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        sys.exit(1)
    
    print("Calculating microhomology for each deletion...")
    # Calculate microhomology for each deletion
    mh_results = []
    for idx, row in deletions_df.iterrows():
        if idx % 1000 == 0:
            print(f"Processing deletion {idx+1}/{len(deletions_df)}")
        
        mh_data = calculate_microhomology(amplicon, cut_site, row['pos'], row['size'])
        
        # Combine original data with MH calculations
        result_row = {
            'pos': row['pos'],
            'size': row['size'],
            'type': row['type'],
            'sample': row['sample'],
            **mh_data
        }
        mh_results.append(result_row)
    
    # Create DataFrame with results
    results_df = pd.DataFrame(mh_results)
    
    print("Calculating sample-level statistics...")
    # Calculate statistics for each sample
    samples = results_df['sample'].unique()
    sample_stats = []
    
    for sample in samples:
        print(f"Processing sample: {sample}")
        
        # Basic statistics
        stats = calculate_sample_statistics(results_df, sample)
        
        # Normalized scores
        norm_scores = calculate_mh_scores(results_df, sample)

        # Combine all statistics
        combined_stats = {**stats, **norm_scores}
        sample_stats.append(combined_stats)
    
    # Create final results DataFrame
    final_stats_df = pd.DataFrame(sample_stats)
    
    print("Saving results...")
    # Save results — drop unwanted columns at save time only,
    # keeping final_stats_df intact for the merge step below
    results_df.to_csv('deletion_microhomology_detailed.csv', index=False)
    drop_unwanted_columns(final_stats_df).to_csv('sample_microhomology_summary.csv', index=False)
    
    print("Detailed results saved to: deletion_microhomology_detailed.csv")
    print("Sample summary saved to: sample_microhomology_summary.csv")
    
    # Print summary
    print("\nSummary:")
    print(f"Total deletions analyzed: {len(results_df)}")
    print(f"Samples processed: {len(samples)}")
    print(f"Deletions with microhomology ≥1: {len(results_df[results_df['mh_length'] >= 1])}")
    print(f"Deletions with microhomology ≥2: {len(results_df[results_df['mh_length'] >= 2])}")
    
    # Merge with step2 summary
    print("Merging with step2 summary and save to: merged_summary.csv")
    try:
        step2_df = pd.read_csv("summary_df.tsv", sep="\t")
        # Merge using the full final_stats_df (f-scores still present, needed by calculate_combined_scores)
        merged_df = pd.merge(step2_df, final_stats_df, on="sample", how="outer")
        
        # Calculate combined scores (requires f-score columns)
        merged_df_c = calculate_combined_scores(merged_df)
        # Drop unwanted columns before saving
        merged_df_c = drop_unwanted_columns(merged_df_c)
        # Save results csv
        merged_df_c.to_csv("merged_summary.csv", index=False)
    
    except Exception as e:
        print(f"Warning: Could not merge with step2 summary: {e}")
    
    print("Analysis complete!")


if __name__ == "__main__":
    main()