#!/usr/bin/env python3
"""
Generate a self-contained HTML report for the CAAVIAR2 pipeline.
All plots are rendered with matplotlib and embedded as base64 PNGs.
"""

import sys
import os
import base64
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from io import BytesIO


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def fig_to_b64(fig):
    """Render a matplotlib figure to a base64-encoded PNG data URI."""
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches='tight', dpi=150)
    buf.seek(0)
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}"


def file_to_b64(filepath):
    """Read a PNG file and return a base64 data URI."""
    if not os.path.exists(filepath):
        return ""
    with open(filepath, "rb") as fh:
        return f"data:image/png;base64,{base64.b64encode(fh.read()).decode()}"


def safe_col(df, col):
    """Return column values filled with 0 if present, else Series of zeros."""
    if col in df.columns:
        return df[col].fillna(0)
    return pd.Series([0.0] * len(df), index=df.index)


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
        print(f"Warning: Could not parse sample order from CSV {csv_path}: {e}")
    return None

# ---------------------------------------------------------------------------
# Individual plot functions — no group logic, one bar per sample
# ---------------------------------------------------------------------------

def plot_basic_stats_panel(df, samples):
    """
    Multi-panel figure of basic statistics (inspired by R lines 56-98).
    Panels: edited_reads_pct, mean_ins_size, mean_del_size,
            mean_mh_del, mean_no_mh_del, deletion_pct, insertion_pct,
            ins1_pct, del1_pct, ins_del_ratio, aavins_pct, mh_2_plus_pct
    """
    # Derive columns that the R script computes on the fly
    if 'ins1_pct' not in df.columns and 'ins1' in df.columns and 'edited_reads' in df.columns:
        df['ins1_pct'] = df['ins1'] / df['edited_reads'] * 100
    if 'del1_pct' not in df.columns and 'del1' in df.columns and 'edited_reads' in df.columns:
        df['del1_pct'] = df['del1'] / df['edited_reads'] * 100

    metrics = [
        ('Raw_Reads',        'Raw Reads'),
        ('total_reads',      'Total Reads'),
        ('edited_reads_pct', 'Editing (%)'),
        ('mean_ins_size',    'Mean Ins Size'),
        ('mean_del_size',    'Mean Del Size'),
        ('deletion_pct',     'Deletions (%)'),
        ('insertion_pct',    'Insertions (%)'),
        ('ins1_pct',         '1bp Ins (%)'),
        ('del1_pct',         '1bp Del (%)'),
        ('ins_del_ratio',    'Ins / Del Ratio'),
        ('aavins_pct',       'AAV Ins (%)'),
        ('mh_2_plus_pct',    'MH >= 2 (%)'),
    ]
    # Keep only metrics present in df
    metrics = [(col, label) for col, label in metrics if col in df.columns]
    ncols = 4
    nrows = int(np.ceil(len(metrics) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = axes.flatten() if nrows > 1 else (axes if hasattr(axes, '__len__') else [axes])
    x = np.arange(len(samples))

    for idx, (col, label) in enumerate(metrics):
        ax = axes[idx]
        vals = safe_col(df, col)
        ax.bar(x, vals, color='#5b9bd5', edgecolor='white', width=0.7)
        ax.set_title(label, fontweight='bold', fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    # Hide unused axes
    for idx in range(len(metrics), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('Basic Statistics per Sample', fontsize=14, fontweight='bold', y=1.01)
    fig.tight_layout()
    return fig


VARIANT_COLORS = {
    'no_variant':          '#d9d9d9',
    'deletion':            '#8da0cb',
    'insertion':           '#fc8d62',
    'insertion_deletion':  '#e78ac3',
    'multiple_deletions':  '#a6d854',
    'multiple_insertions': '#ffd92f',
    'AAV_insertion':       '#e41a1c',
    'SNV':                 '#66c2a5'
}

def get_snv_col(df):
    for col in ['SNV', 'snv', 'SNVs', 'snvs']:
        if col in df.columns and safe_col(df, col).sum() > 0:
            return col
    return None

def plot_variant_type_edited_reads(df, samples):
    """
    Stacked bar chart of variant types normalised to edited reads.
    Includes SNVs if present. Excludes samples with 0 edited reads.
    """
    ed_raw = safe_col(df, 'edited_reads')
    valid_mask = ed_raw > 0
    if not valid_mask.any():
        return None

    sub_df = df[valid_mask].copy()
    sub_samples = [s for s, v in zip(samples, valid_mask) if v]

    ed = safe_col(sub_df, 'edited_reads').replace(0, np.nan)
    aavins = safe_col(sub_df, 'aavins')

    components = {}
    for vtype in ['deletion', 'insertion_deletion', 'multiple_deletions', 'multiple_insertions']:
        if vtype in sub_df.columns:
            components[vtype] = sub_df[vtype] / ed * 100
    if 'insertion' in sub_df.columns:
        components['insertion'] = (sub_df['insertion'] - aavins).clip(lower=0) / ed * 100
    components['AAV_insertion'] = aavins / ed * 100

    snv_col = get_snv_col(sub_df)
    if snv_col:
        components['SNV'] = sub_df[snv_col] / ed * 100

    if not components:
        return None

    comp_df = pd.DataFrame(components, index=sub_df.index).fillna(0)
    row_sums = comp_df.sum(axis=1).replace(0, np.nan)
    comp_norm = comp_df.div(row_sums, axis=0) * 100

    x = np.arange(len(sub_samples))
    fig, ax = plt.subplots(figsize=(max(8, len(sub_samples) * 0.8), 5))
    colors = [VARIANT_COLORS.get(col, '#999999') for col in comp_norm.columns]
    bottom = np.zeros(len(sub_samples))

    for i, col in enumerate(comp_norm.columns):
        vals = comp_norm[col].values
        ax.bar(x, vals, bottom=bottom, label=col, color=colors[i], edgecolor='white', linewidth=0.3)
        bottom += vals

    ax.set_ylabel('Edited Reads (%)')
    ax.set_title('Variant Type Composition (% of edited reads)', fontweight='bold')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(sub_samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def plot_variant_type_all_reads(df, samples):
    """
    Stacked bar chart of variant types including unedited (no variant) reads,
    normalised to total reads.
    Includes SNVs if present. Excludes samples with 0 edited reads (no_variant = 100%).
    """
    ed_raw = safe_col(df, 'edited_reads')
    valid_mask = ed_raw > 0
    if not valid_mask.any():
        return None

    sub_df = df[valid_mask].copy()
    sub_samples = [s for s, v in zip(samples, valid_mask) if v]

    tot = safe_col(sub_df, 'total_reads').replace(0, np.nan)
    aavins = safe_col(sub_df, 'aavins')

    components = {}
    if 'no_variant' in sub_df.columns:
        components['no_variant'] = sub_df['no_variant'] / tot * 100
    elif 'total_reads' in sub_df.columns and 'edited_reads' in sub_df.columns:
        components['no_variant'] = (sub_df['total_reads'] - sub_df['edited_reads']) / tot * 100

    for vtype in ['deletion', 'insertion_deletion', 'multiple_deletions', 'multiple_insertions']:
        if vtype in sub_df.columns:
            components[vtype] = sub_df[vtype] / tot * 100
    if 'insertion' in sub_df.columns:
        components['insertion'] = (sub_df['insertion'] - aavins).clip(lower=0) / tot * 100
    components['AAV_insertion'] = aavins / tot * 100

    snv_col = get_snv_col(sub_df)
    if snv_col:
        components['SNV'] = sub_df[snv_col] / tot * 100

    if not components:
        return None

    comp_df = pd.DataFrame(components, index=sub_df.index).fillna(0)
    row_sums = comp_df.sum(axis=1).replace(0, np.nan)
    comp_norm = comp_df.div(row_sums, axis=0) * 100

    x = np.arange(len(sub_samples))
    fig, ax = plt.subplots(figsize=(max(8, len(sub_samples) * 0.8), 5))
    colors = [VARIANT_COLORS.get(col, '#999999') for col in comp_norm.columns]
    bottom = np.zeros(len(sub_samples))

    for i, col in enumerate(comp_norm.columns):
        vals = comp_norm[col].values
        ax.bar(x, vals, bottom=bottom, label=col, color=colors[i], edgecolor='white', linewidth=0.3)
        bottom += vals

    ax.set_ylabel('Total Reads (%)')
    ax.set_title('Variant Type Composition (% of total reads)', fontweight='bold')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(sub_samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig



def plot_variant_subtype_stacked(df, samples):
    """
    Stacked bar chart of variant sub-types (inspired by R lines 452-511).
    Categories: D1, D2-4, D2-4+MH, D5+, D5++MH, I1, I2-4, I5+, AAV ins
    All as percentage of edited reads.
    """
    ed = safe_col(df, 'edited_reads').replace(0, np.nan)
    aavins = safe_col(df, 'aavins')

    subtypes = {}
    if 'D1' in df.columns:
        subtypes['D1'] = df['D1'] / ed * 100
    if 'D24' in df.columns:
        subtypes['D2-4'] = df['D24'] / ed * 100
    if 'D24MH' in df.columns:
        subtypes['D2-4+MH'] = df['D24MH'] / ed * 100
    if 'D5plus' in df.columns:
        subtypes['D5+'] = df['D5plus'] / ed * 100
    if 'D5plusMH' in df.columns:
        subtypes['D5++MH'] = df['D5plusMH'] / ed * 100
    if 'ins1' in df.columns:
        subtypes['I1'] = df['ins1'] / ed * 100
    if 'ins14' in df.columns and 'ins1' in df.columns:
        subtypes['I2-4'] = (df['ins14'] - df['ins1']) / ed * 100
    if 'ins_gt4' in df.columns:
        subtypes['I5+'] = (df['ins_gt4'] - aavins).clip(lower=0) / ed * 100
    subtypes['AAV ins'] = aavins / ed * 100

    if not subtypes:
        return None

    sub_df = pd.DataFrame(subtypes, index=df.index).fillna(0)
    row_sums = sub_df.sum(axis=1).replace(0, np.nan)
    sub_norm = sub_df.div(row_sums, axis=0) * 100

    x = np.arange(len(samples))
    fig, ax = plt.subplots(figsize=(max(8, len(samples) * 0.8), 5))
    colors = ['#66c2a5', '#fc8d62', '#e78ac3', '#8da0cb', '#a6d854',
              '#ffd92f', '#e5c494', '#b3b3b3', '#e41a1c']
    bottom = np.zeros(len(samples))

    for i, col in enumerate(sub_norm.columns):
        vals = sub_norm[col].values
        c = colors[i % len(colors)]
        ax.bar(x, vals, bottom=bottom, label=col, color=c, edgecolor='white', linewidth=0.3)
        bottom += vals

    ax.set_ylabel('Edited Reads (%)')
    ax.set_title('Variant Sub-type Composition (normalised)', fontweight='bold')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def plot_variant_subtype_heatmap(df, samples):
    """
    Z-scored heatmap of variant sub-types across samples
    (inspired by R lines 401-450, 547-600).
    """
    ed = safe_col(df, 'edited_reads').replace(0, np.nan)
    aavins = safe_col(df, 'aavins')

    subtypes = {}
    if 'D1' in df.columns:
        subtypes['D1'] = df['D1'] / ed * 100
    if 'D24' in df.columns:
        subtypes['D2-4'] = df['D24'] / ed * 100
    if 'D24MH' in df.columns:
        subtypes['D2-4+MH'] = df['D24MH'] / ed * 100
    if 'D5plus' in df.columns:
        subtypes['D5+'] = df['D5plus'] / ed * 100
    if 'D5plusMH' in df.columns:
        subtypes['D5++MH'] = df['D5plusMH'] / ed * 100
    if 'ins1' in df.columns:
        subtypes['I1'] = df['ins1'] / ed * 100
    if 'ins14' in df.columns and 'ins1' in df.columns:
        subtypes['I2-4'] = (df['ins14'] - df['ins1']) / ed * 100
    if 'ins_gt4' in df.columns:
        subtypes['I5+'] = (df['ins_gt4'] - aavins).clip(lower=0) / ed * 100
    subtypes['AAV ins'] = aavins / ed * 100

    if not subtypes or len(samples) < 2:
        return None

    sub_df = pd.DataFrame(subtypes, index=df.index).fillna(0)
    # Z-score per row (variant type) across samples
    mat = sub_df.values.T  # rows = variant types, cols = samples
    row_means = mat.mean(axis=1, keepdims=True)
    row_stds = mat.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1
    mat_z = (mat - row_means) / row_stds

    fig, ax = plt.subplots(figsize=(max(6, len(samples) * 0.6), max(3, len(subtypes) * 0.5)))
    im = ax.imshow(mat_z, aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
    ax.set_yticks(range(len(subtypes)))
    ax.set_yticklabels(list(subtypes.keys()), fontsize=9)
    ax.set_title('Variant Sub-types (Z-score)', fontweight='bold')
    fig.colorbar(im, ax=ax, label='Z-score', shrink=0.8)
    fig.tight_layout()
    return fig


def plot_mh_scores(df, samples):
    """
    MH combined scores bar plots (inspired by R lines 100-143).
    Shows c-scores over deletions and over edits.
    """
    score_cols = ['c25_thr2_universal', 'c25_universal',
                  'c2_thr2_universal',  'c2_universal',
                  'c14_thr2_universal', 'c14_universal',
                  'c24_thr2_universal', 'c24_universal']
    present = [c for c in score_cols if c in df.columns]
    if not present:
        return None, None

    x = np.arange(len(samples))
    # Over deletions (raw score)
    fig_del, axes_del = plt.subplots(1, len(present), figsize=(3.5 * len(present), 4), squeeze=False)
    for i, col in enumerate(present):
        ax = axes_del[0, i]
        ax.bar(x, safe_col(df, col), color='#70ad47', edgecolor='white', width=0.7)
        ax.set_title(col, fontweight='bold', fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=7)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
    fig_del.suptitle('MH Scores (per deletion)', fontsize=12, fontweight='bold', y=1.02)
    fig_del.tight_layout()

    # Over edits (score * total_deletions / edited_reads)
    fig_edit = None
    if 'total_deletions' in df.columns and 'edited_reads' in df.columns:
        ed = safe_col(df, 'edited_reads').replace(0, np.nan)
        td = safe_col(df, 'total_deletions')
        fig_edit, axes_edit = plt.subplots(1, len(present), figsize=(3.5 * len(present), 4), squeeze=False)
        for i, col in enumerate(present):
            ax = axes_edit[0, i]
            vals = safe_col(df, col) * td / ed
            ax.bar(x, vals.fillna(0), color='#4472c4', edgecolor='white', width=0.7)
            ax.set_title(col, fontweight='bold', fontsize=9)
            ax.set_xticks(x)
            ax.set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=7)
            for spine in ['top', 'right']:
                ax.spines[spine].set_visible(False)
        fig_edit.suptitle('MH Scores (normalised over edits)', fontsize=12, fontweight='bold', y=1.02)
        fig_edit.tight_layout()

    return fig_del, fig_edit


def plot_ins_del_sizes(df, samples):
    """
    Paired bar chart comparing mean insertion vs deletion size per sample.
    """
    has_ins = 'mean_ins_size' in df.columns
    has_del = 'mean_del_size' in df.columns
    if not (has_ins or has_del):
        return None

    fig, ax = plt.subplots(figsize=(max(8, len(samples) * 0.8), 5))
    x = np.arange(len(samples))
    width = 0.35

    if has_ins:
        ax.bar(x - width / 2, safe_col(df, 'mean_ins_size'), width, label='Mean Ins Size', color='#ed7d31')
    if has_del:
        ax.bar(x + width / 2, safe_col(df, 'mean_del_size'), width, label='Mean Del Size', color='#4472c4')

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
    ax.set_ylabel('Size (bp)')
    ax.set_title('Mean Insertion vs Deletion Size', fontweight='bold')
    ax.legend()
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def plot_blast_aav(df, samples):
    """
    Combined bar chart: total BLAST AAV hits and AAV integration %.
    """
    has_blast = 'BLAST_AAV' in df.columns
    has_pct = 'aavins_pct' in df.columns
    if not (has_blast or has_pct):
        return None

    x = np.arange(len(samples))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    if has_blast:
        axes[0].bar(x, safe_col(df, 'BLAST_AAV'), color='#c00000', edgecolor='white')
        axes[0].set_title('BLAST AAV Hits (total reads)', fontweight='bold')
        axes[0].set_ylabel('Count')
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
        for spine in ['top', 'right']:
            axes[0].spines[spine].set_visible(False)
    else:
        axes[0].set_visible(False)

    if has_pct:
        axes[1].bar(x, safe_col(df, 'aavins_pct'), color='#c00000', edgecolor='white')
        axes[1].set_title('AAV Integration (% of edited reads)', fontweight='bold')
        axes[1].set_ylabel('Percentage (%)')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
        for spine in ['top', 'right']:
            axes[1].spines[spine].set_visible(False)
    else:
        axes[1].set_visible(False)

    fig.tight_layout()
    return fig


def plot_diversity(df, samples):
    """Diversity metric: number of different deletion types to cover 75% of events."""
    if 'diversity_75' not in df.columns:
        return None

    x = np.arange(len(samples))
    fig, ax = plt.subplots(figsize=(max(8, len(samples) * 0.8), 5))
    ax.bar(x, safe_col(df, 'diversity_75'), color='#7030a0', edgecolor='white', width=0.7)
    ax.set_title('Deletion Diversity (types for 75% of events)', fontweight='bold')
    ax.set_ylabel('Number of Deletion Types')
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig



# ---------------------------------------------------------------------------
# Main report generation
# ---------------------------------------------------------------------------

def generate_report(summary_csv, output_html, params_json=None, samples_csv=None):
    df = pd.read_csv(summary_csv)

    # Load parameters
    params_dict = {}
    if params_json and os.path.exists(params_json):
        try:
            with open(params_json, 'r') as f:
                params_dict = json.load(f)
        except Exception as e:
            print(f"Warning: Could not read parameters JSON: {e}")

    # Resolve sample order from samples CSV if provided or from parameters
    if not samples_csv and params_dict.get('csv'):
        samples_csv = params_dict.get('csv')

    sample_order = load_sample_order(samples_csv)

    # Resolve sample name column
    if 'sample_name' not in df.columns and 'sample' in df.columns:
        df['sample_name'] = df['sample']

    if sample_order and 'sample_name' in df.columns:
        df['sample_name'] = df['sample_name'].astype(str)
        order_map = {name: idx for idx, name in enumerate(sample_order)}
        max_idx = len(sample_order) + 1000
        df['sort_key'] = df['sample_name'].map(lambda x: order_map.get(x, max_idx))
        df = df.sort_values('sort_key').drop(columns=['sort_key']).reset_index(drop=True)

    sample_names = df['sample_name'].astype(str).tolist()

    # ------------------------------------------------------------------
    # Generate all plots
    # ------------------------------------------------------------------
    plot_sections = []  # list of (title, description, b64_img)

    # 1. Basic statistics multi-panel
    fig = plot_basic_stats_panel(df, sample_names)
    plot_sections.append(('Basic Statistics',
        'Key editing metrics per sample: editing frequency, mean indel sizes, '
        'insertion/deletion percentages, AAV integration and microhomology.', fig_to_b64(fig)))
    plt.close(fig)

    # 2a. Variant type composition over all reads (including no variant)
    fig_all = plot_variant_type_all_reads(df, sample_names)
    if fig_all:
        plot_sections.append(('Variant Type Composition (% of total reads)',
            'Stacked bar chart of variant types including unedited (no variant) reads and SNVs, '
            'normalised to total reads. Samples with 0% editing (100% no variant) are excluded.', fig_to_b64(fig_all)))
        plt.close(fig_all)

    # 2b. Variant type composition over edited reads
    fig_ed = plot_variant_type_edited_reads(df, sample_names)
    if fig_ed:
        plot_sections.append(('Variant Type Composition (% of edited reads)',
            'Stacked bar chart of variant types (deletion, insertion, '
            'insertion/deletion, multiple events, AAV insertions, SNVs) normalised to edited reads. '
            'Samples with 0% editing are excluded.', fig_to_b64(fig_ed)))
        plt.close(fig_ed)


    # 3. Variant sub-type composition (stacked)
    fig = plot_variant_subtype_stacked(df, sample_names)
    if fig:
        plot_sections.append(('Variant Sub-type Composition',
            'Detailed breakdown: 1bp deletions, 2-4bp dels ± MH, 5bp+ dels ± MH, '
            '1bp insertions, 2-4bp ins, 5bp+ ins, and AAV insertions.', fig_to_b64(fig)))
        plt.close(fig)

    # 4. Variant sub-type heatmap (z-scored)
    fig = plot_variant_subtype_heatmap(df, sample_names)
    if fig:
        plot_sections.append(('Variant Sub-type Heatmap',
            'Z-scored heatmap of variant sub-types across all samples. '
            'Each row is z-scored independently to highlight relative differences.', fig_to_b64(fig)))
        plt.close(fig)


    # ------------------------------------------------------------------
    # Alignment PNGs
    # ------------------------------------------------------------------
    alignment_files = [f for f in os.listdir('.') if f.endswith('_alignments.png')]
    # Sort alignment files matching sample_names order
    def get_align_sort_key(filename):
        for idx, name in enumerate(sample_names):
            if filename.startswith(f"{name}_") or filename.startswith(f"{name}."):
                return idx
        return len(sample_names) + 100
    alignment_files.sort(key=get_align_sort_key)

    alignment_html = ""
    for f in alignment_files:
        title = f.replace('_alignments.png', '')
        img_b64 = file_to_b64(f)
        alignment_html += (
            f"<h3>{title}</h3>"
            f"<img src='{img_b64}' style='max-width:100%; border:1px solid #ddd; "
            f"border-radius:4px; padding:5px; margin-bottom:20px;'><br>"
        )

    # ------------------------------------------------------------------
    # Parameters table
    # ------------------------------------------------------------------
    params_html = ""
    if params_dict:
        params_html = "<table class='table'><tr><th>Parameter</th><th>Value</th></tr>"
        for k, v in sorted(params_dict.items()):
            params_html += f"<tr><td>{k}</td><td>{v}</td></tr>"
        params_html += "</table>"
    else:
        params_html = "<p><i>No parameters configuration found.</i></p>"

    # ------------------------------------------------------------------
    # Build plot cards HTML
    # ------------------------------------------------------------------
    plots_html = ""
    for title, desc, b64 in plot_sections:
        plots_html += f"""
    <div class="card">
        <h2>{title}</h2>
        <p>{desc}</p>
        <div class="plot-container">
            <img src="{b64}" alt="{title}">
        </div>
    </div>
"""

    # ------------------------------------------------------------------
    # Assemble HTML
    # ------------------------------------------------------------------
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>CAAVIAR2 Pipeline Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; background-color: #f9f9f9; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; border-bottom: 1px solid #eee; padding-bottom: 5px; }}
        h3 {{ color: #555; }}
        .card {{ background: white; padding: 30px; margin-bottom: 30px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
        .plot-container {{ text-align: center; margin-bottom: 20px; }}
        img {{ max-width: 100%; height: auto; }}
        table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
        th, td {{ border: 1px solid #ddd; padding: 10px; text-align: left; }}
        th {{ background-color: #f2f2f2; font-weight: bold; }}
        tr:nth-child(even) {{ background-color: #fafafa; }}
        tr:hover {{ background-color: #f1f1f1; }}
        .table-responsive {{ overflow-x: auto; }}
    </style>
</head>
<body>
    <h1>CAAVIAR2 Pipeline Report</h1>

    <div class="card">
        <h2>Analysis Parameters</h2>
        <p>Configuration parameters used for this pipeline run.</p>
        <div class="table-responsive">
            {params_html}
        </div>
    </div>

{plots_html}

    <div class="card">
        <h2>Alignment Plots</h2>
        <p>CrispRVariants alignment and frequency plots for each sample.</p>
        <div class="plot-container">
            {alignment_html if alignment_files else "<p><i>No alignment PNG files found.</i></p>"}
        </div>
    </div>

    <div class="card">
        <h2>Summary Table</h2>
        <div class="table-responsive">
            {df.to_html(classes='table', index=False, float_format='%.2f')}
        </div>
    </div>
</body>
</html>
'''
    with open(output_html, "w") as fh:
        fh.write(html_content)
    print(f"Report successfully generated: {output_html}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: generate_html_report.py <summary.csv> <output.html> [params.json] [samples.csv]")
        sys.exit(1)

    params_file = sys.argv[3] if len(sys.argv) > 3 else None
    samples_csv_file = sys.argv[4] if len(sys.argv) > 4 else None
    generate_report(sys.argv[1], sys.argv[2], params_file, samples_csv_file)

