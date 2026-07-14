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

    for idx, (col, label) in enumerate(metrics):
        ax = axes[idx]
        vals = safe_col(df, col)
        ax.bar(samples, vals, color='#5b9bd5', edgecolor='white', width=0.7)
        ax.set_title(label, fontweight='bold', fontsize=11)
        ax.tick_params(axis='x', rotation=45, labelsize=8)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    # Hide unused axes
    for idx in range(len(metrics), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('Basic Statistics per Sample', fontsize=14, fontweight='bold', y=1.01)
    fig.tight_layout()
    return fig


def plot_variant_type_stacked(df, samples):
    """
    Stacked bar chart of variant types normalised to edited reads
    (inspired by R lines 317-365).
    """
    # Compute component percentages over edited reads
    ed = safe_col(df, 'edited_reads').replace(0, np.nan)
    aavins = safe_col(df, 'aavins')

    components = {}
    for vtype in ['deletion', 'insertion_deletion', 'multiple_deletions', 'multiple_insertions']:
        if vtype in df.columns:
            components[vtype] = df[vtype] / ed * 100
    # Insertion minus AAV
    if 'insertion' in df.columns:
        components['insertion'] = (df['insertion'] - aavins) / ed * 100
    components['AAV_insertion'] = aavins / ed * 100

    if not components:
        return None

    comp_df = pd.DataFrame(components, index=df.index).fillna(0)
    # Normalise each sample to 100%
    row_sums = comp_df.sum(axis=1).replace(0, np.nan)
    comp_norm = comp_df.div(row_sums, axis=0) * 100

    fig, ax = plt.subplots(figsize=(max(8, len(samples) * 0.8), 5))
    cmap = plt.colormaps.get('Set2', plt.cm.Set2)
    colors = [cmap(i) for i in range(len(comp_norm.columns))]
    bottom = np.zeros(len(samples))

    for i, col in enumerate(comp_norm.columns):
        vals = comp_norm[col].values
        ax.bar(samples, vals, bottom=bottom, label=col, color=colors[i], edgecolor='white', linewidth=0.3)
        bottom += vals

    ax.set_ylabel('Edited Reads (%)')
    ax.set_title('Variant Type Composition (normalised)', fontweight='bold')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
    ax.tick_params(axis='x', rotation=45, labelsize=8)
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

    fig, ax = plt.subplots(figsize=(max(8, len(samples) * 0.8), 5))
    colors = ['#66c2a5', '#fc8d62', '#e78ac3', '#8da0cb', '#a6d854',
              '#ffd92f', '#e5c494', '#b3b3b3', '#e41a1c']
    bottom = np.zeros(len(samples))

    for i, col in enumerate(sub_norm.columns):
        vals = sub_norm[col].values
        c = colors[i % len(colors)]
        ax.bar(samples, vals, bottom=bottom, label=col, color=c, edgecolor='white', linewidth=0.3)
        bottom += vals

    ax.set_ylabel('Edited Reads (%)')
    ax.set_title('Variant Sub-type Composition (normalised)', fontweight='bold')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
    ax.tick_params(axis='x', rotation=45, labelsize=8)
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
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
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

    # Over deletions (raw score)
    fig_del, axes_del = plt.subplots(1, len(present), figsize=(3.5 * len(present), 4), squeeze=False)
    for i, col in enumerate(present):
        ax = axes_del[0, i]
        ax.bar(samples, safe_col(df, col), color='#70ad47', edgecolor='white', width=0.7)
        ax.set_title(col, fontweight='bold', fontsize=9)
        ax.tick_params(axis='x', rotation=45, labelsize=7)
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
            ax.bar(samples, vals.fillna(0), color='#4472c4', edgecolor='white', width=0.7)
            ax.set_title(col, fontweight='bold', fontsize=9)
            ax.tick_params(axis='x', rotation=45, labelsize=7)
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
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
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

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    if has_blast:
        axes[0].bar(samples, safe_col(df, 'BLAST_AAV'), color='#c00000', edgecolor='white')
        axes[0].set_title('BLAST AAV Hits (total reads)', fontweight='bold')
        axes[0].set_ylabel('Count')
        axes[0].tick_params(axis='x', rotation=45, labelsize=8)
        for spine in ['top', 'right']:
            axes[0].spines[spine].set_visible(False)
    else:
        axes[0].set_visible(False)

    if has_pct:
        axes[1].bar(samples, safe_col(df, 'aavins_pct'), color='#c00000', edgecolor='white')
        axes[1].set_title('AAV Integration (% of edited reads)', fontweight='bold')
        axes[1].set_ylabel('Percentage (%)')
        axes[1].tick_params(axis='x', rotation=45, labelsize=8)
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

    fig, ax = plt.subplots(figsize=(max(8, len(samples) * 0.8), 5))
    ax.bar(samples, safe_col(df, 'diversity_75'), color='#7030a0', edgecolor='white', width=0.7)
    ax.set_title('Deletion Diversity (types for 75% of events)', fontweight='bold')
    ax.set_ylabel('Number of Deletion Types')
    ax.tick_params(axis='x', rotation=45, labelsize=8)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Main report generation
# ---------------------------------------------------------------------------

def generate_report(summary_csv, output_html, params_json=None):
    df = pd.read_csv(summary_csv)

    # Load parameters
    params_dict = {}
    if params_json and os.path.exists(params_json):
        try:
            with open(params_json, 'r') as f:
                params_dict = json.load(f)
        except Exception as e:
            print(f"Warning: Could not read parameters JSON: {e}")

    # Resolve sample name column
    if 'sample_name' not in df.columns and 'sample' in df.columns:
        df['sample_name'] = df['sample']
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

    # 2. Variant type composition (stacked)
    fig = plot_variant_type_stacked(df, sample_names)
    if fig:
        plot_sections.append(('Variant Type Composition',
            'Stacked bar chart of variant types (deletion, insertion, '
            'insertion/deletion, multiple events, AAV insertions) normalised to edited reads.', fig_to_b64(fig)))
        plt.close(fig)

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
    alignment_files = sorted([f for f in os.listdir('.') if f.endswith('_alignments.png')])
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
        print("Usage: generate_html_report.py <summary.csv> <output.html> [params.json]")
        sys.exit(1)

    params_file = sys.argv[3] if len(sys.argv) > 3 else None
    generate_report(sys.argv[1], sys.argv[2], params_file)
