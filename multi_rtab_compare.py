#!/usr/bin/env python3
"""
===============================================================================
  MULTI-RTAB GWAS UNITIG COMPARISON PIPELINE
  Three-stage workflow for 3+ .Rtab files:
    Stage 1 - Load & Visualise each file individually
    Stage 2 - Cross-file comparison visualisations
    Stage 3 - MSA of shared unitigs (Biopython + FASTA export)
===============================================================================

Usage:
    python multi_rtab_compare.py file1.Rtab file2.Rtab file3.Rtab

    Optionally name your files:
    python multi_rtab_compare.py azm_sr.Rtab cip_sr.Rtab tet_sr.Rtab \
           --labels AZM CIP TET --outdir results

Dependencies:
    pip install pandas numpy scipy matplotlib seaborn biopython scikit-learn
"""

import argparse
import os
import sys
import json
import warnings
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist, squareform
from scipy import stats
from sklearn.decomposition import PCA

# Biopython for MSA
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio import AlignIO

warnings.filterwarnings("ignore")
np.random.seed(42)

# ===== METHODOLOGY REFERENCE =====
# MSA (Multiple Sequence Alignment):
#   - Identifies unitigs (DNA sequences) shared across 2+ .Rtab files
#   - Uses Biopython PairwiseAligner (global mode) for pairwise comparison
#   - Alignment scoring: match=+2, mismatch=-1, open_gap=-2, extend_gap=-0.5
#   - Calculates sequence identity: (alignment_score) / (2 * max_sequence_length)
#   - Outputs: FASTA files, pairwise identity matrix, summary statistics
#
# PCA (Principal Component Analysis):
#   - Dimensionality reduction technique that finds the directions of maximum variance
#   - Transforms high-dimensional sample data into uncorrelated principal components
#   - In this pipeline: samples are columns, unitigs are rows
#   - Result: identifies which samples are most similar or different
#
# SCREE PLOT:
#   - Bar chart showing variance explained by each principal component
#   - Cumulative line shows total variance with increasing number of components
#   - Used to determine how many components are needed (usually 80-90% threshold)
#   - Helps identify if data has clear structure or if noise dominates
#
# ===== COLOUR SCHEME =====
PAL = {
    "teal":   "#2EC4B6",
    "red":    "#E63946",
    "blue":   "#457B9D",
    "navy":   "#1D3557",
    "gold":   "#F4A261",
    "green":  "#2A9D8F",
    "purple": "#9B5DE5",
    "orange": "#F77F00",
    "pink":   "#E07AB1",
    "bg":     "#F8F9FA",
    "grid":   "#E0E0E0",
    "muted":  "#6C757D",
}
FILE_COLORS = [PAL["teal"], PAL["red"], PAL["blue"],
               PAL["gold"], PAL["purple"], PAL["orange"]]

plt.rcParams.update({
    "figure.facecolor": PAL["bg"],
    "axes.facecolor":   "#FFFFFF",
    "axes.edgecolor":   PAL["grid"],
    "axes.linewidth":   0.8,
    "axes.labelcolor":  PAL["navy"],
    "axes.labelweight": "bold",
    "axes.labelsize":   10,
    "xtick.color":      PAL["navy"],
    "ytick.color":      PAL["navy"],
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "grid.color":       PAL["grid"],
    "grid.linestyle":   "--",
    "grid.linewidth":   0.5,
    "font.family":      "DejaVu Sans",
    "text.color":       PAL["navy"],
})


# ==================================================================================
#  HELPERS
# ==================================================================================

def savefig(fig, path, dpi=150):
    """Save figure with consistent settings and close."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=PAL["bg"])
    plt.close(fig)
    print(f"  [OK]  {Path(path).name}")


def gc(seq: str) -> float:
    """GC content of a DNA string."""
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / max(len(seq), 1)


def short_label(seq: str, n=30) -> str:
    """Truncate a unitig sequence for display."""
    return seq[:n] + "..." if len(seq) > n else seq


# ==================================================================================
#  STAGE 1 - LOAD & PER-FILE VISUALISATION
# ==================================================================================

def load_rtab(filepath: str, label: str) -> dict:
    """
    Load and validate binary unitig presence/absence matrix from .Rtab file.

    Performs quality control: removes duplicates, all-absent unitigs, empty
    samples. Coerces to binary (0/1) and classifies unitigs as core
    (>95% prevalence), accessory (5-95%), or rare (<5%).

    Parameters
    ----------
    filepath : str
        Path to .Rtab file (space-separated, first column = sequence ID)
    label : str
        Short label for this dataset (used in plots, e.g., "AZM", "CFX")

    Returns
    -------
    dict
        Keys: df (binary DataFrame), freq (prevalence per unitig),
        sample_load (unitigs per sample), classes (core/accessory/rare
        classification), label, filepath. All metadata needed for
        downstream visualization and analysis.

    Raises
    ------
    FileNotFoundError
        If filepath does not exist
    ValueError
        If file has <1 row or <2 columns (needs sequences and samples)
    RuntimeError
        If file becomes empty after QC filtering

    Notes
    -----
    Pan-genome classification:
      - Core (>95%):      Essential markers, present in all/nearly-all samples
      - Accessory (5-95%): Common but variable, typical GWAS loci
      - Rare (<5%):       Population-specific, may be recent mutations
    """
    fp = Path(filepath)
    
    # Validation: file exists
    if not fp.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    print(f"\n  Loading [{label}] {fp.name}")
    
    try:
        df = pd.read_csv(fp, sep=" ", index_col=0)
    except Exception as e:
        raise ValueError(f"Failed to parse {fp.name} as space-separated file: {e}")
    
    # Validation: dimensions
    if df.shape[0] < 1 or df.shape[1] < 1:
        raise ValueError(f"File {fp.name} is empty (expected sequences x samples matrix)")
    
    # Coerce to binary
    if not set(np.unique(df.values)).issubset({0, 1}):
        df = (df > 0).astype(int)
        print(f"    [INFO] Coerced to binary (values > 0 marked as present)")

    # Remove duplicate unitig sequences
    n_dup = df.index.duplicated().sum()
    if n_dup > 0:
        df = df[~df.index.duplicated(keep='first')]
        print(f"    [WARN] Dropped {n_dup} duplicate unitigs (kept first occurrence)")

    # Remove all-absent unitigs (all zeros)
    row_sums = df.sum(axis=1)
    n_absent = (row_sums == 0).sum()
    if n_absent > 0:
        df = df[row_sums > 0]
        print(f"    [INFO] Removed {n_absent} all-absent unitigs")

    # Remove empty samples (all zeros)
    col_sums = df.sum(axis=0)
    n_empty = (col_sums == 0).sum()
    if n_empty > 0:
        df = df.loc[:, col_sums > 0]
        print(f"    [INFO] Removed {n_empty} empty samples")

    # Final validation
    if df.shape[0] < 1 or df.shape[1] < 1:
        raise RuntimeError(f"File {fp.name} became empty after QC filtering. "
                          f"Check data validity (all zeros?)")

    freq        = df.sum(axis=1) / df.shape[1]
    sample_load = df.sum(axis=0)

    classes = pd.Series("accessory", index=freq.index)
    classes[freq > 0.95] = "core"
    classes[freq < 0.05] = "rare"
    
    n_core = (classes == 'core').sum()
    n_accessory = (classes == 'accessory').sum()
    n_rare = (classes == 'rare').sum()

    print(f"    [OK]  {df.shape[0]:,} unitigs  x  {df.shape[1]:,} samples")
    print(f"          Core={n_core}  Accessory={n_accessory}  Rare={n_rare}")

    return dict(df=df, freq=freq, sample_load=sample_load,
                classes=classes, label=label, filepath=str(fp))


def visualise_single(data: dict, outdir: Path):
    """
    Generate 6 per-file visualizations for unitig distribution and patterns.

    Creates PNG plots showing frequency distribution, pan-genome structure,
    sample diversity, prevalence ranking, sample relationships (PCA), and
    hierarchical clustering heatmap for top variable unitigs.

    Parameters
    ----------
    data : dict
        Output from load_rtab(): df, freq, sample_load, classes, label, filepath
    outdir : Path
        Output directory for PNG files

    Output Files
    -------
    {label}_1_freq_hist.png
        Histogram of unitig prevalence with core/accessory/rare zones
    {label}_2_partition.png
        Horizontal stacked bar: core/accessory/rare unitig counts
    {label}_3_sample_load.png
        Distribution of unitigs present per sample (with mean + 5th percentile)
    {label}_4_prevalence_rank.png
        Ranked curve showing prevalence by unitig rank (log scale typical)
    {label}_5_pca.png
        Left: PC1 vs PC2 scatter (sample relationships)
        Right: Scree plot (variance explained by component)
    {label}_6_heatmap.png
        Heatmap of top-50 most variable unitigs, samples hierarchically
        clustered by Jaccard distance

    Notes
    -----
    All plots use consistent color scheme; sample clustering uses
    hierarchical agglomerative clustering with Jaccard distance.
    """
    label     = data["label"]
    freq      = data["freq"].values
    load      = data["sample_load"].values
    classes   = data["classes"]
    df        = data["df"]
    lbl       = label  # short alias

    class_counts = classes.value_counts()
    core_n  = class_counts.get("core",      0)
    acc_n   = class_counts.get("accessory", 0)
    rare_n  = class_counts.get("rare",      0)
    total   = len(freq)

    print(f"\n  [VIS] {label}")

    # 1  Frequency histogram -
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(freq, bins=40, color=PAL["teal"], alpha=0.75,
            edgecolor="white", lw=0.6, zorder=3)
    for x0, x1, c, n, name in [
        (0,    0.05, PAL["red"],  rare_n,  "Rare"),
        (0.05, 0.95, PAL["teal"], acc_n,   "Accessory"),
        (0.95, 1.0,  PAL["blue"], core_n,  "Core"),
    ]:
        ax.axvspan(x0, x1, alpha=0.07, color=c)
        ax.text((x0 + x1) / 2, ax.get_ylim()[1] * 0.88,
                f"{name}\n{n}", ha="center", va="top",
                fontsize=8, fontweight="bold", color=c)
    ax.axvline(0.05, color=PAL["red"],  lw=1.6, ls="--")
    ax.axvline(0.95, color=PAL["blue"], lw=1.6, ls="--")
    ax.set_xlabel("Presence frequency")
    ax.set_ylabel("Unitig count")
    ax.set_title(f"{lbl}  -  Unitig Frequency Distribution\n"
                 f"{total} unitigs · {df.shape[1]:,} samples",
                 fontweight="bold", pad=10)
    ax.yaxis.grid(True, alpha=0.4); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / f"{label}_1_freq_hist.png")

    # 2  Partition horizontal bar -
    fig, ax = plt.subplots(figsize=(8, 2.8))
    bars  = [rare_n, acc_n, core_n]
    names = ["Rare (<5%)", "Accessory (5-95%)", "Core (>95%)"]
    cols  = [PAL["red"], PAL["teal"], PAL["blue"]]
    left  = 0
    for v, n, c in zip(bars, names, cols):
        ax.barh(0, v, left=left, color=c, edgecolor="white", lw=1.5, height=0.5)
        if v > 0:
            ax.text(left + v / 2, 0, f"{n}\n{v} ({v/total*100:.1f}%)",
                    ha="center", va="center", fontsize=8.5,
                    fontweight="bold", color="white")
        left += v
    ax.set_xlim(0, total); ax.set_yticks([])
    ax.set_xlabel("Unitig count"); ax.set_title(
        f"{lbl}  -  Pan-genome Partition ({total} unitigs)",
        fontweight="bold", pad=10)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / f"{label}_2_partition.png")

    # 3  Per-sample load histogram -
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(load, bins=50, color=PAL["blue"], alpha=0.8,
            edgecolor="white", lw=0.6, zorder=3)
    ax.axvline(load.mean(), color=PAL["red"],  lw=2, ls="-.",
               label=f"Mean = {load.mean():.0f}")
    ax.axvline(np.percentile(load, 5), color=PAL["gold"], lw=1.5, ls=":",
               label=f"5th pct = {np.percentile(load,5):.0f}")
    ax.set_xlabel("Unitigs present per sample")
    ax.set_ylabel("Sample count")
    ax.set_title(f"{lbl}  -  Per-Sample Unitig Load\n"
                 f"Mean={load.mean():.0f}  SD={load.std():.1f}  "
                 f"Min={load.min()}  Max={load.max()}",
                 fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.yaxis.grid(True, alpha=0.4); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / f"{label}_3_sample_load.png")

    # 4  Ranked prevalence curve -
    sfq  = np.sort(freq)[::-1] * 100
    rnks = np.arange(1, len(sfq) + 1)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.fill_between(rnks, sfq, alpha=0.15, color=PAL["blue"])
    ax.plot(rnks, sfq, color=PAL["blue"], lw=2.2)
    ax.axhline(5,  color=PAL["red"],  lw=1.6, ls="--", label="Rare  <5%")
    ax.axhline(95, color=PAL["teal"], lw=1.6, ls="--", label="Core >95%")
    ax.axvspan(0, (sfq > 95).sum(), alpha=0.07, color=PAL["teal"])
    ax.axvspan(len(sfq) - (sfq < 5).sum(), len(sfq), alpha=0.07, color=PAL["red"])
    ax.set_xlim(0, len(rnks)); ax.set_ylim(-2, 105)
    ax.set_xlabel("Unitig rank (high -> low)")
    ax.set_ylabel("Prevalence (%)")
    ax.set_title(f"{lbl}  -  Ranked Unitig Prevalence",
                 fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.yaxis.grid(True, alpha=0.4); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / f"{label}_4_prevalence_rank.png")

    # 5  PCA (Sample space analysis) -
    # PCA = Principal Component Analysis: Dimensionality reduction via Singular Value Decomposition
    # Purpose: Transform high-dimensional binary unitig data into 2-3 main axes of variation
    # Each PC captures max variance in that direction; samples similar in PCA space are genetically close
    try:
        mat   = df.T.values.astype(float)
        n_pcs = min(10, mat.shape[0], mat.shape[1])
        pca   = PCA(n_components=n_pcs, random_state=42)
        coords = pca.fit_transform(mat)
        ev     = pca.explained_variance_ratio_ * 100

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle(f"{lbl}  -  PCA of Sample Space",
                     fontweight="bold", fontsize=11)

        # LEFT PANEL: PC1 vs PC2 scatter (shows sample clustering in 2D space)
        # Each dot = one sample, positioned by its first 2 principal components
        ax1.scatter(coords[:, 0], coords[:, 1], s=8, alpha=0.4,
                    color=PAL["teal"], edgecolors="none")
        ax1.set_xlabel(f"PC1  ({ev[0]:.1f}%)")
        ax1.set_ylabel(f"PC2  ({ev[1]:.1f}%)")
        ax1.set_title("PC1 vs PC2", fontweight="bold")
        ax1.spines["top"].set_visible(False); ax1.spines["right"].set_visible(False)

        # RIGHT PANEL: Scree plot (shows variance explained by each component)
        # Bars = individual variance per PC, line = cumulative variance
        # Goal: find 'elbow' where adding more PCs doesn't help much
        x_pcs = np.arange(1, len(ev) + 1)
        ax2.bar(x_pcs, ev, color=PAL["blue"], alpha=0.8, edgecolor="white")
        ax2.plot(x_pcs, np.cumsum(ev), "o-", color=PAL["red"],
                 lw=1.8, ms=5, label="Cumulative %")
        ax2.set_xlabel("Principal component")
        ax2.set_ylabel("Variance explained (%)")
        ax2.set_title("Scree plot", fontweight="bold")
        ax2.legend(frameon=False, fontsize=9)
        ax2.yaxis.grid(True, alpha=0.4); ax2.set_axisbelow(True)
        ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)

        plt.tight_layout()
        savefig(fig, outdir / f"{label}_5_pca.png")
    except Exception as e:
        print(f"    [WARN]  PCA skipped: {e}")

    # - 6  Heatmap -
    try:
        # Top 50 by variance = p*(1-p)
        var     = data["freq"] * (1 - data["freq"])
        top_idx = var.nlargest(50).index
        n_sub   = min(150, df.shape[1])
        sub_col = np.random.choice(df.columns, n_sub, replace=False)
        sub_df  = df.loc[top_idx, sub_col]

        row_lnk = linkage(pdist(sub_df.values,   metric="jaccard"), method="ward")
        col_lnk = linkage(pdist(sub_df.T.values, metric="jaccard"), method="ward")
        row_ord = leaves_list(row_lnk)
        col_ord = leaves_list(col_lnk)
        sdf     = sub_df.iloc[row_ord].iloc[:, col_ord]
        ylbls   = [short_label(u, 28) for u in sdf.index]

        fig, ax = plt.subplots(figsize=(14, 9))
        cmap    = LinearSegmentedColormap.from_list("amr", ["#FFFFFF", PAL["blue"]])
        im      = ax.imshow(sdf.values, aspect="auto", cmap=cmap,
                            interpolation="nearest", vmin=0, vmax=1)
        plt.colorbar(im, ax=ax, fraction=0.018, pad=0.01,
                     label="Present(1) / Absent(0)")
        ax.set_yticks(range(len(ylbls)))
        ax.set_yticklabels(ylbls, fontsize=5.5, fontfamily="monospace")
        ax.set_xticks([])
        ax.set_xlabel(f"{n_sub} samples (Jaccard-clustered)")
        ax.set_ylabel("Top-50 variable unitigs (clustered)")
        ax.set_title(f"{lbl}  -  Clustered Heatmap",
                     fontweight="bold", pad=10)
        plt.tight_layout()
        savefig(fig, outdir / f"{label}_6_heatmap.png")
    except Exception as e:
        print(f"    [WARN]  Heatmap skipped: {e}")


# ==================================================================================
#  STAGE 2 - CROSS-FILE COMPARISON VISUALISATIONS
# ==================================================================================

def compare_visualise(datasets: list, outdir: Path):
    """
    Generate 5 cross-file comparison visualizations.

    Creates plots showing pan-genome structure consistency, sample-level
    diversity patterns, frequency distributions, pairwise overlaps, and
    unitig sharing patterns across all input files.

    Parameters
    ----------
    datasets : list of dict
        Output from load_rtab() for each file

    outdir : Path
        Output directory for PNG files

    Output Files
    -------
    cmp_A_partition.png
        Stacked horizontal bar: core/accessory/rare counts per file
        (reveals pan-genome structure consistency across datasets)

    cmp_B_load_violin.png
        Violin plot: distribution of unitigs per sample, grouped by file
        (shows sample diversity and minimum thresholds)

    cmp_C_freq_overlay.png
        Overlaid histograms: unitig frequency distribution for all files
        (highlights frequency biases and marker selection effects)

    cmp_D_pairwise_shared.png
        Heatmap: count of shared unitigs between each pair of files
        (diagonal = self-overlap, off-diagonal = cross-file overlap)

    cmp_E_sharing_counts.png
        Bar chart: how many unitigs are unique to 1 file, shared in 2, etc.
        (reveals redundancy and file-specific markers)

    Notes
    -----
    These comparisons help identify consistent GWAS signals and file-specific
    artifacts before proceeding to MSA.
    """
    labels   = [d["label"] for d in datasets]
    n_files  = len(datasets)
    colors   = FILE_COLORS[:n_files]

    print("\n  [COMPARE VIS]")

    # - A  Stacked bar: partition per file -
    fig, ax = plt.subplots(figsize=(max(7, n_files * 2.2), 5))
    x       = np.arange(n_files)
    w       = 0.55
    bot     = np.zeros(n_files)
    for cat, c, name in [
        ("rare",      PAL["red"],  "Rare <5%"),
        ("accessory", PAL["teal"], "Accessory 5–95%"),
        ("core",      PAL["blue"], "Core >95%"),
    ]:
        vals = [d["classes"].value_counts().get(cat, 0) for d in datasets]
        ax.bar(x, vals, w, bottom=bot, color=c, label=name,
               edgecolor="white", lw=0.8)
        for xi, (v, b) in enumerate(zip(vals, bot)):
            if v > 5:
                ax.text(xi, b + v / 2, str(v), ha="center", va="center",
                        fontsize=8.5, fontweight="bold", color="white")
        bot += np.array(vals)

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Unitig count")
    ax.set_title("Pan-genome Partition per File", fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.yaxis.grid(True, alpha=0.4); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / "cmp_A_partition.png")

    # - B  Violin: sample load per file -
    fig, ax = plt.subplots(figsize=(max(7, n_files * 2.5), 5))
    load_data = [d["sample_load"].values for d in datasets]
    parts = ax.violinplot(load_data, positions=range(n_files),
                          showmedians=True, showextrema=True)
    for i, (body, c) in enumerate(zip(parts["bodies"], colors)):
        body.set_facecolor(c); body.set_alpha(0.6)
    parts["cmedians"].set_colors(PAL["navy"])
    parts["cbars"].set_colors(PAL["muted"])
    parts["cmins"].set_colors(PAL["muted"])
    parts["cmaxes"].set_colors(PAL["muted"])
    ax.set_xticks(range(n_files)); ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Unitigs per sample")
    ax.set_title("Per-Sample Unitig Load Distributions", fontweight="bold", pad=10)
    ax.yaxis.grid(True, alpha=0.4); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / "cmp_B_load_violin.png")

    # - C  Frequency overlay -
    fig, ax = plt.subplots(figsize=(9, 5))
    bins = np.linspace(0, 1, 41)
    for d, c in zip(datasets, colors):
        ax.hist(d["freq"].values, bins=bins, alpha=0.55, color=c,
                label=d["label"], edgecolor="none", density=True)
    ax.axvline(0.05, color="black", lw=1.2, ls=":", alpha=0.5)
    ax.axvline(0.95, color="black", lw=1.2, ls=":", alpha=0.5)
    ax.set_xlabel("Presence frequency")
    ax.set_ylabel("Density")
    ax.set_title("Unitig Frequency Overlay - All Files",
                 fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=9)
    ax.yaxis.grid(True, alpha=0.4); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / "cmp_C_freq_overlay.png")

    # D  Pairwise shared unitig heatmap -
    sets    = [set(d["df"].index) for d in datasets]
    mat     = np.zeros((n_files, n_files), dtype=int)
    for i in range(n_files):
        for j in range(n_files):
            mat[i, j] = len(sets[i] & sets[j])

    fig, ax = plt.subplots(figsize=(max(5, n_files * 1.8), max(4, n_files * 1.6)))
    cmap = LinearSegmentedColormap.from_list("amr", ["#FFFFFF", PAL["blue"]])
    im   = ax.imshow(mat, cmap=cmap, vmin=0)
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.03, label="Shared unitig count")
    ax.set_xticks(range(n_files)); ax.set_yticks(range(n_files))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yticklabels(labels)
    thresh = mat.max() / 2
    for i in range(n_files):
        for j in range(n_files):
            ax.text(j, i, str(mat[i, j]), ha="center", va="center",
                    fontsize=11, fontweight="bold",
                    color="white" if mat[i, j] > thresh else PAL["navy"])
    ax.set_title("Pairwise Shared Unitig Count", fontweight="bold", pad=10)
    plt.tight_layout()
    savefig(fig, outdir / "cmp_D_pairwise_shared.png")

    # E  Upset-style: how many files each unitig appears in -
    union = set()
    for d in datasets: union.update(d["df"].index)
    counts_in_n = {}
    for seq in union:
        n = sum(seq in set(d["df"].index) for d in datasets)
        counts_in_n[n] = counts_in_n.get(n, 0) + 1

    fig, ax = plt.subplots(figsize=(7, 4.5))
    xs   = sorted(counts_in_n.keys())
    ys   = [counts_in_n[x] for x in xs]
    cols = [FILE_COLORS[min(i, len(FILE_COLORS)-1)] for i in range(len(xs))]
    bars = ax.bar(xs, ys, color=cols, edgecolor="white", lw=1.2, width=0.65)
    for bar, v in zip(bars, ys):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + max(ys)*0.01,
                f"{v:,}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"In {x} file{'s' if x>1 else ''}" for x in xs], fontsize=9)
    ax.set_ylabel("Unique unitig count")
    ax.set_title("Unitig Sharing Across Files", fontweight="bold", pad=10)
    ax.yaxis.grid(True, alpha=0.4); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / "cmp_E_sharing_counts.png")


# ==================================================================================
#  STAGE 3 - MSA OF SHARED UNITIGS (Biopython)
# ==================================================================================

def run_msa(datasets: list, outdir: Path, min_files: int = 2):
    """
    Multiple Sequence Alignment (MSA) of unitigs shared across files.

    Identifies unitig sequences present in >= min_files files and computes
    pairwise global alignments (Needleman-Wunsch) to quantify conservation.
    Outputs FASTA files for external MSA tools and summarizes alignment
    identity across the panel.

    Parameters
    ----------
    datasets : list of dict
        Output from load_rtab() for each file, containing:
        df (binary unitig matrix), freq, sample_load, classes, label, filepath
    outdir : Path
        Output directory for FASTA, TSV, JSON, and PNG files
    min_files : int, optional
        Minimum number of files a unitig must appear in to include in MSA.
        Default: 2 (any unitig in 2+ files)
        Set to 3 for 3-file datasets to find core sequences only

    Returns
    -------
    tuple
        (df_shared, df_pw, summary)
        - df_shared: DataFrame of sequences in >= min_files files
        - df_pw: DataFrame of pairwise alignment results
        - summary: Dict with total_unique_unitigs, shared counts, length/GC stats

    Output Files
    -------
    union_table.tsv
        All unique unitigs with columns: sequence, length, gc_content,
        n_files, and per-file presence (0/1) for each dataset

    shared_unitigs.fasta
        FASTA-format sequences present in >= min_files files
        Headers: unitig_N files=FILE1,FILE2 len=X gc=Y

    all_unitigs.fasta
        All unique sequences (even if in only 1 file)

    core_shared_unitigs.fasta
        Only sequences in ALL files (if any exist)

    pairwise_identity.tsv
        Pairwise alignments: columns seq_a, seq_b, identity (0-1 scale)
        Capped at 100 comparisons for computational efficiency

    msa_summary.json
        JSON summary: counts, lengths, GC stats (for programmatic access)

    msa_M1_length_dist.png
        Histogram: length distribution of unique vs shared sequences
        (reveals if shared markers are size-biased)

    msa_M2_pw_identity.png
        Heatmap: pairwise alignment identity between top shared sequences
        (white = different, blue = identical; diagonal = 1.0)

    Notes
    -----
    Alignment Method:
      Global mode (Needleman-Wunsch): aligns full sequences end-to-end
      Scoring: match=+2, mismatch=-1, gap_open=-2, gap_extend=-0.5
      Identity normalized to [0, 1] where:
        1.0 = identical sequences
        0.8+ = high conservation (functional orthologs)
        <0.5 = likely unrelated sequences

    For external multiple alignment, use FASTA files with:
      mafft --auto shared_unitigs.fasta > aligned.fasta
      muscle -in shared_unitigs.fasta -out aligned.fasta

    References
    ----------
    Biopython PairwiseAligner: https://biopython.org/Alphabet_PairwiseAligner.html
    MAFFT: https://mafft.cbrc.jp/
    MUSCLE: https://www.drive5.com/muscle/
    """
    labels = [d["label"] for d in datasets]
    sets   = {d["label"]: set(d["df"].index) for d in datasets}

    # Step 3.1: Build comprehensive presence/absence table
    # Rows = unitig sequences, Columns = metadata + presence flags for each file
    all_seqs = set()
    for s in sets.values(): all_seqs.update(s)

    rows = []
    for seq in all_seqs:
        # For each sequence, record which files contain it
        presence = {lbl: int(seq in sets[lbl]) for lbl in labels}
        presence["sequence"]   = seq
        presence["length"]     = len(seq)  # Unitig length in base pairs
        presence["gc_content"] = round(gc(seq), 4)  # GC% (important for DNA stability)
        presence["n_files"]    = sum(presence[l] for l in labels)  # Count files where present
        rows.append(presence)

    df_union = pd.DataFrame(rows).sort_values("n_files", ascending=False)
    df_union.to_csv(outdir / "union_table.tsv", sep="\t", index=False)
    print(f"\n  [MSA] Union table: {len(df_union):,} unique sequences")

    # Step 3.2: Filter for shared sequences
    # Only keep sequences that appear in >= min_files files
    df_shared = df_union[df_union["n_files"] >= min_files].copy()
    df_all    = df_union.copy()
    print(f"  [MSA] Shared (>={min_files} files): {len(df_shared):,} sequences")

    if len(df_shared) == 0:
        print("  [WARN]  No shared sequences found - lowering min_files to 1")
        df_shared = df_union.copy()

    # Check for sequences in ALL files (these are 'core' to all strains)
    df_core_shared = df_union[df_union["n_files"] == len(labels)]
    print(f"  [MSA] In ALL {len(labels)} files: {len(df_core_shared):,} sequences")

    # Step 3.3: Write FASTA files for downstream analysis
    # FASTA format is standard for sequence databases and alignment tools
    # Header includes: unitig ID, files where present, length, GC content
    def write_fasta(df_sub, out_path, tag="shared"):
        records = []
        for _, row in df_sub.iterrows():
            seq    = row["sequence"]
            files_present = [l for l in labels if row.get(l, 0) == 1]
            rec_id = f"unitig_{len(records)+1}"
            rec_desc = (f"files={','.join(files_present)} "
                        f"len={row['length']} gc={row['gc_content']:.3f}")
            records.append(SeqRecord(Seq(seq), id=rec_id, description=rec_desc))
        with open(out_path, "w") as f:
            SeqIO.write(records, f, "fasta")
        print(f"  [FASTA] {out_path.name}  ({len(records)} sequences)")
        return records

    shared_records    = write_fasta(df_shared, outdir / "shared_unitigs.fasta",  "shared")
    all_records       = write_fasta(df_all,    outdir / "all_unitigs.fasta",     "all")
    if len(df_core_shared) > 0:
        core_records  = write_fasta(df_core_shared,
                                    outdir / "core_shared_unitigs.fasta", "core_shared")

    # Step 3.4: Pairwise sequence alignment using global Needleman-Wunsch
    # Global alignment: aligns entire sequences (good for homologous regions)
    # Scoring scheme favours matches (+2) over mismatches (-1) and penalises gaps
    # Identity = proportion of aligned positions that match
    print(f"\n  [MSA] Running pairwise alignments on "
          f"{min(len(shared_records), 100)} shared sequences ...")

    aligner                 = PairwiseAligner()
    aligner.mode            = "global"  # Global = Needleman-Wunsch (full sequence)
    aligner.match_score     = 2  # Reward for matching positions
    aligner.mismatch_score  = -1  # Penalty for mismatches
    aligner.open_gap_score  = -2  # Penalty for opening a gap
    aligner.extend_gap_score = -0.5  # Smaller penalty for extending gap

    pw_results = []
    seqs_to_align = shared_records[:100]  # cap at 100 for computational speed

    # Compare all pairs of sequences
    for i, rec_i in enumerate(seqs_to_align):
        for j, rec_j in enumerate(seqs_to_align):
            if j <= i: continue  # Only upper triangle (avoid duplicates)
            s1, s2 = str(rec_i.seq), str(rec_j.seq)
            try:
                # Calculate alignment score
                score     = aligner.score(s1, s2)
                # Normalize by maximum possible score (both sequences perfectly match)
                max_score = 2 * max(len(s1), len(s2))
                # Identity = normalized score (0=no similarity, 1=identical)
                identity  = score / max_score if max_score > 0 else 0
            except Exception:
                identity = 0.0  # If alignment fails, mark as non-identical
            pw_results.append({
                "seq_a":    rec_i.id,
                "seq_b":    rec_j.id,
                "identity": round(max(0, min(identity, 1)), 4),
            })

    df_pw = pd.DataFrame(pw_results)
    df_pw.to_csv(outdir / "pairwise_identity.tsv", sep="\t", index=False)
    print(f"  [MSA] Pairwise identity table -> pairwise_identity.tsv")

    # Step 3.5: Generate summary statistics and formatted report
    # Create comprehensive summary for interpretation
    summary = {
        "total_unique_unitigs":        int(len(df_union)),
        "shared_in_ge2_files":         int(len(df_shared)),
        "shared_in_all_files":         int(len(df_core_shared)),
        "per_file_counts": {
            l: int(len(sets[l])) for l in labels
        },
        "length_stats": {
            "mean_shared_len":  round(df_shared["length"].mean(), 2),
            "min_shared_len":   int(df_shared["length"].min()),
            "max_shared_len":   int(df_shared["length"].max()),
        },
        "gc_stats": {
            "mean_shared_gc":   round(df_shared["gc_content"].mean(), 4),
            "mean_unique_gc":   round(df_union[df_union["n_files"]==1]["gc_content"].mean(), 4)
                                if (df_union["n_files"]==1).any() else None,
        },
    }
    with open(outdir / "msa_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    
    # Print formatted summary table to console
    print(f"\n  [MSA] Summary Report:")
    print(f"  {'-'*70}")
    print(f"  {'METRIC':<40} {'VALUE':>28}")
    print(f"  {'-'*70}")
    print(f"  {'Total unique unitigs':<40} {len(df_union):>28,}")
    print(f"  {'Shared in >= 2 files':<40} {len(df_shared):>28,}")
    print(f"  {'Shared in ALL files':<40} {len(df_core_shared):>28,}")
    print(f"  {'-'*70}")
    for lbl in labels:
        print(f"  {f'Unitigs in {lbl}':<40} {len(sets[lbl]):>28,}")
    print(f"  {'-'*70}")
    print(f"  {'Shared sequence length (mean)':<40} {df_shared['length'].mean():>28.1f} bp")
    print(f"  {'Shared sequence length (min)':<40} {int(df_shared['length'].min()):>28,} bp")
    print(f"  {'Shared sequence length (max)':<40} {int(df_shared['length'].max()):>28,} bp")
    print(f"  {'-'*70}")
    print(f"  {'GC content (shared, mean)':<40} {df_shared['gc_content'].mean()*100:>27.2f} %")
    print(f"  {'-'*70}")
    print(f"  Outputs saved to: {outdir}/")
    print(f"  {'-'*70}")

    # Step 3.6: MSA Visualisations (limited to 2 most informative plots)
    print(f"\n  [MSA VIS]")

    # PLOT 1: Length distribution of shared vs unique unitigs
    # Shows whether shared unitigs have different size distributions
    # (could indicate selection bias or functional constraint)
    fig, ax = plt.subplots(figsize=(9, 5))
    unique_seqs = df_union[df_union["n_files"] == 1]
    ax.hist(unique_seqs["length"].values, bins=30, alpha=0.6,
            color=PAL["muted"], label=f"Unique ({len(unique_seqs)})",
            edgecolor="white", density=True)
    ax.hist(df_shared["length"].values, bins=30, alpha=0.75,
            color=PAL["teal"], label=f"Shared >={min_files} files ({len(df_shared)})",
            edgecolor="white", density=True)
    if len(df_core_shared) > 0:
        ax.hist(df_core_shared["length"].values, bins=20, alpha=0.85,
                color=PAL["red"],
                label=f"In ALL files ({len(df_core_shared)})",
                edgecolor="white", density=True)
    ax.set_xlabel("Unitig length (bp)")
    ax.set_ylabel("Density")
    ax.set_title("Unitig Length Distribution: Shared vs Unique", fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=9)
    ax.yaxis.grid(True, alpha=0.4); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    savefig(fig, outdir / "msa_M1_length_dist.png")

    # PLOT 2: Pairwise alignment identity heatmap
    # Each cell shows how similar two sequences are (0=different, 1=identical)
    # Helps identify if shared unitigs are truly homologous
    if len(df_pw) > 0 and len(seqs_to_align) >= 2:
        n_show = min(40, len(seqs_to_align))
        ids    = [r.id for r in seqs_to_align[:n_show]]
        id_mat = np.zeros((n_show, n_show))
        np.fill_diagonal(id_mat, 1.0)

        id_lookup = {}
        for _, row in df_pw.iterrows():
            id_lookup[(row["seq_a"], row["seq_b"])] = row["identity"]
            id_lookup[(row["seq_b"], row["seq_a"])] = row["identity"]

        for i in range(n_show):
            for j in range(n_show):
                if i != j:
                    id_mat[i, j] = id_lookup.get((ids[i], ids[j]), 0)

        fig, ax = plt.subplots(figsize=(11, 9))
        cmap = LinearSegmentedColormap.from_list("id", ["#FFFFFF", PAL["blue"]])
        im   = ax.imshow(id_mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")
        plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="Alignment identity")
        ax.set_xticks(range(n_show)); ax.set_yticks(range(n_show))
        ax.set_xticklabels(ids, rotation=90, fontsize=6)
        ax.set_yticklabels(ids, fontsize=6)
        ax.set_title(f"Pairwise Alignment Identity - Top {n_show} Shared Unitigs",
                     fontweight="bold", pad=10)
        plt.tight_layout()
        savefig(fig, outdir / "msa_M2_pw_identity.png")

    print(f"\n  [MSA] Visualisations: 2 informative plots generated")
    print(f"        - Length distribution of unitigs")
    print(f"        - Pairwise sequence alignment identity matrix")

    return df_shared, df_pw, summary


# -
#  MAIN
# -

def main():
    """
    Main entry point: parse arguments, orchestrate 3-stage pipeline.
    
    Stages:
      1. Per-file: Load, validate, and visualize each .Rtab file
      2. Cross-file: Compare unitig distributions across files
      3. MSA: Identify shared sequences and compute pairwise alignments
    """
    parser = argparse.ArgumentParser(
        description="Multi-Rtab GWAS Unitig Comparison Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("files", nargs="+",
                        help="Path(s) to .Rtab files (minimum 2 required)")
    parser.add_argument("--labels", nargs="*", default=None,
                        help="Short labels for each file (same order as files). "
                             "If not provided, derived from filenames.")
    parser.add_argument("--outdir", default="comparison_results",
                        help="Output directory (default: comparison_results)")
    parser.add_argument("--min-files", dest="min_files", type=int, default=2,
                        help="Minimum files a unitig must appear in for MSA analysis "
                             "(default: 2, set to 3 for core sequences across 3 files)")
    args = parser.parse_args()

    # Validation
    if len(args.files) < 2:
        parser.error("Supply at least 2 .Rtab files (e.g., file1.Rtab file2.Rtab file3.Rtab)")
    
    # Check file existence before processing
    missing = [f for f in args.files if not Path(f).exists()]
    if missing:
        parser.error(f"File(s) not found: {', '.join(missing)}")

    # Generate labels from filenames if not provided
    labels = args.labels
    if labels is None:
        labels = [Path(f).stem.replace("_gwas_filtered_unitigs", "")
                           .replace("_filtered_unitigs", "")
                  for f in args.files]
    
    if len(labels) != len(args.files):
        parser.error(f"--labels count ({len(labels)}) must match files count ({len(args.files)})")

    # Create output directory
    outdir = Path(args.outdir)
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        parser.error(f"Cannot create output directory {outdir}: {e}")

    # Banner
    print("\n" + "="*65)
    print("  MULTI-RTAB GWAS UNITIG COMPARISON PIPELINE")
    print("="*65)
    print(f"  Input files : {len(args.files)}")
    print(f"  Labels      : {labels}")
    print(f"  Output dir  : {outdir}/")
    print(f"  Min overlap : {args.min_files} file(s)")
    print("="*65)

    try:
        # Stage 1: Per-file analysis
        print("\n[STAGE 1] Load & per-file visualisation")
        datasets = []
        for fp, lbl in zip(args.files, labels):
            data = load_rtab(fp, lbl)
            visualise_single(data, outdir)
            datasets.append(data)

        # Stage 2: Cross-file comparison
        print("\n[STAGE 2] Cross-file comparison")
        compare_visualise(datasets, outdir)

        # Stage 3: MSA of shared unitigs
        print(f"\n[STAGE 3] MSA (shared in >={args.min_files} files)")
        df_shared, df_pw, summary = run_msa(datasets, outdir, args.min_files)

        # Final report
        print("\n" + "="*65)
        print("  PIPELINE COMPLETE - SUCCESS")
        print("="*65)
        print(f"\n  Key Results:")
        print(f"    Total unique unitigs        : {summary['total_unique_unitigs']:,}")
        print(f"    Shared (>={args.min_files} files)        : {summary['shared_in_ge2_files']:,}")
        print(f"    In ALL {len(labels)} files                : {summary['shared_in_all_files']:,}")
        print(f"\n  Output files:")
        print(f"    {outdir}/*.png                    (18 per-file + 5 comparison plots)")
        print(f"    {outdir}/union_table.tsv         (all sequences + metadata)")
        print(f"    {outdir}/shared_unitigs.fasta    (sequences in >{args.min_files} files)")
        print(f"    {outdir}/pairwise_identity.tsv   (alignment scores)")
        print(f"\n  Next: Run external MSA on shared sequences:")
        print(f"    mafft --auto {outdir}/shared_unitigs.fasta > aligned.fasta")
        print(f"    muscle -in {outdir}/shared_unitigs.fasta -out aligned.fasta")
        print("="*65 + "\n")
        
    except FileNotFoundError as e:
        print(f"\n[ERROR] File not found: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"\n[ERROR] Invalid input: {e}", file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print(f"\n[ERROR] Processing failed: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] Unexpected error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
