# Multi-Rtab GWAS Unitig Comparison Pipeline

A Python toolkit for comparative analysis of GWAS unitig panels across multiple antibiotic resistance phyla. Performs dimensionality reduction, pan-genome characterization, and pairwise sequence alignment of shared markers.

**Status**: Production-ready | **Python**: 3.8+ | **License**: MIT

---

## Overview

This pipeline enables systematic comparison of unitig presence/absence matrices (.Rtab format) across multiple GWAS datasets. Common use case: analyzing resistance markers for different antibiotics (e.g., Azithromycin, Cefixime, Ciprofloxacin) in bacterial populations.

### Key Features

- **Per-file analysis**: Frequency distributions, pan-genome partitioning (core/accessory/rare), sample load profiles
- **Cross-file comparison**: Shared unitig identification, pairwise overlap matrices, frequency overlays
- **Multiple Sequence Alignment (MSA)**: Biopython Needleman-Wunsch alignment of shared sequences with identity matrices
- **Dimensionality reduction**: PCA with scree plots for sample space visualization
- **Publication-quality plots**: Matplotlib/Seaborn outputs (PNG, 150 DPI)

---

## Installation

### Requirements
- Python ≥ 3.8
- Dependencies listed in `requirements.txt`

### Setup

```bash
git clone https://github.com/yourusername/multi-rtab-compare.git
cd multi-rtab-compare
pip install -r requirements.txt
```

---

## Usage

### Basic Command

```bash
python multi_rtab_compare.py file1.Rtab file2.Rtab file3.Rtab
```

### With Custom Labels & Output Directory

```bash
python multi_rtab_compare.py \
  azm_sr.Rtab cip_sr.Rtab cfx_sr.Rtab \
  --labels AZM CIP CFX \
  --outdir results_2024
```

### Require Stricter Sharing (All 3 Files)

```bash
python multi_rtab_compare.py file1.Rtab file2.Rtab file3.Rtab --min-files 3
```

---

## Input Format: .Rtab Files

Space-separated matrix of binary unitig presence/absence:

```
unitig_sequence      sample_1  sample_2  sample_3  ...
ATCGATCGATCG...      1         0         1
GGCTAGGCTAG...       1         1         0
```

**Requirements**:
- First column = unitig DNA sequence (no spaces)
- Remaining columns = binary (0/1) presence per sample
- Header row with sample names
- Consistent column count across all files (sample count must match)

---

## Output Structure

```
comparison_results/
├── Stage 1 (Per-file analysis)
│   ├── azm_sr_1_freq_hist.png          # Unitig frequency distribution
│   ├── azm_sr_2_partition.png          # Core/accessory/rare breakdown
│   ├── azm_sr_3_sample_load.png        # Unitigs per sample
│   ├── azm_sr_4_prevalence_rank.png    # Ranked prevalence curve
│   ├── azm_sr_5_pca.png                # PCA + scree plot
│   ├── azm_sr_6_heatmap.png            # Clustered heatmap (top 50)
│   └── [repeat for cfx_sr, cip_sr]
│
├── Stage 2 (Cross-file comparison)
│   ├── cmp_A_partition.png             # Stacked bar chart
│   ├── cmp_B_load_violin.png           # Per-sample load distributions
│   ├── cmp_C_freq_overlay.png          # Frequency histogram overlay
│   ├── cmp_D_pairwise_shared.png       # Shared unitig heatmap
│   └── cmp_E_sharing_counts.png        # Unitig presence distribution
│
├── Stage 3 (MSA analysis)
│   ├── union_table.tsv                 # All unique unitigs + metadata
│   ├── shared_unitigs.fasta            # Sequences in >=2 files (FASTA)
│   ├── all_unitigs.fasta               # All unique sequences (FASTA)
│   ├── core_shared_unitigs.fasta       # Core sequences if any (FASTA)
│   ├── pairwise_identity.tsv           # Pairwise alignment scores
│   ├── msa_summary.json                # Summary statistics
│   ├── msa_M1_length_dist.png          # Length distribution comparison
│   └── msa_M2_pw_identity.png          # Identity heatmap
└── comparison_results/
```

---

## Methodology

### Stage 1: Per-File Unitig Analysis

For each .Rtab file:
1. **Loading & Validation**
   - Parse space-separated matrix
   - Coerce to binary (>0 = present)
   - Remove duplicates and all-absent unitigs
   
2. **Classification**
   - Core: prevalence > 95% (essential markers)
   - Accessory: 5-95% (common but variable)
   - Rare: < 5% (population-specific)
   
3. **Visualization**
   - Frequency histogram (log scale)
   - Pan-genome proportional stacked bar
   - Per-sample unitig load distribution
   - Ranked prevalence curve
   - PCA with scree plot (explains sample relationships)
   - Clustered heatmap (top 50 variable unitigs by Jaccard distance)

### Stage 2: Cross-File Comparison

1. **Partition overlay** - Compare core/accessory/rare ratios across files
2. **Sample load violin plot** - Distribution of unitig counts per sample by file
3. **Frequency overlay** - Histogram comparison across all files
4. **Pairwise shared count** - Matrix of overlap between files
5. **Sharing distribution** - How many files each unitig appears in

### Stage 3: Multiple Sequence Alignment (MSA)

For sequences shared in ≥ min_files (default: 2):

1. **Sequence Union**
   - Extract all unique sequences across files
   - Record presence/absence per file
   - Calculate GC content and length for each

2. **Pairwise Alignment** (Biopython PairwiseAligner)
   - Mode: Global (Needleman-Wunsch)
   - Scoring: match=+2, mismatch=-1, gap_open=-2, gap_extend=-0.5
   - Identity = (alignment_score) / (2 × max_seq_length)
   - Range: 0 (completely different) to 1 (identical)

3. **Output**
   - FASTA files for external tools (MAFFT, MUSCLE)
   - Pairwise identity matrix
   - Length distribution comparison
   - Alignment identity heatmap

---

## Example: Analyzing Antibiotic Resistance Markers

```bash
# Compare resistance unitigs across 3 antibiotics
python multi_rtab_compare.py \
  azm_sr_gwas.Rtab \      # Azithromycin resistance markers
  cfx_sr_gwas.Rtab \      # Cefixime resistance markers  
  cip_sr_gwas.Rtab \      # Ciprofloxacin resistance markers
  --labels AZM CFX CIP \
  --outdir resistance_comparison

# Question: Are overlapping markers true orthologs?
# → Check msa_M2_pw_identity.png for conservation patterns
# → High identity (>0.95) = likely same gene
# → Low identity (<0.80) = may be phenotypic convergence
```

---

## Technical Details

### Dimensionality Reduction (PCA)

Principal Component Analysis reduces samples from thousands of binary dimensions down to 2-3 principal components while preserving maximum variance. Interprets as:

- **PC1 & PC2 scatter plot**: Samples positioned by genetic similarity. Close clusters = genetically similar strains; distant clusters = divergent lineages.
- **Scree plot**: Shows variance explained %. Typically, 2-3 components capture 60-80% of variance in bacterial populations.
- **Use case**: Identify population structure without phylogenetic tree inference.

### Alignment Identity Calculation

For sequences S1 and S2:
```
alignment_score = Σ(matches × +2) + Σ(mismatches × -1) + Σ(gaps × -2 or -0.5)
max_score = 2 × length(longer_sequence)
identity = alignment_score / max_score
```

Values normalize to [0, 1]:
- **1.0** = identical sequences
- **0.8+** = highly conserved (likely functional orthologs)
- **0.5-0.8** = moderate divergence (possible gene family)
- **<0.5** = independent/unrelated

---

## Troubleshooting

### Error: "Dropped N duplicate unitigs"
**Cause**: Same sequence present multiple times in one .Rtab file (artifact of k-mer overlap)  
**Fix**: Use `--keep-duplicates` flag if intentional, otherwise normal behavior

### Error: "No shared sequences found"
**Cause**: Files have <2 unitigs in common (expected for diverse resistance profiles)  
**Fix**: Lower `--min-files` to 1, or check file validity with `head comparison_results/union_table.tsv`

### PCA skipped with error
**Cause**: Fewer samples than unitigs, or all unitigs identical across samples  
**Fix**: Check data quality; may indicate low phenotypic diversity

### Graph output is blank/corrupted
**Cause**: Matplotlib backend rendering issue on headless systems  
**Fix**: Already set to "Agg" backend (non-interactive); check disk space for PNG output

---

## Performance

Typical runtimes (3 files, ~500-9000 unitigs, ~4000 samples):

| Stage | Time | Notes |
|-------|------|-------|
| Stage 1 (per-file) | 30-60 sec | Dominated by heatmap clustering |
| Stage 2 (comparison) | 10-20 sec | Fast matrix operations |
| Stage 3 (MSA) | 20-120 sec | Scales with pairwise comparisons |
| **Total** | **2-4 min** | Linear in files, quadratic in shared sequences |

Memory usage: ~500 MB - 2 GB depending on file sizes.

---

## Dependencies

See `requirements.txt` for exact versions. Key packages:

- **pandas** - Data frame manipulation
- **numpy** - Numerical operations
- **matplotlib/seaborn** - Visualization
- **scipy** - Clustering (hierarchical dendrogram)
- **scikit-learn** - PCA implementation
- **biopython** - Sequence alignment (PairwiseAligner)

---

## Citation

If you use this pipeline, please cite:

```bibtex
@software{raut2026multirtab,
  author = {Raut, Meetrayu},
  title = {Multi-Rtab GWAS Unitig Comparison Pipeline},
  year = {2026},
  url = {https://github.com/yourusername/multi-rtab-compare}
}
```

---

## License

MIT License - See [LICENSE](LICENSE) file

---

## Contributing

Bug reports and feature requests welcome. Pull requests appreciated.

```bash
# Local development setup
git clone <repo>
pip install -r requirements.txt
python -m pytest tests/  # When available

# Run on test data
python multi_rtab_compare.py sample_data/*.Rtab --outdir test_output
```

---

## FAQs

**Q: Can I use .Rtab files from non-resistance GWAS?**  
A: Yes - any binary presence/absence matrix works. Workflow is agnostic to phenotype.

**Q: Should I filter rare unitigs before input?**  
A: No - pipeline handles classification automatically. Keep all sequences as-is from GWAS caller.

**Q: What if samples don't match across files?**  
A: Pipeline will fail if column counts differ. Ensure same sample set across all inputs.

**Q: Can I extend this for >3 files?**  
A: Yes - loops support arbitrary file count. Linear performance up to ~10 files; test beyond that.

**Q: How do I use FASTA outputs for external MSA?**  
A: 
```bash
mafft --auto comparison_results/shared_unitigs.fasta > aligned.fasta
# or
muscle -in comparison_results/shared_unitigs.fasta -out aligned.fasta
```

---

## Contact

Meetrayu Raut | [email protected]

---

**Last Updated**: April 2026  
**Tested on**: Python 3.10-3.13, Windows 10/11 + Linux
