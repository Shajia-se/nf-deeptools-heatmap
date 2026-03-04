# nf-deeptools-heatmap

Nextflow DSL2 module for the deepTools mean-track heatmap workflow used in this ChIP-seq pipeline.

This version follows the analysis logic below:

1. collect all enabled non-control sample peaks from `nf-macs3/strict_q0.01`
2. merge them into one union peak set: `all_peaks_merged.bed`
3. count reads-in-peaks for each clean BAM
4. calculate per-sample scaling factors using the minimum reads-in-peaks sample as denominator
5. generate scaled bigWig files with `bamCoverage`
6. compute per-condition mean bigWig tracks from the two replicates
7. run `computeMatrix` on DiffBind gain/loss peak sets
8. produce heatmap and profile plots

## Default Inputs

The module is designed to run automatically from `samples_master.csv` plus upstream module outputs.

Required upstream outputs:
- `nf-chipfilter/chipfilter_output`
- `nf-macs3/macs3_output/strict_q0.01`
- `nf-diffbind/diffbind_output`

Required metadata:
- `--samples_master`

Expected assumptions:
- enabled non-control samples are ChIP samples
- each condition has exactly 2 replicates
- DiffBind contrast file names follow:
  - `condition_unique_up.<treatment>.vs.<reference>.bed`
  - `condition_unique_down.<treatment>.vs.<reference>.bed`

## Main Parameters

- `--samples_master`
- `--chipfilter_output`
- `--macs3_output`
- `--deeptools_macs3_profile`
- `--diffbind_output`
- `--diffbind_contrast`
- `--reference_condition`
- `--treatment_condition`
- `--blacklist` optional

Plot/scaling parameters:
- `--bin_size`
- `--extend_reads`
- `--smooth_length`
- `--center_reads`
- `--ignore_duplicates`
- `--before_region_start_length`
- `--after_region_start_length`
- `--color_map`
- `--plot_title`

## Outputs

Top-level outputs in `deeptools_heatmap_output/`:
- `all_peaks_merged.bed`
- `reads_in_peaks.tsv`
- `scaling_factors.tsv`
- `matrix_meanTracks.gz`
- `regions_sorted.bed`
- `Heatmap_<REF>_vs_<TREAT>_meanTracks.png`
- `Heatmap_<REF>_vs_<TREAT>_meanTracks.pdf`
- `Profile_<REF>_vs_<TREAT>_meanTracks.png`
- `Profile_<REF>_vs_<TREAT>_meanTracks.pdf`

Subdirectories:
- `deeptools_heatmap_output/scaled_bigwig/`
- `deeptools_heatmap_output/mean_tracks/`

## Example

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output \
  --macs3_output /path/to/nf-macs3/macs3_output \
  --diffbind_output /path/to/nf-diffbind/diffbind_output \
  --reference_condition WT \
  --treatment_condition TG
```
