# nf-deeptools-heatmap

Nextflow DSL2 module for deepTools matrix/heatmap/profile plotting from ChIP-seq bigWig signals.

## What It Does

For each region set (BED), this module runs:

- `computeMatrix` (`reference-point` or `scale-regions`)
- `plotHeatmap` (PNG + PDF)
- `plotProfile` (PNG + PDF)

## Inputs

Required:

- `--bigwig_pattern` (glob for `.bw` files) OR `--group_pairs_sheet` (group means)
- one regions mode:
  - `--regions_bed` (single BED)
  - `--regions_pattern` (glob of BED files)
  - `--regions_sheet` (CSV with `set_name,bed`)
  - or `--samples_master` auto regions mode (`<condition>_idr.sorted.chr.bed` from `idr_output`)

Optional:

- `--samples_master` (auto-select bigWigs by sample_id; applies when not using `group_pairs_sheet`)
- `--bigwig_input_dir` (directory used with `samples_master` mode)
- `--deeptools_include_controls` (default `false`)
- `--idr_output` + `--deeptools_region_suffix` for auto regions mode
- `--samples_label` comma-separated labels in bigWig order
- `--group_pairs_sheet` CSV (`group_name,bw1,bw2`) to build mean bigWig per group (via `bigwigCompare --operation mean`)
- matrix settings (`--matrix_mode`, `--reference_point`, `--before_region_start_length`, `--after_region_start_length`, `--region_body_length`, `--bin_size`)
- plot style (`--color_map`, `--dpi`, `--heatmap_height`, `--heatmap_width`, `--profile_height`, `--profile_width`)

## Output

Output root: `${project_folder}/${deeptools_output}`

Per region set (`<set_name>/`):

- `<set_name>.matrix.gz`
- `<set_name>.matrix.tab`
- `<set_name>.sorted_regions.bed`
- `<set_name>.heatmap.png`
- `<set_name>.heatmap.pdf`
- `<set_name>.profile.png`
- `<set_name>.profile.pdf`

If `--group_pairs_sheet` is used:

- `${deeptools_output}/group_mean_bw/<group_name>.mean.bw`

## Run

Default (IDR BED pattern):

```bash
nextflow run main.nf -profile hpc
```

Single BED:

```bash
nextflow run main.nf -profile hpc \
  --regions_bed /ictstr01/groups/idc/projects/uhlenhaut/jiang/pipelines/nf-idr/idr_output/WT_idr.sorted.chr.bed
```

Regions sheet:

```bash
nextflow run main.nf -profile hpc \
  --regions_sheet regions_sheet.example.csv
```

`samples_master`-filtered bigWigs:

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --bigwig_input_dir /path/to/nf-bamcoverage/bamcoverage_output/bigwig \
  --regions_sheet regions_sheet.example.csv
```

Full auto from `samples_master` (bigWigs + regions):

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --bigwig_input_dir /path/to/nf-bamcoverage/bamcoverage_output/bigwig \
  --idr_output /path/to/nf-idr/idr_output
```

Group-mean mode (WT/TG means first, then heatmap/profile):

```bash
nextflow run main.nf -profile hpc \
  --group_pairs_sheet group_pairs_sheet.example.csv \
  --regions_sheet regions_sheet.example.csv \
  --samples_label WT TG
```

Scale-regions mode:

```bash
nextflow run main.nf -profile hpc \
  --matrix_mode scale-regions \
  --region_body_length 5000
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume
```
