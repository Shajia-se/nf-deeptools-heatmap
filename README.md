# nf-deeptools-heatmap

Nextflow DSL2 module for deepTools matrix/heatmap/profile plotting from ChIP-seq bigWig signals.

## What It Does

For each region set (BED), this module runs:

- `computeMatrix` (`reference-point` or `scale-regions`)
- `plotHeatmap` (PNG + PDF)
- `plotProfile` (PNG + PDF)

## Inputs

Required:

- `--bigwig_pattern` (glob for `.bw` files)
- one regions mode:
  - `--regions_bed` (single BED)
  - `--regions_pattern` (glob of BED files)
  - `--regions_sheet` (CSV with `set_name,bed`)

Optional:

- `--samples_label` comma-separated labels in bigWig order
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
