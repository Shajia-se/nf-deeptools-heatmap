# nf-deeptools-heatmap

`nf-deeptools-heatmap` creates scaled BigWig tracks, mean tracks, and heatmap/profile plots with deepTools.

This module is optional and usually runs after `nf-chipfilter`, `nf-macs3`, and `nf-diffbind`.

## Required Inputs

```bash
--samples_master /path/to/samples_master.csv
--chipfilter_output /path/to/chipfilter_output
--macs3_output /path/to/macs3_output
```

If using DiffBind gain/loss regions, also provide:

```bash
--diffbind_output /path/to/diffbind_output
```

## Output

```text
${project_folder}/${deeptools_output}/
```

Includes scaled BigWigs, mean tracks, heatmaps, and profile plots.

By default, `reference_condition`, `treatment_condition`, and `diffbind_contrast` are inferred from `samples_master`. Set them explicitly only when you need to override the automatic two-condition comparison.

## Run

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --macs3_output /path/to/macs3_output \
  --diffbind_output /path/to/diffbind_output \
  --project_folder /path/to/output_project
```

Actual execution should be tested where Nextflow is installed.
