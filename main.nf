#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def outdir = params.deeptools_output ?: 'deeptools_heatmap_output'

process deeptools_matrix_plot {
  tag "${set_name}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${outdir}/${set_name}", mode: 'copy', overwrite: true

  input:
    tuple val(set_name), path(regions_bed), path(bigwigs)

  output:
    path("${set_name}.matrix.gz")
    path("${set_name}.matrix.tab")
    path("${set_name}.heatmap.png")
    path("${set_name}.heatmap.pdf")
    path("${set_name}.profile.png")
    path("${set_name}.profile.pdf")
    path("${set_name}.sorted_regions.bed")

  script:
    def bw_args = bigwigs.collect { "${it}" }.join(' ')
    def label_args = ''
    if (params.samples_label) {
      def labels = params.samples_label.toString().split(',').collect { it.trim() }.findAll { it }
      if (labels) {
        label_args = '--samplesLabel ' + labels.join(' ')
      }
    }

    def mode = (params.matrix_mode ?: 'reference-point').toString()
    def matrix_cmd = ''

    if (mode == 'scale-regions') {
      matrix_cmd = """
      computeMatrix scale-regions \\
        -S ${bw_args} \\
        -R ${regions_bed} \\
        -o ${set_name}.matrix.gz \\
        --outFileNameMatrix ${set_name}.matrix.tab \\
        --outFileSortedRegions ${set_name}.sorted_regions.bed \\
        -a ${params.after_region_start_length} \\
        -b ${params.before_region_start_length} \\
        -m ${params.region_body_length} \\
        --binSize ${params.bin_size} \\
        --skipZeros \\
        --numberOfProcessors ${task.cpus}
      """
    } else {
      matrix_cmd = """
      computeMatrix reference-point \\
        --referencePoint ${params.reference_point} \\
        -S ${bw_args} \\
        -R ${regions_bed} \\
        -o ${set_name}.matrix.gz \\
        --outFileNameMatrix ${set_name}.matrix.tab \\
        --outFileSortedRegions ${set_name}.sorted_regions.bed \\
        -a ${params.after_region_start_length} \\
        -b ${params.before_region_start_length} \\
        --binSize ${params.bin_size} \\
        --skipZeros \\
        --numberOfProcessors ${task.cpus}
      """
    }

  """
  set -euo pipefail

  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp
  export MPLCONFIGDIR=\$PWD/tmp/mpl
  mkdir -p \"\$MPLCONFIGDIR\"

  ${matrix_cmd}

  plotHeatmap \\
    -m ${set_name}.matrix.gz \\
    -out ${set_name}.heatmap.png \\
    --dpi ${params.dpi} \\
    --colorMap ${params.color_map} \\
    --heatmapHeight ${params.heatmap_height} \\
    --heatmapWidth ${params.heatmap_width} \\
    ${label_args}

  plotHeatmap \\
    -m ${set_name}.matrix.gz \\
    -out ${set_name}.heatmap.pdf \\
    --colorMap ${params.color_map} \\
    --heatmapHeight ${params.heatmap_height} \\
    --heatmapWidth ${params.heatmap_width} \\
    ${label_args}

  plotProfile \\
    -m ${set_name}.matrix.gz \\
    -out ${set_name}.profile.png \\
    --dpi ${params.dpi} \\
    --plotHeight ${params.profile_height} \\
    --plotWidth ${params.profile_width} \\
    ${label_args}

  plotProfile \\
    -m ${set_name}.matrix.gz \\
    -out ${set_name}.profile.pdf \\
    --plotHeight ${params.profile_height} \\
    --plotWidth ${params.profile_width} \\
    ${label_args}
  """
}

workflow {
  if (!params.bigwig_pattern) {
    exit 1, 'ERROR: Missing --bigwig_pattern'
  }

  def ch_bigwigs = Channel
    .fromPath(params.bigwig_pattern, checkIfExists: true)
    .ifEmpty { exit 1, "ERROR: No bigWig files found with pattern: ${params.bigwig_pattern}" }
    .collect()

  def ch_regions

  if (params.regions_bed) {
    ch_regions = Channel
      .fromPath(params.regions_bed, checkIfExists: true)
      .ifEmpty { exit 1, "ERROR: regions_bed not found: ${params.regions_bed}" }
      .map { bed -> tuple(bed.baseName, bed) }
  } else if (params.regions_pattern) {
    ch_regions = Channel
      .fromPath(params.regions_pattern, checkIfExists: true)
      .ifEmpty { exit 1, "ERROR: No region BED files found with pattern: ${params.regions_pattern}" }
      .map { bed -> tuple(bed.baseName, bed) }
  } else if (params.regions_sheet) {
    ch_regions = Channel
      .fromPath(params.regions_sheet, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
        assert row.set_name && row.bed : 'regions_sheet must contain: set_name,bed'
        def bed = file(row.bed.toString())
        assert bed.exists() : "BED not found for ${row.set_name}: ${bed}"
        tuple(row.set_name.toString().trim(), bed)
      }
  } else {
    exit 1, 'ERROR: Provide one of --regions_bed, --regions_pattern, or --regions_sheet'
  }

  ch_regions
    .combine(ch_bigwigs)
    .set { ch_tasks }

  deeptools_matrix_plot(ch_tasks)
}
