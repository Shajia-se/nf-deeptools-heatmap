#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def outdir = params.deeptools_output ?: 'deeptools_heatmap_output'

process deeptools_group_mean_bw {
  tag "${group_name}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${outdir}/group_mean_bw", mode: 'copy', overwrite: true

  input:
    tuple val(group_name), path(bw1), path(bw2)

  output:
    tuple val(group_name), path("${group_name}.mean.bw")

  script:
  """
  set -euo pipefail

  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp
  export MPLCONFIGDIR=\$PWD/tmp/mpl
  mkdir -p "\$MPLCONFIGDIR"

  bigwigCompare \
    -b1 ${bw1} \
    -b2 ${bw2} \
    --operation mean \
    --skipNonCoveredRegions \
    --binSize ${params.bin_size} \
    --numberOfProcessors ${task.cpus} \
    -o ${group_name}.mean.bw
  """
}

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
  def ch_bigwigs
  if (params.group_pairs_sheet) {
    def ch_pairs = Channel
      .fromPath(params.group_pairs_sheet, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
        assert row.group_name && row.bw1 && row.bw2 : 'group_pairs_sheet must contain: group_name,bw1,bw2'
        def b1 = file(row.bw1.toString())
        def b2 = file(row.bw2.toString())
        assert b1.exists() : "bw1 not found for ${row.group_name}: ${b1}"
        assert b2.exists() : "bw2 not found for ${row.group_name}: ${b2}"
        tuple(row.group_name.toString().trim(), b1, b2)
      }

    ch_bigwigs = deeptools_group_mean_bw(ch_pairs)
      .map { group_name, bw -> bw }
      .collect()
      .map { bw_list -> tuple(bw_list) }
  } else {
    if (params.samples_master) {
      def master = file(params.samples_master)
      assert master.exists() : "samples_master not found: ${params.samples_master}"

      def header = null
      def records = []
      master.eachLine { line, n ->
        if (!line?.trim()) return
        def cols = line.split(',', -1)*.trim()
        if (n == 1) {
          header = cols
        } else {
          def rec = [:]
          header.eachWithIndex { h, i -> rec[h] = i < cols.size() ? cols[i] : '' }
          records << rec
        }
      }

      assert header : "samples_master header not found: ${params.samples_master}"
      assert header.contains('sample_id') : "samples_master missing required column: sample_id"

      def isEnabled = { rec ->
        def v = rec.enabled?.toString()?.trim()?.toLowerCase()
        (v == null || v == '' || v == 'true')
      }
      def isControl = { rec ->
        rec.is_control?.toString()?.trim()?.toLowerCase() == 'true'
      }

      def includeControls = (params.deeptools_include_controls == null) ? false : params.deeptools_include_controls
      def bwDir = file(params.bigwig_input_dir)
      assert bwDir.exists() : "bigwig_input_dir not found: ${params.bigwig_input_dir}"

      def bwRows = records
        .findAll { rec -> isEnabled(rec) }
        .findAll { rec -> includeControls ? true : !isControl(rec) }
        .collect { rec ->
          def sid = rec.sample_id?.toString()?.trim()
          if (!sid) return null
          def hits = bwDir.listFiles()?.findAll { f ->
            f.isFile() && f.name.endsWith('.bw') && (f.name == "${sid}.bw" || f.name.startsWith("${sid}_"))
          } ?: []
          if (hits.isEmpty()) throw new IllegalArgumentException("No bigWig found for sample_id '${sid}' under: ${params.bigwig_input_dir}")
          if (hits.size() > 1) throw new IllegalArgumentException("Multiple bigWigs matched sample_id '${sid}': ${hits*.name.join(', ')}")
          file(hits[0].toString())
        }
        .findAll { it != null }

      ch_bigwigs = Channel
        .fromList(bwRows)
        .ifEmpty { exit 1, "ERROR: No bigWigs generated from samples_master: ${params.samples_master}" }
        .collect()
        .map { bw_list -> tuple(bw_list) }
    } else {
      if (!params.bigwig_pattern) {
        exit 1, 'ERROR: Missing --bigwig_pattern'
      }
      ch_bigwigs = Channel
        .fromPath(params.bigwig_pattern, checkIfExists: true)
        .ifEmpty { exit 1, "ERROR: No bigWig files found with pattern: ${params.bigwig_pattern}" }
        .collect()
        .map { bw_list -> tuple(bw_list) }
    }
  }

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
  } else if (params.samples_master) {
    def master = file(params.samples_master)
    assert master.exists() : "samples_master not found: ${params.samples_master}"

    def header = null
    def records = []
    master.eachLine { line, n ->
      if (!line?.trim()) return
      def cols = line.split(',', -1)*.trim()
      if (n == 1) {
        header = cols
      } else {
        def rec = [:]
        header.eachWithIndex { h, i -> rec[h] = i < cols.size() ? cols[i] : '' }
        records << rec
      }
    }
    assert header : "samples_master header not found: ${params.samples_master}"
    assert header.contains('condition') : "samples_master missing required column: condition"

    def isEnabled = { rec ->
      def v = rec.enabled?.toString()?.trim()?.toLowerCase()
      (v == null || v == '' || v == 'true')
    }
    def isControl = { rec ->
      rec.is_control?.toString()?.trim()?.toLowerCase() == 'true'
    }

    def conds = records
      .findAll { rec -> isEnabled(rec) && !isControl(rec) }
      .collect { rec -> rec.condition?.toString()?.trim() }
      .findAll { it }
      .unique()

    def regionRows = conds.collect { cond ->
      def bed = file("${params.idr_output}/${cond}${params.deeptools_region_suffix}")
      assert bed.exists() : "Auto regions BED not found for condition '${cond}': ${bed}"
      tuple(cond, bed)
    }

    ch_regions = Channel
      .fromList(regionRows)
      .ifEmpty { exit 1, "ERROR: No regions generated from samples_master. Provide --regions_sheet or check samples_master/idr_output." }
  } else {
    exit 1, 'ERROR: Provide one of --regions_bed, --regions_pattern, --regions_sheet, or --samples_master (auto regions).'
  }

  ch_tasks = ch_regions
    .combine(ch_bigwigs)
    .map { items ->
      def set_name = items[0]
      def regions_bed = items[1]
      def bw_list = items.size() > 2 ? items.subList(2, items.size()) : []
      tuple(set_name, regions_bed, bw_list)
    }

  deeptools_matrix_plot(ch_tasks)
}
