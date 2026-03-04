#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def outdir = params.deeptools_output ?: 'deeptools_heatmap_output'

process merge_master_peakset {
  tag "master_peakset"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${outdir}", mode: 'copy', overwrite: true

  container "${params.container_bedtools}"

  input:
    path peak_files

  output:
    path("all_peaks_merged.bed"), emit: master_bed

  script:
  def peakArgs = peak_files.collect { "${it}" }.join(' ')
  """
  set -euo pipefail

  cat ${peakArgs} \
    | awk 'BEGIN{OFS="\\t"} NF >= 3 {print \$1,\$2,\$3}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > all_peaks_merged.bed
  """
}

process count_reads_in_peaks {
  tag "${sample_id}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${outdir}", mode: 'copy', overwrite: true, pattern: "*.reads_in_peaks.tsv"

  container "${params.container_bedtools}"

  input:
    tuple val(sample_id), val(condition), val(replicate), path(bam), path(master_bed)

  output:
    path("${sample_id}.reads_in_peaks.tsv"), emit: reads_tsv

  script:
  """
  set -euo pipefail

  reads_in_peaks=\$(bedtools intersect -u -a ${bam} -b ${master_bed} | samtools view -c -)

  printf "sample\\tcondition\\treplicate\\treads_in_peaks\\n%s\\t%s\\t%s\\t%s\\n" \
    "${sample_id}" "${condition}" "${replicate}" "\${reads_in_peaks}" > ${sample_id}.reads_in_peaks.tsv
  """
}

process compute_scaling_factors {
  tag "scaling_factors"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${outdir}", mode: 'copy', overwrite: true

  container "${params.container_bedtools}"

  input:
    path(read_count_files)

  output:
    path("reads_in_peaks.tsv"), emit: reads_in_peaks_tsv
    path("scaling_factors.tsv"), emit: scaling_factors_tsv

  script:
  def countArgs = read_count_files.collect { "\"${it}\"" }.join(' ')
  """
  set -euo pipefail

  {
    printf "sample\\tcondition\\treplicate\\treads_in_peaks\\n"
    for f in ${countArgs}; do
      awk 'NR > 1 {print}' "\$f"
    done
  } > reads_in_peaks.tsv

  lowest=\$(awk 'NR > 1 {print \$4}' reads_in_peaks.tsv | sort -n | head -n 1)
  if [[ -z "\$lowest" || "\$lowest" == "0" ]]; then
    echo "ERROR: lowest reads_in_peaks is empty or zero" >&2
    exit 1
  fi

  awk -F'\\t' -v OFS='\\t' -v low="\$lowest" '
    BEGIN { print "sample","condition","replicate","reads_in_peaks","scaling_factor" }
    NR > 1 { printf "%s\\t%s\\t%s\\t%s\\t%.8f\\n", \$1,\$2,\$3,\$4,(low/\$4) }
  ' reads_in_peaks.tsv > scaling_factors.tsv
  """
}

process deeptools_scaled_bw {
  tag "${sample_id}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${outdir}/scaled_bigwig", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), val(condition), val(replicate), path(bam), path(scale_table)

  output:
    path("${sample_id}.scaled.bw"), emit: scaled_bw

  script:
  def blacklistArg = params.blacklist ? "--blackListFileName ${params.blacklist}" : ""
  def centerArg = params.center_reads ? "--centerReads" : ""
  def ignoreDupArg = params.ignore_duplicates ? "--ignoreDuplicates" : ""
  """
  set -euo pipefail

  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp
  export MPLCONFIGDIR=\$PWD/tmp/mpl
  mkdir -p "\$MPLCONFIGDIR"

  bam_local="${bam}"
  bam_abs=\$(readlink -f "\$bam_local" || echo "\$bam_local")
  if [[ -f "\${bam_local}.bai" ]]; then
    :
  elif [[ -f "\${bam_local%.bam}.bai" ]]; then
    ln -sf "\${bam_local%.bam}.bai" "\${bam_local}.bai"
  elif [[ -f "\${bam_abs}.bai" ]]; then
    ln -sf "\${bam_abs}.bai" "\${bam_local}.bai"
  elif [[ -f "\${bam_abs%.bam}.bai" ]]; then
    ln -sf "\${bam_abs%.bam}.bai" "\${bam_local}.bai"
  else
    echo "ERROR: BAM index not found for \$bam_local" >&2
    exit 1
  fi

  scale_factor=\$(awk -F'\\t' 'NR > 1 && \$1 == "${sample_id}" {print \$5; exit}' ${scale_table})
  if [[ -z "\$scale_factor" ]]; then
    echo "ERROR: scaling factor not found for ${sample_id}" >&2
    exit 1
  fi

  bamCoverage \
    -b ${bam} \
    -o ${sample_id}.scaled.bw \
    --binSize ${params.bin_size} \
    --scaleFactor "\$scale_factor" \
    --normalizeUsing None \
    --extendReads ${params.extend_reads} \
    --smoothLength ${params.smooth_length} \
    ${centerArg} \
    ${ignoreDupArg} \
    ${blacklistArg} \
    --numberOfProcessors ${task.cpus}
  """
}

process deeptools_group_mean_bw {
  tag "group_mean"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${outdir}/mean_tracks", mode: 'copy', overwrite: true

  input:
    path(scaled_bws)
    val(sample_manifest_tsv)

  output:
    path("*.mean.bw"), emit: mean_bws

  script:
  """
  set -euo pipefail

  cat > sample_manifest.tsv <<'TSV'
${sample_manifest_tsv}
TSV

  while read -r condition; do
    mapfile -t samples < <(awk -F'\\t' -v c="\$condition" 'NR > 1 && \$2 == c {print \$1}' sample_manifest.tsv)
    if [[ \${#samples[@]} -ne 2 ]]; then
      echo "ERROR: expected exactly 2 replicates for condition '\$condition', found \${#samples[@]}" >&2
      exit 1
    fi

    bigwigCompare \
      -b1 "\${samples[0]}.scaled.bw" \
      -b2 "\${samples[1]}.scaled.bw" \
      --operation mean \
      --skipNonCoveredRegions \
      --binSize ${params.bin_size} \
      --numberOfProcessors ${task.cpus} \
      -o "\${condition}.mean.bw"
  done < <(awk -F'\\t' 'NR > 1 {print \$2}' sample_manifest.tsv | sort -u)
  """
}

process deeptools_matrix_plot {
  tag "${reference_condition}_vs_${treatment_condition}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${outdir}", mode: 'copy', overwrite: true

  input:
    path(mean_bws)
    path(up_bed)
    path(down_bed)
    val(reference_condition)
    val(treatment_condition)

  output:
    path("matrix_meanTracks.gz")
    path("regions_sorted.bed")
    path("Heatmap_*.png")
    path("Heatmap_*.pdf")
    path("Profile_*.png")
    path("Profile_*.pdf")

  script:
  def heatmapPng = "Heatmap_${reference_condition}_vs_${treatment_condition}_meanTracks.png"
  def heatmapPdf = "Heatmap_${reference_condition}_vs_${treatment_condition}_meanTracks.pdf"
  def profilePng = "Profile_${reference_condition}_vs_${treatment_condition}_meanTracks.png"
  def profilePdf = "Profile_${reference_condition}_vs_${treatment_condition}_meanTracks.pdf"
  def title = params.plot_title ?: "Mean Tracks ${reference_condition} vs ${treatment_condition} on Gain/Loss Peaks"
  def gainLabel = params.gain_label ?: "${treatment_condition} up"
  def lossLabel = params.loss_label ?: "${reference_condition} up"
  """
  set -euo pipefail

  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp
  export MPLCONFIGDIR=\$PWD/tmp/mpl
  mkdir -p "\$MPLCONFIGDIR"

  computeMatrix reference-point \
    --referencePoint ${params.reference_point} \
    -S ${reference_condition}.mean.bw ${treatment_condition}.mean.bw \
    -R ${up_bed} ${down_bed} \
    --regionsLabel "${gainLabel}" "${lossLabel}" \
    --samplesLabel "${reference_condition} (mean)" "${treatment_condition} (mean)" \
    -b ${params.before_region_start_length} \
    -a ${params.after_region_start_length} \
    --skipZeros \
    --binSize ${params.bin_size} \
    --numberOfProcessors ${task.cpus} \
    -o matrix_meanTracks.gz \
    --outFileSortedRegions regions_sorted.bed

  plotHeatmap \
    -m matrix_meanTracks.gz \
    -out ${heatmapPng} \
    --dpi ${params.dpi} \
    --colorMap ${params.color_map} \
    --heatmapWidth ${params.heatmap_width} \
    --heatmapHeight ${params.heatmap_height} \
    --plotTitle "${title}"

  plotHeatmap \
    -m matrix_meanTracks.gz \
    -out ${heatmapPdf} \
    --colorMap ${params.color_map} \
    --heatmapWidth ${params.heatmap_width} \
    --heatmapHeight ${params.heatmap_height} \
    --plotTitle "${title}"

  plotProfile \
    -m matrix_meanTracks.gz \
    -out ${profilePng} \
    --dpi ${params.dpi} \
    --plotWidth ${params.profile_width} \
    --plotHeight ${params.profile_height} \
    --plotTitle "${title}"

  plotProfile \
    -m matrix_meanTracks.gz \
    -out ${profilePdf} \
    --plotWidth ${params.profile_width} \
    --plotHeight ${params.profile_height} \
    --plotTitle "${title}"
  """
}

workflow {
  assert params.samples_master : "Missing --samples_master"

  def master = file(params.samples_master)
  assert master.exists() : "samples_master not found: ${params.samples_master}"

  def records = parseSamplesMaster(master)
  def chipRecords = records
    .findAll { rec -> isEnabled(rec) && !isControl(rec) }
    .sort { a, b ->
      def c = (a.condition ?: '') <=> (b.condition ?: '')
      if (c != 0) return c
      (a.replicate ?: Integer.MAX_VALUE) <=> (b.replicate ?: Integer.MAX_VALUE)
    }
  assert !chipRecords.isEmpty() : "No enabled non-control samples found in samples_master"

  def conditionMap = chipRecords.groupBy { it.condition }
  conditionMap.each { cond, list ->
    assert list.size() == 2 : "deepTools mean-track mode expects exactly 2 replicates per condition. Condition '${cond}' has ${list.size()}."
  }

  def sortedConditions = conditionMap.keySet().toList().sort()
  def referenceCondition = params.reference_condition ?: (sortedConditions.contains('WT') ? 'WT' : sortedConditions[0])
  def treatmentCondition = params.treatment_condition ?: (sortedConditions.contains('TG') ? 'TG' : sortedConditions[-1])
  assert referenceCondition != treatmentCondition : "reference_condition and treatment_condition must differ"

  def contrastName = params.diffbind_contrast ?: "${treatmentCondition}.vs.${referenceCondition}"
  def diffbindDir = file(params.diffbind_output)
  assert diffbindDir.exists() : "diffbind_output not found: ${params.diffbind_output}"
  def upBed = file("${params.diffbind_output}/${params.diffbind_up_prefix ?: 'condition_unique_up'}.${contrastName}.bed")
  def downBed = file("${params.diffbind_output}/${params.diffbind_down_prefix ?: 'condition_unique_down'}.${contrastName}.bed")
  assert upBed.exists() : "DiffBind up BED not found: ${upBed}"
  assert downBed.exists() : "DiffBind down BED not found: ${downBed}"

  def bamDir = file(params.chipfilter_output)
  assert bamDir.exists() : "chipfilter_output not found: ${params.chipfilter_output}"

  def peakDir = file("${params.macs3_output}/${params.deeptools_macs3_profile}")
  assert peakDir.exists() : "MACS3 peak directory not found: ${peakDir}"

  def bamRows = chipRecords.collect { rec ->
    def bam = resolveUniqueFile(bamDir, { f ->
      f.name.endsWith('.clean.bam') && (f.name == "${rec.sample_id}.clean.bam" || f.name.startsWith("${rec.sample_id}_"))
    }, "clean BAM", rec.sample_id)
    tuple(rec.sample_id, rec.condition, rec.replicate, bam)
  }

  def peakFiles = chipRecords.collect { rec ->
    def peak = file("${peakDir}/${rec.sample_id}_peaks.${params.deeptools_peak_ext}")
    assert peak.exists() : "MACS3 peak file not found for ${rec.sample_id}: ${peak}"
    peak
  }

  def sampleManifestTsv = (['sample\tcondition\treplicate'] + chipRecords.collect { rec ->
    "${rec.sample_id}\t${rec.condition}\t${rec.replicate}"
  }).join('\n')

  def masterPeakset = merge_master_peakset(Channel.fromList(peakFiles).collect())

  def countInputs = Channel
    .fromList(bamRows)
    .combine(masterPeakset.master_bed)
    .map { sample_id, condition, replicate, bam, masterBed ->
      tuple(sample_id, condition, replicate, bam, masterBed)
    }

  def readCounts = count_reads_in_peaks(countInputs)
  def scaling = compute_scaling_factors(readCounts.reads_tsv.collect())

  def scaledInputs = Channel
    .fromList(bamRows)
    .combine(scaling.scaling_factors_tsv)
    .map { sample_id, condition, replicate, bam, scaleTable ->
      tuple(sample_id, condition, replicate, bam, scaleTable)
    }

  def scaledBws = deeptools_scaled_bw(scaledInputs)
  def meanTracks = deeptools_group_mean_bw(scaledBws.scaled_bw.collect(), Channel.value(sampleManifestTsv))

  deeptools_matrix_plot(
    meanTracks.mean_bws.collect(),
    Channel.value(upBed),
    Channel.value(downBed),
    Channel.value(referenceCondition),
    Channel.value(treatmentCondition)
  )
}

def parseSamplesMaster(master) {
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
      records << [
        sample_id: rec.sample_id?.toString()?.trim(),
        condition: rec.condition?.toString()?.trim(),
        replicate: (rec.replicate?.toString()?.isInteger() ? rec.replicate.toInteger() : Integer.MAX_VALUE),
        library_type: rec.library_type?.toString()?.trim()?.toLowerCase(),
        is_control: rec.is_control?.toString()?.trim()?.toLowerCase(),
        enabled: rec.enabled?.toString()?.trim()?.toLowerCase()
      ]
    }
  }
  assert header : "samples_master header not found: ${master}"
  assert header.contains('sample_id') : "samples_master missing required column: sample_id"
  assert header.contains('condition') : "samples_master missing required column: condition"
  records.findAll { it.sample_id && it.condition }
}

def isEnabled(rec) {
  rec.enabled == null || rec.enabled == '' || rec.enabled == 'true'
}

def isControl(rec) {
  rec.is_control == 'true' || rec.library_type == 'input'
}

def resolveUniqueFile(dir, matcher, label, sampleId) {
  def hits = dir.listFiles()?.findAll { f -> f.isFile() && matcher(f) } ?: []
  if (hits.isEmpty()) {
    throw new IllegalArgumentException("No ${label} found for sample_id '${sampleId}' under: ${dir}")
  }
  if (hits.size() > 1) {
    throw new IllegalArgumentException("Multiple ${label} files matched sample_id '${sampleId}': ${hits*.name.join(', ')}")
  }
  file(hits[0].toString())
}
