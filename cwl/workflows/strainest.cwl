class: Workflow
cwlVersion: v1.0
id: strainest
label: Strainest
$namespaces:
  sbg: 'https://www.sevenbridges.com'
inputs:
  - id: reference_dir
    type: Directory
    doc: Directory containing the bowtie-indexed database
    'sbg:x': -795.4100341796875
    'sbg:y': -361
  - id: reference_basename
    type: string
    doc: bowtie indexed database basename
    'sbg:x': -931.6915893554688
    'sbg:y': -252.3832244873047
  - id: bowtie2_read2
    type: File
    'sbg:x': -898.4100341796875
    'sbg:y': -71
  - id: bowtie2_read1
    type: File
    'sbg:x': -771.4100341796875
    'sbg:y': 51
  - id: strainest_est_dgrp
    type: File
    doc: SNV file
    'sbg:x': -353
    'sbg:y': -287
  - id: strainest_est_output_dir_name
    type: string
    doc: Output directory
    'sbg:x': -170.3577117919922
    'sbg:y': -409.4101257324219
outputs:
  - id: strainest_est_output_dir
    outputSource:
      - strainest_est/strainest_est_output_dir
    type: Directory
    'sbg:x': 10.544271469116211
    'sbg:y': -123.31122589111328
steps:
  - id: samtools_view
    in:
      - id: input
        source: bowtie2/output
      - id: input_1
        default: reads.bam
    out:
      - id: output
    run: ../tools/samtools-view.cwl
    label: samtools-view
    'sbg:x': -579
    'sbg:y': 68
  - id: samtools_sort
    in:
      - id: input
        source: samtools_view/output
      - id: input_1
        default: reads.sorted.bam
    out:
      - id: output
    run: ../tools/samtools-sort.cwl
    label: samtools-sort
    'sbg:x': -498
    'sbg:y': -170
  - id: samtools_index
    in:
      - id: input
        source: samtools_sort/output
    out:
      - id: output
    run: ../tools/samtools-index.cwl
    label: samtools-index
    'sbg:x': -377
    'sbg:y': 25
  - id: bowtie2
    in:
      - id: bowtie2_read1
        source: bowtie2_read1
      - id: bowtie2_read2
        source: bowtie2_read2
      - id: bowtie2_output
        default: bowtie_output
      - id: reference_dir
        source: reference_dir
      - id: reference_basename
        source: reference_basename
    out:
      - id: output
    run: ../tools/bowtie2.cwl
    label: bowtie2
    'sbg:x': -703.40625
    'sbg:y': -129
  - id: strainest_est
    in:
      - id: strainest_est_dgrp
        source: strainest_est_dgrp
      - id: strainest_est_bam
        source: samtools_index/output
      - id: strainest_est_output_dir_name
        source: strainest_est_output_dir_name
    out:
      - id: strainest_est_output_dir
    run: ../tools/strainest-est.cwl
    label: Strainest-est
    'sbg:x': -145
    'sbg:y': -39
requirements: []
