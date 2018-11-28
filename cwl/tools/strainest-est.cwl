class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: strainest_est
baseCommand:
  - strainest
  - est
inputs:
  - id: strainest_est_dgrp
    type: File
    inputBinding:
      position: 0
    label: dgrp
    doc: dgrp file containing the reference SNP profiles
  - id: strainest_est_bam
    type: File
    inputBinding:
      position: 1
    label: bam
    doc: Sorted bam file
    secondaryFiles:
      - .bai
  - id: strainest_est_output_dir_name
    type: string
    inputBinding:
      position: 2
      shellQuote: false
outputs:
  - id: strainest_est_output_dir
    type: Directory
    outputBinding:
      glob: $(inputs.strainest_est_output_dir_name)
label: Strainest-est
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: compmetagen/strainest
  - class: InlineJavascriptRequirement
