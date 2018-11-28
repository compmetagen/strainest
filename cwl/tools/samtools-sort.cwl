class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: samtools_sort
baseCommand:
  - samtools
  - sort
inputs:
  - id: input
    type: File
    inputBinding:
      position: 1
  - id: input_1
    type: string
    inputBinding:
      position: 0
      prefix: '-o'
outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.input_1)
label: samtools-sort
requirements:
  - class: DockerRequirement
    dockerPull: compmetagen/strainest
  - class: InlineJavascriptRequirement
