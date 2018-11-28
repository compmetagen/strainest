class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: samtools_index
baseCommand:
  - samtools
  - index
inputs:
  - id: input
    type: File
    inputBinding:
      position: 0
outputs:
  - id: output
    type: File
    secondaryFiles: .bai
    outputBinding:
      glob: $(inputs.input.basename)
label: samtools-index
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.input)
  - class: InlineJavascriptRequirement
