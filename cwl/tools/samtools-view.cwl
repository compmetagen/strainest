class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: samtools_view
baseCommand:
  - samtools
  - view
inputs:
  - id: input
    type: File
    inputBinding:
      position: 2
  - id: input_1
    type: string?
    inputBinding:
      position: 1
      prefix: '-o'
outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.input_1)
label: samtools-view
arguments:
  - position: 0
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: '-b'
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: compmetagen/strainest
  - class: InlineJavascriptRequirement
