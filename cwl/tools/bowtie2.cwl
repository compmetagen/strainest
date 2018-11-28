class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: bowtie2
baseCommand:
  - bowtie2
inputs:
  - id: bowtie2_read1
    type: File
    inputBinding:
      position: 3
      prefix: '-1'
  - id: bowtie2_read2
    type: File
    inputBinding:
      position: 4
      prefix: '-2'
  - id: bowtie2_output
    type: string
    inputBinding:
      position: 5
      prefix: '-S'
  - id: reference_dir
    type: Directory
  - id: reference_basename
    type: string
outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.bowtie2_output)
label: bowtie2
arguments:
  - position: 0
    prefix: ''
    valueFrom: '--very-fast'
  - position: 2
    prefix: '-x'
    valueFrom: $(inputs.reference_dir.path)/$(inputs.reference_basename)
  - position: 1
    prefix: ''
    valueFrom: '--no-unal'
requirements:
  - class: DockerRequirement
    dockerPull: compmetagen/strainest
  - class: InlineJavascriptRequirement
