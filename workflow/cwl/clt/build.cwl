cwlVersion: v1.2
class: CommandLineTool
requirements:
  InitialWorkDirRequirement:
    listing: $(inputs.src)
baseCommand: [make]
inputs:
  src:
    type:
      type: array
      items: [File, Directory]
  output_path: string
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_path)
  