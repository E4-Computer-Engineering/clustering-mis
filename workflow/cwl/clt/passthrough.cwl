# clt/passthrough.cwl
cwlVersion: v1.2
class: CommandLineTool
baseCommand: [sh, -c]
inputs:
  exe:
    type: File
    inputBinding:
      position: 1
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.exe.basename)
arguments:
  - valueFrom: "chmod +x $(inputs.exe.path) && cp $(inputs.exe.path) $(inputs.exe.basename)"
