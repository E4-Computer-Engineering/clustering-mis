cwlVersion: v1.2
class: CommandLineTool
requirements:
  ToolTimeLimit:
    timelimit: 300
baseCommand: [mpirun, --bind-to, core:overload-allowed]
arguments:
  - position: 4
    valueFrom: output.txt
  - position: 5
    valueFrom: indices.txt
inputs:
  clustering:
    type: File
    inputBinding:
      position: 2
  points:
    type: File
    inputBinding:
      position: 3
  processes:
    type: int
    inputBinding:
      position: 1
      prefix: -n
outputs:
  indices:
    type: File
    outputBinding:
      glob: indices.txt
  output:
    type: File
    outputBinding:
      glob: output.txt
