cwlVersion: v1.2
class: CommandLineTool
requirements:
  ToolTimeLimit:
    timelimit: 300
arguments:
  - position: 3
    valueFrom: output.txt
inputs:
  annealing:
    type: File
    inputBinding:
      position: 1
  qubo:
    type: File
    inputBinding:
      position: 2
outputs:
  output:
    type: File
    outputBinding:
      glob: output.txt
