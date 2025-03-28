cwlVersion: v1.2
class: CommandLineTool
requirements:
  ToolTimeLimit:
    timelimit: 300
arguments:
  - position: 5
    valueFrom: output.txt
inputs:
  annealing:
    type: File
    inputBinding:
      position: 4
  indices:
    type: File
    inputBinding:
      position: 3
  points:
    type: File
    inputBinding:
      position: 2
  silhouette:
    type: File
    inputBinding:
      position: 1
outputs:
  output:
    type: float
    outputBinding:
      glob: output.txt
      loadContents: true
      outputEval: $(parseFloat(self[0].contents))
