cwlVersion: v1.2
class: Workflow
requirements:
  StepInputExpressionRequirement: {}
  ToolTimeLimit:
    timelimit: 300
inputs:
  annealing_src:
    type:
      type: array
      items: [File, Directory]
  clustering_src:
    type:
      type: array
      items: [File, Directory]
  points: File
  processes:
    type: int
    default: 3
outputs:
  annealing:
    type: File
    outputSource: annealing/output
steps:
  build-clustering:
    run: clt/build.cwl
    in:
      src: clustering_src
      output_path:
        valueFrom: build/bin/clustering
    out: [output]
  clustering:
    run: clt/clustering.cwl
    in:
      clustering: build-clustering/output
      points: points
      processes: processes
    out: [indices, output]
  build-annealing:
    run: clt/build.cwl
    in:
      src: annealing_src
      output_path:
        valueFrom: build/bin/simAnnSingle.out
    out: [output]
  annealing:
    run: clt/annealing.cwl
    in:
      annealing: build-annealing/output
      qubo: clustering/output
    out: [output]
