cwlVersion: v1.2
class: Workflow

$namespaces:
  cwltool: "http://commonwl.org/cwltool#"

requirements:
  StepInputExpressionRequirement: {}
  ToolTimeLimit:
    timelimit: 300
inputs:
  annealing:
    type:
      type: array
      items: [File, Directory]
  clustering:
    type:
      type: array
      items: [File, Directory]
  points: File
  processes: int?
  threshold: float?
  annealing_bin: File
  clustering_bin: File
  silhouette_bin: File
outputs:
  output:
    type: File
    outputSource: loop/output
  score:
    type: float
    outputSource: loop/score
steps:
  build-clustering:
    run: clt/passthrough.cwl
    in:
      exe: clustering_bin
    out: [output]
  build-annealing:
    run: clt/passthrough.cwl
    in:
      exe: annealing_bin
    out: [output]
  build-silhouette:
    run: clt/passthrough.cwl
    in:
      exe: silhouette_bin
    out: [output]
  loop:
    requirements:
      InlineJavascriptRequirement: {}
      cwltool:Loop:
        loopWhen: $(inputs.score < inputs.threshold)
        loop:
          score: score
          seed:
            valueFrom: $(inputs.seed + 1)
        outputMethod: last
    run:
      class: Workflow
      inputs:
        annealing: File
        clustering: File
        points: File
        processes:
          type: int
          default: 3
        seed: int
        silhouette: File
      outputs:
        output:
          type: File
          outputSource: annealing/output
        score:
          type: float
          outputSource: silhouette/output
      steps:
        clustering:
          run: clt/clustering.cwl
          in:
            clustering: clustering
            points: points
            processes: processes
            seed: seed
          out: [indices, output]
        annealing:
          run: clt/annealing.cwl
          in:
            annealing: annealing
            qubo: clustering/output
          out: [output]
        silhouette:
          run: clt/silhouette.cwl
          in:
            silhouette: silhouette
            points: points
            indices: clustering/indices
            annealing: annealing/output
          out: [output]
    in:
      annealing: build-annealing/output
      clustering: build-clustering/output
      points: points
      processes: processes
      seed:
        default: 0
      score:
        default: 0.0
      silhouette: build-silhouette/output
      threshold:
        default: 0.8
    out: [output, score]
