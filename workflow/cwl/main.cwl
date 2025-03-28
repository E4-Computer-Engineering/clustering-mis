cwlVersion: v1.2
class: Workflow

$namespaces:
  cwltool: "http://commonwl.org/cwltool#"

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
  processes: int?
  threshold: float?
outputs:
  output:
    type: File
    outputSource: loop/output
  score:
    type: float
    outputSource: loop/score
steps:
  build-clustering:
    run: clt/build.cwl
    in:
      src: clustering_src
      output_path:
        valueFrom: build/bin/clustering
    out: [output]
  build-annealing:
    run: clt/build.cwl
    in:
      src: annealing_src
      output_path:
        valueFrom: build/bin/simAnnSingle.out
    out: [output]
  build-silhouette:
    run: clt/build.cwl
    in:
      src: clustering_src
      output_path:
        valueFrom: build/bin/silhouette
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
        annealing_script: File
        clustering_script: File
        points: File
        processes:
          type: int
          default: 3
        seed: int
        silhouette_script: File
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
            clustering: clustering_script
            points: points
            processes: processes
            seed: seed
          out: [indices, output]
        annealing:
          run: clt/annealing.cwl
          in:
            annealing: annealing_script
            qubo: clustering/output
          out: [output]
        silhouette:
          run: clt/silhouette.cwl
          in:
            silhouette: silhouette_script
            points: points
            indices: clustering/indices
            annealing: annealing/output
          out: [output]
    in:
      annealing_script: build-annealing/output
      clustering_script: build-clustering/output
      points: points
      processes: processes
      seed:
        default: 0
      score:
        default: 0.0
      silhouette_script: build-silhouette/output
      threshold:
        default: 0.5
    out: [output, score]
