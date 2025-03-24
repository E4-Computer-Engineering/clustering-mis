# Hybrid Classical-Quantum clustering aggregation

This work is based upon the approach from [A clustering aggregation algorithm on neutral-atoms and annealing quantum processors](https://arxiv.org/pdf/2412.07558).

## How to run

Make sure to have a working MPI installation available. Its include path should either be added to `$INCLUDE` or `$MPI_INC`.

The code can be compiled using `make`. Relevant executables will be under the build/bin directory.

### Clustering

The clustering executable can be run as follows:

```bash
mpirun -n 3 build/bin/clustering data/input/cluster_points_article.csv
```

You can optionally add another argument to save the output matrix to file:

```bash
mpirun -n 3 build/bin/clustering data/input/cluster_points_article.csv example_output.txt
```

You can add one more optional argument to save the indices of points that comprise each cluster:

```bash
mpirun -n 3 build/bin/clustering data/input/cluster_points_article.csv example_output.txt cluster_indices.txt
```

### Silhouette

The Silhouette score computation executable can be run as follows:

```bash
build/bin/silhouette data/input/cluster_points_article.csv cluster_indices.txt quantum_job_output.txt
```

You can add another argument if you want to save the score to a specific file, e.g.:

```bash
build/bin/silhouette data/input/cluster_points_article.csv cluster_indices.txt quantum_job_output.txt s-score.txt
```

## Expected output

### Clustering

Running the `clustering` executable will create an overlap matrix in the following form:

```
-1 8 8 8 0 0 0 0
0 -1 0 8 0 0 0 0
0 0 -1 8 8 0 8 0
0 0 0 -1 0 8 0 8
0 0 0 0 -1 0 8 0
0 0 0 0 0 -1 0 8
0 0 0 0 0 0 -1 8
0 0 0 0 0 0 0 -1
```

Each column/row represent a possible cluster. The diagonal terms are equal to -1, the off-diagonal ones are either 0 or a positive integer $\lambda$. Positive values denote overlaps between clusters. The value of $\lambda$ is defined as the number of different clusters, in this case 8, in order to prevent the selection of overlapping clusters.

If you choose to also save the points of each cluster, they will be in this form:

```
0,1,3,4
2,5,7
6,8,9
```

Each line corresponds to a different cluster. Each of its comma-separated values corresponds to a point from the original input file.

### Silhouette

Running the `silhouette` executable will compute and output the Silhouette score of a given clustering as a decimal number between -1 and 1. 

## Workflow run

It is also possible to run the whole Classical-Quantum pipeline (clustering + simulated annealing) as a workflow using the [StreamFlow](https://streamflow.di.unito.it) WMS. To do that, you need to clone this repository and all the included submodules, as follows:

```bash
git clone --recurse-submodules git@github.com:E4-Computer-Engineering/clustering-mis.git
```

The StreamFlow WMS requires Python 3.9 or newer. It can easily be installed as a Python package using the following commands:

```bash
python -m venv venv
source venv/bin/activate
pip install streamflow[report]==0.2.0.dev12
```

The workflow configuration is expressed in a declarative `streamflow.yml` file. An [example](workflow/streamflow.yml) targeting the [CINECA@Leonardo](https://leonardo-supercomputer.cineca.eu/) HPC facility is included in this repository. Modify it by adding your credentials (`username` and `sshKey`) and a path to a working directory in a shared portion of the Leonardo filesystem (e.g., in your `$HOME` folder).

At this point, simply run the workflow using this command:

```bash
streamflow run --name smart-hpc-qc workflow/streamflow.yml
```

When the workflow completes succesfully, you should find an `output.txt` file containing the results of the simulated annealing phase. In addition, the following command generates a report of the workflow run:

```bash
streamflow report --file workflow/streamflow.yml smart-hpc-qc
```

## TODO

- [ ] Add brief description with images

## Suggested dev setup

It recommended to use [VS Code](https://code.visualstudio.com/).

### Linting and autocompletion

IntelliSense from the Microsoft-provided C++ and Makefile extensions reports errors even if the code compiles.
It is recommended to use the [clangd extension](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd) instead.

Install the clangd extension and allow it to disable IntelliSense. Install the clangd language server if prompted.

Then install [bear](https://github.com/rizsotto/Bear) and, from the project root directory, run the following:

```bash
make clean; bear -- make
```

This will create a `compile_commands.json` file that is used by clangd to correctly inspect code.
Run "clangd: Restart language server" from the Command Palette (Ctrl+Shift+P) to read the newly created file.

You need to execute again the commands from above and restart the language server after each Makefile change

### Formatting

The clangd extension from the previous section can also format C/C++ code. Invoke it from Command Palette -> Format Document.
