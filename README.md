# Hybrid Classical-Quantum clustering aggregation

This work is based upon the approach from [A clustering aggregation algorithm on neutral-atoms and annealing quantum processors](https://arxiv.org/pdf/2412.07558).

## How to run
Make sure to have a working MPI installation available. Its include path should either be added to $INCLUDE or $MPI_INC.

The code can be compiled using `make`. The newly built executable will be under the build/bin directory.

The code can be run as follows:
```bash
mpirun -n 3 build/bin/clustering data/input/cluster_points_article.csv
```
You can optionally add another argument to save the output matrix to file:
```bash
mpirun -n 3 build/bin/clustering data/input/cluster_points_article.csv example_output.txt
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
