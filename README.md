# Hybrid Classical-Quantum clustering aggregation

This work is based upon the approach from [A clustering aggregation algorithm on neutral-atoms and annealing quantum processors](https://arxiv.org/pdf/2412.07558).

## TODO
- [ ] Add license
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
