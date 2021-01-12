# Flow Explorer

_an interactive explorer for single-cell flow cytometry data_

Flow Explorer is an interactive data explorer for `FCS` datasets and their accompanying FlowJo workspace files, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). This tool is built in Julia and allows biologists to interactively inter-operate between FlowJo and modern computational libraries like `GigaSOM.jl`

# Getting started

### Quick start

To install Flow Explorer you need Julia 1.5+. Follow [your platform specific instructions](https://julialang.org/downloads/platform/). In order to make the `julia` executable available from the command line / terminal follow the instructions on adding Julia to your `PATH`/enviroment.

[Open a command line / terminal in the folder](https://www.groovypost.com/howto/open-command-window-terminal-window-specific-folder-windows-mac-linux) that you want to install Flow Explorer in, then type

```bash
git clone https://github.com/gszep/flow-explorer.git
cd flow-explorer
julia run.jl
```