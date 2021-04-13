# Flow Atlas

_an interactive explorer for single-cell flow cytometry data_

Flow Atlas is an interactive data explorer for `FCS` datasets and their accompanying FlowJo workspace files, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). This tool is built in Julia and allows biologists to interactively inter-operate between FlowJo and modern computational libraries like [`GigaSOM.jl`](https://github.com/LCSB-BioCore/GigaSOM.jl)

# Getting started

### Quick start

To install Flow Atlas you need Julia 1.5+. Follow [your platform specific instructions](https://julialang.org/downloads/platform/). In order to make the `julia` executable available from the command line / terminal follow the instructions on adding Julia to your `PATH`/enviroment.

[Open a command line / terminal in a folder](https://www.groovypost.com/howto/open-command-window-terminal-window-specific-folder-windows-mac-linux) of your choice where you want to install Flow Atlas in. Make sure Git is installed by simply typing `git`. If the command is not recognised then [download and install it](https://git-scm.com/downloads). Now you are ready to download and install Flow Atlas, type

```bash
git clone https://github.com/gszep/FlowAtlas.jl.git # create new folder and download
cd FlowAtlas # navigate to new folder
julia src/server.jl # launch Flow Atlas
```
~~The script will build a system image `build/cytometry.so` first time you launch it which may take up to 30mins. This only needs to be built once; subsequent executions of `julia server.jl` will open a browser window instantly. Updates to Flow Atlas require re-building `build/cytometry.so`; simply delete that file to re-build it.~~
