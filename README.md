# Flow Atlas

_an interactive explorer for single-cell flow cytometry data_

Flow Atlas is an interactive data explorer for `FCS` datasets and their accompanying FlowJo workspace files, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). This tool is built in Julia and allows biologists to interactively inter-operate between FlowJo and modern computational libraries like [`GigaSOM.jl`](https://github.com/LCSB-BioCore/GigaSOM.jl)

[![Build Status](https://travis-ci.com/gszep/FlowAtlas.jl.svg?branch=master)](https://travis-ci.com/gszep/FlowAtlas.jl)
[![Coverage](https://codecov.io/gh/gszep/FlowAtlas.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gszep/FlowAtlas.jl)

## Installation

To install Flow Atlas you need Julia 1.5+. Follow [your platform specific instructions](https://julialang.org/downloads/platform/). In order to make the `julia` executable available from the command line / terminal follow the instructions on adding Julia to your `PATH`/enviroment. Open a julia shell then type `] add FlowAtlas` and then hit ⏎ Return at the REPL. You should see `pkg> add FlowAtlas`

## Basic Usage
> :warning: **FCS files under a workspace must have unique names**. This limitation will be removed in future versions

The `FlowAtlas.run` method runs the web application in your default broswer
### Positional Arguments

### Keyword Arguments

```julia
using FlowAtlas

############################# paths to workspace and fcs files
workspace = "workspace.wsp"
files = glob"workspace/*.fcs"

################################ channel map for concatinating dataframes
channelMap = Dict([

    "FJComp-355 820_60-A"=>"CD4",
    "FJComp-355 670_30-A"=>"CD4",
    "Foxp3"=>"Foxp3-IgM",
    ...
])

################################ list of intermediate gate labels to drop
drop = ["CD4","CD3","CD4 | Memory","Th17-Th22","non-Tregs","non-B", ... ]

################################ launch webapp
FlowAtlas.run(workspace, files; channelMap=channelMap, drop=drop)
```