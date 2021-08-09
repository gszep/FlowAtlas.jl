# Flow Atlas

_an interactive explorer for flow cytometry data_

Flow Atlas is an interactive data explorer for `FCS` datasets and their accompanying FlowJo workspace files, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). This tool is built in Julia and allows biologists to interactively inter-operate between FlowJo and modern computational libraries like [`GigaSOM.jl`](https://github.com/LCSB-BioCore/GigaSOM.jl)

[![Build Status](https://travis-ci.com/gszep/FlowAtlas.jl.svg?branch=master)](https://travis-ci.com/gszep/FlowAtlas.jl)
[![Coverage](https://codecov.io/gh/gszep/FlowAtlas.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gszep/FlowAtlas.jl)

## Installation

To install Flow Atlas, Julia 1.5+ is required. Follow [platform specific instructions](https://julialang.org/downloads/platform/). In order to make the `julia` executable available from the command prompt / terminal follow the instructions on adding Julia to the `PATH`/enviroment; make sure this has been done according to the instructions linked above.

Once julia has been added to the environment, open a command prompt / terminal and type `julia`. This opens a julis shell; you should see `julia>`. Then type `]` to enter the package manager, and type `add FlowAtlas` and then hit Return ⏎. You should see `pkg> add FlowAtlas`.

This installs the most recent version of Flow Atlas and makes it availalbe in your global Julia enviroment. Finally to exit the package manager hit Backspace ⌫; you should see `julia>` again.

## Update
Once inside a julia shell, type `]` to enter the package manager and then type `update FlowAtlas` then hit Return ⏎

## Basic Usage
> :warning: **FCS files under a workspace must have unique names**. This limitation will be removed in future versions

Open a julia shell in the directory that contains your workspace file and import Flow Altas using `using FlowAtlas`. The `FlowAtlas.run` method runs the web application in your default broswer. To launch the web application use the `run` method
```julia
using FlowAtlas
FlowAtlas.run( "workspace.wsp"; files="workspace/*.fcs")
```
### Keyword Arguments
