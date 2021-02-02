using PackageCompiler

create_sysimage([
    :AbstractPlotting,
    :DataFrames,
    :EzXML,
    :GigaSOM,
    :GigaScatter,
    :Glob,
    :GraphRecipes,
    :Images,
    :JSServe,
    :KernelDensity,
    :LaTeXStrings,
    :LightGraphs,
    :MetaGraphs,
    :Observables,
    :Plots,
    :PolygonOps,
    :Setfield,
    :StaticArrays,
    :TSne,
    :WGLMakie    
],

sysimage_path="build/cytometry.so",
precompile_execution_file=["lib/components.jl","lib/gates.jl","lib/utils.jl","cytometry.jl"])