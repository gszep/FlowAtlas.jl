using PackageCompiler

create_sysimage([
    :AbstractPlotting,
    :DataFrames,
    :EzXML,
    :GigaSOM,
    :GigaScatter,
    :GraphRecipes,
    :Images,
    :JSServe,
    :KernelDensity,
    :LaTeXStrings,
    :LightGraphs,
    :MetaGraphs,
    :Plots,
    :PolygonOps,
    :StaticArrays,
    :WGLMakie    
    ],

sysimage_path="build/cytometry.so",
precompile_execution_file=["lib/gates.jl","lib/utils.jl","cytometry.jl"])