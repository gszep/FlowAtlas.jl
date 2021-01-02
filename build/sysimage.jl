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

sysimage_path="cytometry.so",
precompile_execution_file=["lib/utils.jl","cytometry.jl"])