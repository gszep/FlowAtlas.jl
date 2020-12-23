### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ c49493f2-4366-11eb-2969-b3fff0c37e7b
begin
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
		
	using AbstractPlotting.MakieLayout
	using WGLMakie,AbstractPlotting,LaTeXStrings
	using GraphRecipes,Plots
	
	AbstractPlotting.inline!(true)
	using PolygonOps,StaticArrays,LightGraphs,MetaGraphs,EzXML
	include("lib/utils.jl")
	
	using GigaSOM,GigaScatter
	using Images,DataFrames
end

# ╔═╡ ba294c1c-43ab-11eb-0695-1147232dc62a
md"""
## New Tools for Exploring Flow Cytometry Data
In this interactive article we will demonstrate flexible tools in Julia that allow 
us to interact with FCS files and FlowJo Workspaces. Features include but are not 
limited to clustering and dimensionality reduction with `GigaSOM.jl`, automated population label assignment, threshold detection and visualisation with dot plots and violin plots.
"""

# ╔═╡ 9dbf3d0c-4432-11eb-1170-4f59655b7820
md"""
#### Loading and Pre-processing
First we import an `fcs` file and its corresponding gating strategy that lives in
the FlowJo Workspace `wsp` file. Typically a biexponential-type transformation is applied to raw fluorescence values ``\mathbf{x}``. For simplicity we choose `arcsinh` with some cofactor ``\beta=250`` that scales the argument like
"""

# ╔═╡ b80a14f4-4435-11eb-07bd-7fcc8ca10325
L"
\mathbf{x}\rightarrow \mathrm{arcsinh}(\mathbf{x}/\beta)
"

# ╔═╡ b7ff517c-450c-11eb-2054-61f8572bbccf
begin
	cofactor = 250
	path = "data/403C/Donor_BM_024.cleaned.fcs"
	workspace = "data/workspace.wsp"
end

# ╔═╡ 477abee2-4367-11eb-003d-792fed6546ca
begin
	fcsdata,gating = load(path; gating=workspace, cofactor=cofactor)
	fcsdata
end

# ╔═╡ 10458cda-450c-11eb-2a3f-c168e35db194
md"""
The gating strategy is stored as a directed graph with each vertex having `:polygon` attribute which stores the ``\mathrm{arcsinh}`` transformed vertecies of the gating polygons in the chosen `:channels` and `:name` storing the name of each sub/population.
"""

# ╔═╡ 12e279d2-4477-11eb-0f4e-510c1329a935
begin	
	default(size=(2000, 2000))
	graphplot(gating,method=:tree,nodeshape=:rect,
		names=[ get_prop(gating,k,:name) for k ∈ 1:nv(gating)],
		fontsize=24,nodesize=0.01,linewidth=5)
end

# ╔═╡ 199c6556-4525-11eb-154b-034d1b0e0692
md"""
The gating graph is then applied to the imported data `fcsdata` and we obtain a boolean dataframe that tells us whether a cell has been labelled by the gate and its parents or not. We can drop columns representing intermediate gates with `select!` and only view final populations
"""

# ╔═╡ 81118ba6-4522-11eb-3de4-7f385efb6375
begin
	labels = gate(fcsdata,gating)
	select!(labels,
		Not(["CD4","CD8","Memory","Th17 | Th22","CD127- CD25+","non-Tregs"]))
end

# ╔═╡ 1973cde6-4523-11eb-3739-f9a9f97f197e
md"""
#### Cluster Exploration
In progress....
"""

# ╔═╡ aae5b128-436a-11eb-092b-0fc350961437
begin
	channels = ("CD45RA","CCR7")
	
	###################################################### scene construction
	scene, layout = layoutscene(resolution = (500,500))
	ax = layout[1,1] = LAxis(scene,
		
		xpanlock = true, xzoomlock = true,
		ypanlock = true, yzoomlock = true,
		
		xlabel=first(channels), ylabel=last(channels) )

	empty!(scene.events.mousedrag.listeners)
	mouseevents = MakieLayout.addmouseevents!(ax.scene)
	population = fcsdata[labels[!,"CD4 | EMRA"],:]

	AbstractPlotting.scatter!( ax,
		population[!,first(channels)], population[!,last(channels)],
		color=:darkred, markersize=1 )
	
	population = fcsdata[labels[!,"CD4 | Naive"],:]

	AbstractPlotting.scatter!( ax,
		population[!,first(channels)], population[!,last(channels)],
		color=:gold, markersize=1 )
	
	population = fcsdata[labels[!,"CD4 | EM"],:]

	AbstractPlotting.scatter!( ax,
		population[!,first(channels)], population[!,last(channels)],
		color=:lightblue, markersize=1 )
	
	population = fcsdata[labels[!,"CD4 | CM"],:]

	AbstractPlotting.scatter!( ax,
		population[!,first(channels)], population[!,last(channels)],
		color=:darkblue, markersize=1 )
	
# 	###################################################### bind input events	
# 	on(scene.events.keyboardbuttons) do button
		
# 		if ispressed(button, Keyboard._1)
# 			μ₁[] = mouseevents.obs.val.data
			
# 		elseif ispressed(button, Keyboard._2)
# 			μ₂[] = mouseevents.obs.val.data
# 		end
# 	end
	
    RecordEvents( scene, "output" )
    scene
end

# ╔═╡ 0b758ce0-4528-11eb-1df4-e97707ec4f1c
md"""
#### Dimensionality Reduction
Sample dimensionality reduction using `EmbedSOM` Colour by `CD4` channel
"""

# ╔═╡ 8b928400-43b4-11eb-1882-8b4d1ceff395
begin
	som = initGigaSOM(fcsdata, 20, 20)    # random initialization of the SOM codebook
	som = trainGigaSOM(som, fcsdata)      # SOM training
	clusters = mapToGigaSOM(som, fcsdata) # extraction of per-cell cluster IDs
	e = embedGigaSOM(som, fcsdata)        # EmbedSOM projection to 2D

	raster = rasterize( (900, 900), Matrix{Float64}(e'),
		expressionColors(
			
			scaleNorm(Array{Float64}(fcsdata[!,"CD4"])),
			expressionPalette(10, alpha=0.25)
		)
	)
	
	colorview(RGBA,rasterKernelCircle(1,raster))
end

# ╔═╡ Cell order:
# ╟─c49493f2-4366-11eb-2969-b3fff0c37e7b
# ╟─ba294c1c-43ab-11eb-0695-1147232dc62a
# ╟─9dbf3d0c-4432-11eb-1170-4f59655b7820
# ╟─b80a14f4-4435-11eb-07bd-7fcc8ca10325
# ╠═b7ff517c-450c-11eb-2054-61f8572bbccf
# ╟─477abee2-4367-11eb-003d-792fed6546ca
# ╟─10458cda-450c-11eb-2a3f-c168e35db194
# ╟─12e279d2-4477-11eb-0f4e-510c1329a935
# ╟─199c6556-4525-11eb-154b-034d1b0e0692
# ╠═81118ba6-4522-11eb-3de4-7f385efb6375
# ╟─1973cde6-4523-11eb-3739-f9a9f97f197e
# ╟─aae5b128-436a-11eb-092b-0fc350961437
# ╟─0b758ce0-4528-11eb-1df4-e97707ec4f1c
# ╟─8b928400-43b4-11eb-1882-8b4d1ceff395
