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

	AbstractPlotting.inline!(true)
	using JSServe,JSServe.DOM
	using JSServe: @js_str
	
	style = DOM.style("""

		div { 
			font-family: "Alegreya Sans", sans-serif;
			font-weight: 400;
			font-size:   1em;
			color:       hsl(0, 0%, 25%);
		}
		
		div > * {
			vertical-align: middle;
			horizontal-align: middle;
		    display:        inline-block;
		}
	
	""",type="text/css")

	using Plots: default
	using GraphRecipes

	include("lib/utils.jl")
	using GigaSOM,Images,GigaScatter
	using MetaGraphs,DataFrames
	
	using StaticArrays,PolygonOps
	using KernelDensity
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

# ╔═╡ 12e279d2-4477-11eb-0f4e-510c1329a935
begin
	using MetaGraphs: nv
	default(size=(500,500))

	graphplot(gating,method=:tree,nodeshape=:rect,
		names=[ get_prop(gating,k,:name) for k ∈ 1:nv(gating)],
		fontsize=6,nodesize=0.12,linewidth=3,markercolor=:lightblue)
end

# ╔═╡ 10458cda-450c-11eb-2a3f-c168e35db194
md"""
The gating strategy is stored as a directed graph with each vertex having `:polygon` attribute which stores the ``\mathrm{arcsinh}`` transformed vertecies of the gating polygons in the chosen `:channels` and `:name` storing the name of each sub/population.
"""

# ╔═╡ 199c6556-4525-11eb-154b-034d1b0e0692
md"""
The gating graph is then applied to the imported data `fcsdata` and we obtain a boolean dataframe that tells us whether a cell has been labelled by the gate and its parents or not. We can drop columns representing intermediate gates with `select!` and only view final populations
"""

# ╔═╡ 81118ba6-4522-11eb-3de4-7f385efb6375
begin
	labels = gate(fcsdata,gating)
	select!(labels, Not(["CD4","CD8","Memory","Th17 | Th22",
				         "CD127- CD25+","non-Tregs"]) )
end

# ╔═╡ 20596a96-4708-11eb-0b61-539eba64e3fd
begin
	#################### train self-organised map
	som = initGigaSOM(fcsdata,20,20)
	som = trainGigaSOM(som,fcsdata)
	
	######################## extract clusters and embedding
	clusters = mapToGigaSOM(som,fcsdata) 
	embedding = embedGigaSOM(som,fcsdata)

	embedding .-= 10
	println("Training Self-organised Map : Done")
end

# ╔═╡ 0b758ce0-4528-11eb-1df4-e97707ec4f1c
md"""
#### Cluster Exploration
We can explore the clusters in the embedded space given by `EmbedSOM`. We can colour clusters by expression levels or labels to verify the identity of known clusters and discover novel clusters
"""

# ╔═╡ aae5b128-436a-11eb-092b-0fc350961437
begin ###################################################### scene construction

	channels, populations = names(fcsdata),names(labels)
	channel = Observable(first(channels))
	fluorescence = Observable(false)
	
	########################################### transformation to unit interval
	function channelTransform(x::AbstractArray)
		return ( x .- minimum(x) )/( maximum(x) - minimum(x) )
	end
	
	scaled = combine(fcsdata,[ col => channelTransform => col
			for col ∈ names(fcsdata) ] )
	
	################################################## colormaps
	colors = Dict([ name=> Observable("#EEEEEE") for name ∈ populations ])
	fluorescenceMap = Dict([
		name => @lift( colorview(RGBA, expressionColors(
						
			scaled[labels[!,name],$channel], $fluorescence ?
			expressionPalette(10,alpha=1/4) : 
						[parse(RGBA,$(colors[name]))]
		)))
	for name ∈ populations ])
	
	####################################################
	scene, layout = layoutscene(resolution = (500,500))
	ax = layout[1,1] = LAxis(scene,title="Self-organised Map Embedding")

	empty!(scene.events.mousedrag.listeners)
	mouseevents = MakieLayout.addmouseevents!(ax.scene)
	
	#################################################### scatterplot	
	for name ∈ populations
		scatter!( ax, embedding[labels[!,name],:],
			color=fluorescenceMap[name], markersize=1 )
	end
	
	#################################################### polygon selection
	polygon = Node([ SVector(0.0,0.0),SVector(5.0,0.0),SVector(5.0,5.0) ])
	closedPolygon = @lift([$polygon; [first($polygon)]])
	lines!(ax, closedPolygon, color=RGBA(1/2,1,1/2,1))
	
	on(ax.scene.events.mousebuttons) do buttons
	   if ispressed(ax.scene, Mouse.left)
		   polygon[] = push!(polygon[], mouseevents.obs.val.data)
	   end
	   return

	end
	nothing
end

# ╔═╡ 7cdadea8-4715-11eb-220c-475f60a98543
begin
	JSServe.with_session() do session, request
		return DOM.div( style, DOM.br(), DOM.br(),
			
		######################################## channel intensity
		DOM.div( DOM.input(type="checkbox",
			onchange=js"update_obs($fluorescence,this.checked);"),

		DOM.label("Fluorescence Intensity"), DOM.select(DOM.option.(channels),
		onchange=js"update_obs($channel,this.options[this.selectedIndex].text)")),

		########################## embedding scatterplot
		DOM.div(scene),

		######################################### population color interactions
		DOM.div( DOM.label("Populations"), DOM.br(),
		[ ( DOM.input(type="color",name=name,id=name,value=colors[name].val,

				onchange=js"update_obs(($colors)[name],this.value)"),
			DOM.label(name), DOM.br() ) for name ∈ populations ]),
		)
	end
end

# ╔═╡ 24a08ba4-4a15-11eb-1107-2dc357228e44
md"""
* constran zoom so that user does not fly off plot
* sliders for zoom, markersize, markeralpha
"""

# ╔═╡ 9be70f5e-4a14-11eb-3635-bd9fef2ea09a
md"""
Fluorescence distributions of selected populations can be compared across conditions. The conditions are extracted from groups defined in the FlowJo workspace.
"""

# ╔═╡ b2a63c7a-49f7-11eb-300d-b7e283639a39
begin
	#################################################### density estimation	
	selected = @lift(fcsdata[ map( (x,y)->inpolygon(SVector(x,y),
		$closedPolygon; in=true,on=false,out=false),
		embedding[:,1], embedding[:,2] ),:])
	
	channelRange = range(-2,7,length=50)
	density(x::AbstractVector) = map(SVector,kde(x,channelRange).density,channelRange)
	densities = @lift(combine( $selected, [ x => density => x for x ∈ channels ]) )

	################################################################# violin plot
	violins, violinsLayout = layoutscene(resolution = (650,300))
	violinsAx = violinsLayout[1,1] = LAxis(violins,title="Fluorescence")
	
	for (k,channel) ∈ enumerate(channels)
		maxDensity = @lift( maximum(x->x[1], $densities[!,channel]) )
		
		lines!(violinsAx, @lift(map( x-> x./SVector(2.5*$maxDensity,1)+SVector(k,0),
			$densities[!,channel]) ), color=RGBA(0,1,0,0.2))
		
		lines!(violinsAx, @lift(map( x-> x./SVector(-2.5*$maxDensity,1)+SVector(k,0), 
			$densities[!,channel])), color=RGBA(0,1,0,0.2))
	end
	
	violinsAx.xticks = (1:length(channels), channels)
	violinsAx.xticklabelrotation = Float32(π/2)
	violinsAx.xticklabelalign = (:top,:center)
	violins
end

# ╔═╡ 0867a38c-4a15-11eb-0583-21b96f8009c3
md"""
* dropdown menues for filter and group by (tissue,label,patient)
* displaying stacked violin plots
* category columns must be parsed from workspace groups
"""

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
# ╟─20596a96-4708-11eb-0b61-539eba64e3fd
# ╟─0b758ce0-4528-11eb-1df4-e97707ec4f1c
# ╟─aae5b128-436a-11eb-092b-0fc350961437
# ╟─7cdadea8-4715-11eb-220c-475f60a98543
# ╟─24a08ba4-4a15-11eb-1107-2dc357228e44
# ╟─9be70f5e-4a14-11eb-3635-bd9fef2ea09a
# ╟─b2a63c7a-49f7-11eb-300d-b7e283639a39
# ╟─0867a38c-4a15-11eb-0583-21b96f8009c3
