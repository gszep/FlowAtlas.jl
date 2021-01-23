### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ c49493f2-4366-11eb-2969-b3fff0c37e7b
begin
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()

	using Glob
	using AbstractPlotting.MakieLayout
	using WGLMakie,AbstractPlotting,LaTeXStrings

	AbstractPlotting.inline!(true)
	using JSServe,JSServe.DOM
	using JSServe: @js_str
	
	style = """

		font-family: "Alegreya Sans", sans-serif;
		font-weight: 400;
		font-size:   1em;
		color:       hsl(0, 0%, 25%);
		
		vertical-align: middle;
		horizontal-align: middle;
		display:        inline-block;
	
	"""

	using Plots: default
	using GraphRecipes

	include("lib/gates.jl")
	include("lib/utils.jl")

	using GigaSOM,Images,GigaScatter
	using MetaGraphs,DataFrames
	
	using StaticArrays,PolygonOps
	using LinearAlgebra,KernelDensity
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
Files matching the `glob` string are imported with their corresponding metadata that live in the FlowJo Workspace `wsp` file. A biexponential-type transformation is applied to raw fluorescence ``\mathbf{x}``. For simplicity we choose `arcsinh` with some cofactor ``\beta=250`` that scales the argument like
"""

# ╔═╡ b80a14f4-4435-11eb-07bd-7fcc8ca10325
L"
\mathbf{x}\rightarrow \mathrm{arcsinh}(\mathbf{x}/\beta)
"

# ╔═╡ 8110ec1e-54a9-11eb-05fb-d7e5281f6236
md"Since each `fcs` file may have different channel names, the user should provide a `channelMap` to map the channels into a common set of names so that the data can be stored in `DataFrame` types"

# ╔═╡ b7ff517c-450c-11eb-2054-61f8572bbccf
begin
	
	##################### paths to workpace and data
	workspace = "data/workspace.wsp"
	files = @glob_str("data/*/*.cleaned.fcs")
	
	#################### cofactor for biexp transformation
	cofactor = 250
	
	#################### drop intermediate gates from strategy
	selectLabels = Not([
			"CD4","CD8","Memory","Th17 | Th22",
			"CD127- CD25+","non-Tregs"
	])

	#################### channel map in case FCS files do not share channel names
	channelMap = Dict([

		"FJComp-355 379_28-A"=>"CD3", 
		"FJComp-355 560_40-A"=>"CD8", 

		"FJComp-355 820_60-A"=>"CD4",
		"FJComp-355 670_30-A"=>"CD4",

		"FJComp-640 780_60-A"=>"CCR7",
		"FJComp-405 780_60-A"=>"CD45RA", 

		"FJComp-561 780_60-A"=>"CD127", 
		"FJComp-640 670_30-A"=>"CD25", 

		"FJComp-561 610_20-A"=>"Helios", 
		"FJComp-561 585_15-A"=>"Foxp3", 
		"Foxp3-IgM"=>"Foxp3",

		"FJComp-405 710_40-A"=>"PD-1", 
		"FJComp-640 730_35-A"=>"CXCR5", 

		"FJComp-405 670_30-A"=>"CCR6", 
		"FJComp-488 715_30-A"=>"CXCR3", 

		"FJComp-405 605_40-A"=>"CCR4", 
		"FJComp-488 525_50-A"=>"CCR10", 

		"FJComp-405 450_50-A"=>"CD103", 
		"FJComp-355 740_35-A"=>"CD69",
		"FJComp-405 515_20-A"=>"HLA-DR"
	])
end

# ╔═╡ 477abee2-4367-11eb-003d-792fed6546ca
begin
	fcsdata,labels,groups,gating = load( files;
		workspace=workspace, cofactor=cofactor, channelMap=channelMap)

	select!( labels, selectLabels)
	nothing
end

# ╔═╡ 12e279d2-4477-11eb-0f4e-510c1329a935
begin
	using MetaGraphs: nv
	default(size=(500,500))
	exampleStrategy = last(first(gating))

	graphplot( exampleStrategy, method=:tree, nodeshape=:rect,
		names=[ get_prop(exampleStrategy,k,:name) for k ∈ 1:nv(exampleStrategy) ],
		fontsize=6, nodesize=0.12, linewidth=3, markercolor=:lightblue)
end

# ╔═╡ 10458cda-450c-11eb-2a3f-c168e35db194
md"""
The gating strategy is stored as a directed graph with each vertex having `:polygon` attribute which stores the ``\mathrm{arcsinh}`` transformed vertecies of the gating polygons in the chosen `:channels` and `:name` storing the name of each sub/population.
"""

# ╔═╡ 199c6556-4525-11eb-154b-034d1b0e0692
md"""
The gating graph is applied to the imported data `fcsdata` and we obtain a boolean dataframe `labels` that tells us whether a cell has been labelled by the gate and its parents or not. We can drop columns representing intermediate gates by passing them to `selectLabels`
"""

# ╔═╡ 232bb380-5439-11eb-0bf0-517df30fd027
md"Group information is also parsed from the workspace file to the variable `groups` . The user should use FlowJo sample groups to specify patient ids and condition names. These are used in interactive plots"

# ╔═╡ 0b758ce0-4528-11eb-1df4-e97707ec4f1c
md"""
#### Cluster Exploration
We can explore the clusters in the embedded space given by `EmbedSOM`. **Select colours for each population label** from the right-hand pannel or colour cells by fluoresence intensity of a **selected channel from the dropwon menu**. Verify the identity of known populations and discover novel ones.
"""

# ╔═╡ 20596a96-4708-11eb-0b61-539eba64e3fd
begin
	clusters,embedding = som(fcsdata,path="data/workspace.som")
	nothing
end

# ╔═╡ 1811f6f6-5439-11eb-33a5-11a16ce5ce76
md"""
**Note:**  Calculating the embedding for the first time takes a couple of minutes. Once, calculated it will be stored in as a `.som` file. The scatter plot will display once the cell above has finished executing. Feel free to drink more tea!
"""

# ╔═╡ aae5b128-436a-11eb-092b-0fc350961437
begin ###################################################### scene construction

	channels, populations, conditions = names(fcsdata),names(labels),names(groups)
	channel = Observable(first(channels))
	fluorescence = Observable(true)
	
	doneCalculations = Observable(true)
	compareGates = JSServe.Button("Compare Gates",style=style)
	
	markersize = JSServe.Slider(range(1,20,length=100))
	markeralpha = JSServe.Slider(range(0,1,length=100))
	
	markersize.value[] = 1.0
	markeralpha.value[] = 1.0
	
	########################################### transformation to unit interval
	function channelTransform(x::AbstractArray)
		return ( x .- minimum(x) )/( maximum(x) - minimum(x) )
	end
	
	scaled = combine(fcsdata,[ col => channelTransform => col
			for col ∈ names(fcsdata) ] )
	
	################################################## colormaps
	colors = Dict([ name=> Observable("#EEEEEE00") for name ∈ populations ])	
	fluorescenceMap = Dict([
			
		name => @lift( colorview(RGBA, expressionColors(	
			scaled[labels[!,name],$channel],
						
			$fluorescence ? expressionPalette(10,alpha=$markeralpha) : 
				[parse(RGBA,$(colors[name]))]
		)))
	for name ∈ populations ])
	
	####################################################
	embedScatter, embedScatterLayout = layoutscene(resolution = (500,500))
	embedScatterAx = embedScatterLayout[1,1] = AbstractPlotting.Axis(embedScatter,
		title="Self-organised Map Embedding")
	
	#################################################### scatterplot	
	for name ∈ populations#, condition ∈ conditions
		mask = labels[!,name] #.& groups[!,condition]

		scatter!( embedScatterAx, embedding[mask,:],
			color=fluorescenceMap[name],
			markersize=markersize)
	end
	
	#################################################### polygon selection
	left = Gate([
		MVector(0.0,0.0),MVector(1.5,1.5),MVector(3,1.5),
		MVector(2.0,-1.0),MVector(1.0,-1.8),MVector(0.0,-2.0)
	])

	right = Gate([
		MVector(6.0,-1.5),MVector(7.0,0.0),MVector(5.0,1.0),
		MVector(3,1.5),MVector(2.0,-1.0),MVector(4.0,-2.0)
	])
	
	gates = [left,right]

	lines!(embedScatterAx, left.polygon,
		linewidth=@lift( $(left.selected) & ($doneCalculations) ? 3 : 1 ),
		color=@lift( $doneCalculations ? RGBA(1/2,1,1/2,1) : RGBA(0,0,0,0.2) ))
	
	lines!(embedScatterAx, right.polygon,
		linewidth=@lift( $(right.selected) & ($doneCalculations) ? 3 : 1 ),
		color=@lift( $doneCalculations ? RGBA(1,1/4,0,1) : RGBA(0,0,0,0.2) ))
	
	downDisplacement = Observable(mouseposition(embedScatterAx.scene))
	vertexIdx = Observable(0)
	
	################################################ polygon interactions
	on(embedScatterAx.scene.events.mousebuttons) do events
		for gate ∈ [left,right]
			
			if ispressed(embedScatterAx.scene, Mouse.left) #### selection
				gate.selected[] = inpolygon( mouseposition(embedScatterAx.scene),
					gate.polygon[]; in=true,on=true,out=false)
			end
		end
	end
	
	empty!(embedScatter.events.mousedrag.listeners)
	on(embedScatterAx.scene.events.mousedrag) do drag
		
		################################ detect selection
		idx = findfirst(x->x.selected[],gates)
		if ~isnothing(idx) & ispressed(embedScatterAx.scene,Mouse.left)
			gate = gates[idx]
			
			################################## move selected
			if drag==Mouse.down
				
				position = mouseposition(embedScatterAx.scene)
				downDisplacement[] = position - PolygonOps.centroid(gate.polygon[])
				
				selectionThreshold = 0.12*sqrt(abs(PolygonOps.area(gate.polygon[])))
				vertex = @. norm([position]-gate.polygon[].x) < selectionThreshold
				vertexIdx[] = isnothing(findfirst(vertex)) ? 0 : findfirst(vertex)
			
			elseif (drag==Mouse.pressed) & (vertexIdx[]==0)
				
				displacement = mouseposition(embedScatterAx.scene) - 
					PolygonOps.centroid(gate.polygon[]) -
					downDisplacement[]

				polygon = gate.polygon[] .+ [displacement]
				gate.polygon[] = ClosedVector(polygon[Not(end)])
				
			elseif (drag==Mouse.pressed) & (vertexIdx[]>0)
				
				gate.polygon[][vertexIdx[]] = mouseposition(embedScatterAx.scene)
				gate.polygon[] = gate.polygon[]
			end
		end
	end

	############################################### zoom constraint
	on(embedScatterAx.scene.events.scroll) do scroll
		if maximum(embedScatterAx.limits[].widths) > 24
			limits!(embedScatterAx,-12,12,-12,12)
		end
	end
	nothing
end

# ╔═╡ 7cdadea8-4715-11eb-220c-475f60a98543
begin
	JSServe.App() do session, request
		return DOM.div( style=style,
						
		DOM.div(
			DOM.label(style="font-weight: bold","Colour:"),
			
			DOM.label("Population"),
			DOM.input(type="checkbox",checked=true,
				onchange=js"update_obs($fluorescence,this.checked);"),

			DOM.label("Fluorescence Intensity"),
			DOM.select(DOM.option.(channels),
			onchange=js"update_obs($channel,this.options[this.selectedIndex].text)")
		),
			
		DOM.div( style="margin-top: 10px;",
			DOM.label(style="font-weight: bold","Markers:"),
			DOM.label("Size"), markersize,
			DOM.label("Opacity"), markeralpha,
		),
			
		DOM.div( style="margin-top: 15px;",
			DOM.label(style="font-weight: bold","Population Labelling:"),
			DOM.label("Automatic"),
			DOM.label("Manual")
		),

		########################## embedding scatterplot
		DOM.div(style=style,embedScatter),

		######################################### interactive legend
		DOM.div(style=style,
			DOM.label(style="font-weight: bold","Populations"), DOM.br(),
				
			DOM.div( style="margin-top: 10px;",
				[( DOM.input(type="color",name=name,id=name,value="#EEEEEE",
				   onchange=js"update_obs(($colors)[name],this.value)"),
				   DOM.label(name), DOM.br() ) for name ∈ populations ]),
				
			DOM.br(),
			compareGates
		))

# 		DOM.div(style=style,[
# 			(DOM.input(type="checkbox",checked=true,
# 			onchange=js"update_obs($fluorescence,this.checked);"),
# 			DOM.label(name) )
				
				
# 			for name in ["BM","Spleen"]
# 		]))
	end
end

# ╔═╡ 21033346-4a17-11eb-22fe-8b938ed61ece
channelRange = range(-2,7,length=50)

# ╔═╡ b2a63c7a-49f7-11eb-300d-b7e283639a39
begin
	selectedLeft = Observable(fcsdata[ map( (x,y)->inpolygon(SVector(x,y),
				left.polygon[]; in=true,on=false,out=false),
				embedding[:,1], embedding[:,2] ),:])
	
	selectedRight = Observable(fcsdata[ map( (x,y)->inpolygon(SVector(x,y),
				right.polygon[]; in=true,on=false,out=false),
				embedding[:,1], embedding[:,2] ),:])
	
	################################## update
	on(compareGates) do button
		doneCalculations[] = false

		selectedLeft[] =
			fcsdata[ map( (x,y)->inpolygon(SVector(x,y),
			left.polygon[]; in=true,on=false,out=false),
			embedding[:,1], embedding[:,2] ),:]

		selectedRight[] =
			fcsdata[ map( (x,y)->inpolygon(SVector(x,y),
			right.polygon[]; in=true,on=false,out=false),
			embedding[:,1], embedding[:,2] ),:]
	end
	
	#################################################### density estimation	
	density(x::AbstractVector) = kde(x,channelRange).density
	densitiesLeft = @lift(combine( $selectedLeft,
			[ x => density => x for x ∈ channels ]) )

	densitiesRight = @lift(combine( $selectedRight,
			[ x => density => x for x ∈ channels ]) )
	########## todo maybe not the most efficient way to select.......
	
	###################################################### violin plot
	violins,violinsLayout = layoutscene(resolution = (700,250))
	violinsAx = violinsLayout[1,1] = AbstractPlotting.Axis(violins,
		
		xgridvisible=false, ygridvisible=false, 
		xticks = (1:length(channels),channels) )

	violinsAx.xticklabelrotation = Float32(π/2)
	violinsAx.xticklabelalign = (:top,:center)

	####################################### per channel
	for (k,channel) ∈ enumerate(channels)

		########################## left split
		maxDensityLeft = @lift( 5/2*maximum(x->x[1], $densitiesLeft[:,channel]) )
		pointsLeft = @lift( map( Point,
			k .- $densitiesLeft[:,channel]/$maxDensityLeft, channelRange ) )

		violinMeshLeft = @lift(AbstractPlotting.triangle_mesh($pointsLeft))
		mesh!( violinsAx, pointsLeft; shading=false,
			color=@lift( $doneCalculations ? RGBA(0,1,0,0.2) : RGBA(0,0,0,0.2) ) )

		########################## right split
		maxDensityRight = @lift( 5/2*maximum(x->x[1], $densitiesRight[:,channel]) )
		pointsRight = @lift( map( Point,
			k .+ $densitiesRight[:,channel]/$maxDensityRight, channelRange ) )

		violinMeshRight = @lift(AbstractPlotting.triangle_mesh($pointsRight))
		mesh!( violinsAx, violinMeshRight; shading=false,
			color=@lift( $doneCalculations ? RGBA(1,1/4,0,0.2) : RGBA(0,0,0,0.2) ) )
	end
	
	on(densitiesRight) do event
		doneCalculations[] = true
	end

	#################################### remove mouse interactions
	empty!(violins.events.mousedrag.listeners)
	empty!(violins.events.scroll.listeners)
	violins
end

# ╔═╡ 9be70f5e-4a14-11eb-3635-bd9fef2ea09a
md"""
Fluorescence distributions of selected populations can be compared across conditions. The conditions are extracted from groups defined in the FlowJo workspace.
"""

# ╔═╡ 0867a38c-4a15-11eb-0583-21b96f8009c3
md"""
* use checkboxes next to population names for filtering in violins
* click drag on 6 vertex polygons
* select polygons

* dropdown menues for filter and group by (tissue,label,patient)
* displaying stacked violin plots
* category columns must be parsed from workspace groups
"""

# ╔═╡ 02ed06fa-54b4-11eb-3e51-2908afcd617f
md"""
#### Labelling Agreement
By looking at the co-occurance matrix or the confusion matrix we get an overview of the agreement of population labelling between the expert and an automated proeedure. Furthermore this can be used to measure the agreement between experts and consistency of the expert with themselves.
"""

# ╔═╡ af681848-4da2-11eb-0607-1d22b3763245
begin
	labels.cell, clusters.cell = 1:size(labels,1), 1:size(clusters,1)
	multilabels = stack(labels,Not(:cell),:cell,variable_name=:label)
	select!(labels,Not(:cell))
	
	filter!(:value => x->x, multilabels)
	select!(multilabels,[:cell,:label])
	
	categorical!(multilabels,:label,compress=true)
	sort!(multilabels,:cell)
	
	overlap = combine(
		groupby(leftjoin(multilabels,clusters,on=:cell),Not(:cell)),
		nrow => :cells)
	
	sort!(overlap,:label)
	overlap = unstack(overlap,:label,:cluster,:cells)
	#transform!(overlap, ["1","2"] .=> (x->x/sum(skipmissing(x))) .=> ["1","2"] )
	#select!(overlap,"2")
	# skipmissing(x)/sum(skipmissing(x))
end

# ╔═╡ Cell order:
# ╟─c49493f2-4366-11eb-2969-b3fff0c37e7b
# ╟─ba294c1c-43ab-11eb-0695-1147232dc62a
# ╟─9dbf3d0c-4432-11eb-1170-4f59655b7820
# ╟─b80a14f4-4435-11eb-07bd-7fcc8ca10325
# ╟─8110ec1e-54a9-11eb-05fb-d7e5281f6236
# ╠═b7ff517c-450c-11eb-2054-61f8572bbccf
# ╟─477abee2-4367-11eb-003d-792fed6546ca
# ╟─10458cda-450c-11eb-2a3f-c168e35db194
# ╟─12e279d2-4477-11eb-0f4e-510c1329a935
# ╟─199c6556-4525-11eb-154b-034d1b0e0692
# ╟─232bb380-5439-11eb-0bf0-517df30fd027
# ╟─0b758ce0-4528-11eb-1df4-e97707ec4f1c
# ╟─20596a96-4708-11eb-0b61-539eba64e3fd
# ╟─1811f6f6-5439-11eb-33a5-11a16ce5ce76
# ╟─aae5b128-436a-11eb-092b-0fc350961437
# ╟─7cdadea8-4715-11eb-220c-475f60a98543
# ╟─b2a63c7a-49f7-11eb-300d-b7e283639a39
# ╟─21033346-4a17-11eb-22fe-8b938ed61ece
# ╟─9be70f5e-4a14-11eb-3635-bd9fef2ea09a
# ╟─0867a38c-4a15-11eb-0583-21b96f8009c3
# ╟─02ed06fa-54b4-11eb-3e51-2908afcd617f
# ╠═af681848-4da2-11eb-0607-1d22b3763245
