### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ d04f55b2-64b8-11eb-1eac-53cb43f9b117
begin
	using JSServe
	Page()
end

# ╔═╡ c49493f2-4366-11eb-2969-b3fff0c37e7b
begin
	using WGLMakie
	using Plots: default
	using GraphRecipes

	style = JSServe.Asset(joinpath(@__DIR__, "assets/style.css"))
	include("lib/components.jl")

	include("lib/gates.jl")
	include("lib/utils.jl")

	using GigaSOM,Images,GigaScatter
	using MetaGraphs,DataFrames
	using Setfield,Glob
	
	using StaticArrays,PolygonOps
	using LinearAlgebra,KernelDensity
	using LaTeXStrings
end

# ╔═╡ ca622d78-64b8-11eb-3668-974053d5f9c2
begin
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
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
	files = @glob_str("data/403C/*.cleaned.fcs")
	
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
	clusters,embedding = embed(fcsdata,path="data/workspace.som",
		perplexity=300, maxIter=10000);                  nothing
end

# ╔═╡ 1811f6f6-5439-11eb-33a5-11a16ce5ce76
md"""
**Note:**  Calculating the embedding for the first time takes a couple of minutes. Once, calculated it will be stored in as a `.som` file. The scatter plot will display once the cell above has finished executing. Feel free to drink more tea!
"""

# ╔═╡ aae5b128-436a-11eb-092b-0fc350961437
begin ###################################################### scene construction

	doneCalculations = Observable(true)
	fluorescence, automatic = Observable(true), Observable(true)
	compareGates = JSServe.Button("Compare Gates")
	
	markersize = JSServe.Slider(range(1,20,length=100),value=1.0)
	markeralpha = JSServe.Slider(range(0,1,length=100),value=1.0)
	
	################################################## colormaps
	channelRange, nlevels = range(-2,7,length=50), 10
	toIndex(x::AbstractVector{<:Number})=toIndex(x,channelRange;nlevels=nlevels)
	colorIndex = combine(fcsdata,[ col => toIndex => col for col ∈ names(fcsdata) ] )
	
	segments = fill( parse(RGBA,"#EEEEEE00"), size(fcsdata,1) )
	palette = convert( Vector{RGBA{Normed{UInt8,8}}},
		cgrad(:curl,nlevels,categorical=true,rev=true).colors.colors)

	channel = Observable(first(names(fcsdata)))
	color = Observable(palette[ colorIndex[!,channel[]] ])
	
	legend = convert( Vector{RGB},
		cgrad(:Accent_8,size(labels,2),categorical=true).colors.colors)

	legend = [ map( name -> Group(
		name,"#"*hex(popfirst!(legend))), names(labels) );
		map( name -> Group(name,true), names(groups) )
	]
	
	############################# filter codes
	codes = select( hcat(labels,groups), 
		AsTable(:) => ByRow(encode∘values) => "encoding")[:,:encoding]

	labelCode = encode([i ≤ size(labels,2) for i ∈ 1:length(legend)])
	groupCode = encode([i > size(labels,2) for i ∈ 1:length(legend)])
	
	############################# update on channel change
	on(channel) do channel
		if fluorescence[]
			
			selected = codes .& encode(map(x->x.selected[],legend))
			selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
			
			colors = palette[ colorIndex[!,channel] ]
			colors[.~selected] .= parse(RGBA,"#EEEEEE00")
			color[] = colors
		end
	end
	
	on(fluorescence) do fluorescence

		selected = codes .& encode(map(x->x.selected[],legend))
		selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
		
		colors = fluorescence ? palette[ colorIndex[!,channel[]] ] : segments
		colors[.~selected] .= parse(RGBA,"#EEEEEE00")
		color[] = colors
	end
	
	############################# update legend interactions
	for label ∈ legend
		
		mask=label.name ∈ names(labels) ? labels[!,label.name] : groups[!,label.name]
		if label.name ∈ names(labels) segments[mask] .= parse(RGBA,label.color[]) end

		on(label.color) do labelColor
			segments[mask] .= parse(RGBA,labelColor)
			
			selected = codes .& encode(map(x->x.selected[],legend))
			selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
			
			colors = fluorescence[] ? palette[ colorIndex[!,channel[]] ] : segments
			colors[.~selected] .= parse(RGBA,"#EEEEEE00")
			color[] = colors
		end
		
		on(label.selected) do selected
			if selected segments[mask] .= parse(RGBA,label.color[]) end
			
			selected = codes .& encode(map(x->x.selected[],legend))
			selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
			
			colors = fluorescence[] ? palette[ colorIndex[!,channel[]] ] : segments
			colors[.~selected] .= parse(RGBA,"#EEEEEE00")
			color[] = colors
		end
	end

	####################################################
	embedScatter = Figure(resolution = (510,500))
	embedScatterAx = embedScatter[1,1] = AbstractPlotting.Axis(embedScatter,
		title="Self-organised Map Embedding")
	
	#################################################### scatterplot	
	scatterPlot = scatter!(embedScatterAx,embedding,color=color,markersize=markersize)
	scatterPlot.colorrange = extrema(channelRange)
	embedScatter[:,2] = Colorbar(embedScatter,scatterPlot,
		colormap=cgrad(palette),width=10)
	
	#################################################### polygon selection	
	gates = [
			Gate([
				MVector(0.0,0.0),MVector(0.15,0.15),MVector(0.3,0.15),
				MVector(0.2,-0.1),MVector(0.1,-0.18),MVector(0.0,-0.2)
			],"#80FF80"),
		
			Gate([
				MVector(0.6,-0.15),MVector(0.7,0.0),MVector(0.5,0.1),
				MVector(0.3,0.15),MVector(0.2,-0.1),MVector(0.4,-0.2)
			],"#FF4000")
	]
	
	for gate ∈ gates
		
		lines!(embedScatterAx, gate.polygon,
			linewidth=@lift( $(gate.selected) & ($doneCalculations) ? 3 : 1 ),
			
			color=@lift( $doneCalculations ? parse(RGBA,$(gate.color)) : 
				RGBA(0,0,0,0.2)))
	end
	
	downDisplacement = Observable(mouseposition(embedScatterAx.scene))
	vertexIdx = Observable(0)
	
	################################################ polygon interactions
	on(embedScatterAx.scene.events.mousebuttons) do events
		for gate ∈ gates
			
			if ispressed(embedScatterAx.scene, Mouse.left) #### selection
				gate.selected[] = inpolygon( mouseposition(embedScatterAx.scene),
					gate.polygon[]; in=true,on=true,out=false)
			end
		end
	end
	
	empty!(embedScatter.scene.events.mousedrag.listeners)
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
	xmin,xmax = extrema(embedding[:,1])
	ymin,ymax = extrema(embedding[:,2])
	
	on(embedScatterAx.scene.events.scroll) do scroll
		if maximum(embedScatterAx.limits[].widths) > max(xmax-xmin,ymax-ymin)
			limits!(embedScatterAx,xmin,xmax,ymin,ymax)
		end
	end
	nothing
end

# ╔═╡ 7cdadea8-4715-11eb-220c-475f60a98543
begin
	JSServe.App() do session::Session
		return DOM.div( style,
						
		DOM.div( class="container",
			DOM.label(class="bold","Colour:"),
				
			DOM.label("Population"),
			DOM.label(class="switch",
					
				DOM.input(type="checkbox",checked=true,
				onchange=js"JSServe.update_obs($fluorescence,this.checked);"),
				DOM.span(class="slider round")
			),

			DOM.label("Fluorescence Intensity"),
			DOM.select(DOM.option.(names(fcsdata)),
			onchange=js"JSServe.update_obs($channel,this.options[this.selectedIndex].text)")
		),
			
		DOM.div( class="container",
			DOM.label(class="bold","Markers:"),
			DOM.label("Size"), markersize,
			DOM.label("Opacity"), markeralpha,
		),
			
		DOM.div( class="container",
			DOM.label(class="bold","Population Labelling:"),
				
			DOM.label("Gating Strategy"),
			DOM.label(class="switch",
					
				DOM.input(type="checkbox",checked=true,
				onchange=js"JSServe.update_obs($automatic,this.checked);"),
				DOM.span(class="slider round")
			),
			DOM.label("Self-organised Map")
		),

		########################## embedding scatterplot
		DOM.div( class="container",
			embedScatter.scene,

		######################################### interactive legend
			DOM.span( DOM.label(class="switch", style="width:70px",	
				DOM.input(type="checkbox", checked=true,
					onchange=js"""$(map(legend->legend.selected,
						filter( label->label.name ∈ names(labels), legend))
						).map(x=>JSServe.update_obs(x,this.checked))"""
				),
						
				DOM.span(class="slider text",style="font-weight:bold","Populations")),
				map( legend -> DOM.div( class="container",
							
					DOM.input(type="color",name=legend.name,id=legend.name,
						value=legend.color, onchange=js"""
							JSServe.update_obs($(legend.color),this.value)"""),
					DOM.label(class="switch", style="width:100%",
								
					DOM.input(type="checkbox",
						checked=legend.selected, onchange=js"""
							JSServe.update_obs($(legend.selected),this.checked)"""),
					DOM.span(class="slider text",legend.name))

				), filter( label->label.name ∈ names(labels), legend) ),
			),
				
			DOM.span( DOM.label(class="switch", style="width:70px",	
				DOM.input(type="checkbox", checked=true, 					
					onchange=js"""$(map(legend->legend.selected,
						filter( label->label.name ∈ names(groups), legend))
						).map(x=>JSServe.update_obs(x,this.checked))"""
				),
						
				DOM.span(class="slider text",style="font-weight:bold","Groups")),
				map( legend -> DOM.div( class="container",
							
					DOM.input(type="color",name=legend.name,id=legend.name,
						value=legend.color, onchange=js"""
							JSServe.update_obs($(legend.color),this.value)"""),
					DOM.label(class="switch", style="width:100%",
								
					DOM.input(type="checkbox",
						checked=legend.selected, onchange=js"""
							JSServe.update_obs($(legend.selected),this.checked)"""),
					DOM.span(class="slider text",legend.name))

				), filter( label->label.name ∈ names(groups), legend) ),
			),
		)
	)
	end
end

# ╔═╡ 628db7e6-6533-11eb-0afe-4973d1320a9f
compareGates

# ╔═╡ b2a63c7a-49f7-11eb-300d-b7e283639a39
begin
	#################################################### density estimation	
	density(x::AbstractVector) = kde(x,channelRange;
		bandwidth=step(channelRange)/2).density
	
	selected = codes .& encode(map(x->x.selected[],legend))
	selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
	
	selections = map( gate-> fcsdata[
		map( (x,y)->inpolygon(SVector(x,y),
						
			gate.polygon[]; in=true,on=false,out=false),
			embedding[:,1], embedding[:,2] ) .& selected,:], gates)
	
	densities = map( selection -> Observable( combine( selection,
		[ x => density => x for x ∈ names(fcsdata) ] )), selections )
	
	################################## update
	doneCalculations[] = true
	on(compareGates) do button
		
		doneCalculations[] = false
		for i ∈ 1:length(densities)
			
			selection = fcsdata[
				map( (x,y)->inpolygon(SVector(x,y),
						
				gates[i].polygon[]; in=true,on=false,out=false),
				embedding[:,1], embedding[:,2] ) .& selected,:]
			
			densities[i][] = combine( selection,
				[ x => density => x for x ∈ names(fcsdata) ] )
		end
	end
	########## todo maybe not the most efficient way to select.......
	
	###################################################### violin plot
	violins,violinsLayout = layoutscene(resolution = (700,250))
	violinsAx = violinsLayout[1,1] = AbstractPlotting.Axis(violins,
		
		xgridvisible=false, ygridvisible=false, 
		xticks = (1:size(fcsdata,2),names(fcsdata)) )

	violinsAx.xticklabelrotation = Float32(π/2)
	violinsAx.xticklabelalign = (:top,:center)

	####################################### per channel
	for (k,channel) ∈ enumerate(names(fcsdata))
		for i ∈ 1:length(densities)
			
			maxDensity = @lift( 5/2*maximum(x->x[1], $(densities[i])[:,channel]) )
			points=@lift(map(Point,
				k .+ (-1)^i * $(densities[i])[:,channel] / $maxDensity,
			channelRange))

			violinMesh = @lift(AbstractPlotting.triangle_mesh($points))
			mesh!( violinsAx, violinMesh; shading=false,
				
				color=@lift( $doneCalculations ? parse(RGBA,$(gates[i].color)) : 
					RGBA(0,0,0,0.2)))
		end
	end
	
	on(densities[end]) do event
		doneCalculations[] = true
	end

	#################################### remove mouse interactions
	empty!(violins.events.mousedrag.listeners)
	empty!(violins.events.scroll.listeners)
	violins
end

# ╔═╡ c010c502-6538-11eb-3f14-c9ea81b48b74
md"""
* gate buttons for both gates
* select population count to divide stats through. frequency by group or population
* Visualise confusion matrix
"""

# ╔═╡ 80238db4-6533-11eb-2db5-db6650a2af55
begin
	leftMask = map( (x,y)->inpolygon(SVector(x,y),
			left.polygon[]; in=true,on=false,out=false),
			embedding[:,1], embedding[:,2] )
	
	combine( groups[leftMask,:], names(groups) .=> sum .=> names(groups))
end

# ╔═╡ 2507bf12-653e-11eb-193e-2fef057dd4ec
begin
	combine( labels[leftMask,:], names(labels) .=> (x->sum(x)/564) .=> names(labels))
end

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
# ╟─ca622d78-64b8-11eb-3668-974053d5f9c2
# ╟─d04f55b2-64b8-11eb-1eac-53cb43f9b117
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
# ╠═20596a96-4708-11eb-0b61-539eba64e3fd
# ╟─1811f6f6-5439-11eb-33a5-11a16ce5ce76
# ╠═aae5b128-436a-11eb-092b-0fc350961437
# ╟─7cdadea8-4715-11eb-220c-475f60a98543
# ╟─628db7e6-6533-11eb-0afe-4973d1320a9f
# ╟─b2a63c7a-49f7-11eb-300d-b7e283639a39
# ╠═c010c502-6538-11eb-3f14-c9ea81b48b74
# ╟─2507bf12-653e-11eb-193e-2fef057dd4ec
# ╟─80238db4-6533-11eb-2db5-db6650a2af55
# ╟─02ed06fa-54b4-11eb-3e51-2908afcd617f
# ╠═af681848-4da2-11eb-0607-1d22b3763245
