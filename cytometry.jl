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

# ╔═╡ 48d505a4-469b-11eb-35b5-397a9a5498c2
begin ######## convert one-hot encoding to string labels
# 	labels.cell, fcsdata.cell = 1:size(fcsdata,1), 1:size(fcsdata,1)
# 	populations = stack(labels,Not(:cell),:cell,variable_name=:label)
	
# 	filter!(:value => x->x, populations)
# 	select!(populations,[:cell,:label])
	
# 	categorical!(populations,:label,compress=true)
# 	sort!(populations,:cell)
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

	channels = names(fcsdata)
	channel = Observable(first(channels))
	fluorescence = Observable(false)
	
	populations = names(labels)
	colors = Dict([ name=> Observable("#EEEEEE") for name ∈ populations ])
	
	scene, layout = layoutscene(resolution = (500,500))
	ax = layout[1,1] = LAxis(scene, 
		title="Self-organised Map Embedding")

	empty!(scene.events.mousedrag.listeners)
	mouseevents = MakieLayout.addmouseevents!(ax.scene)
	
	for name ∈ populations
		scatter!( ax, embedding[labels[!,name],:], markersize=1,
			color=@lift( $fluorescence ? "#ABCDEF" : $(colors[name]) ) )
	end
		
	# @lift( colorview(RGBA, expressionColors(
	# 	  scaleNorm(fcsdata[:,$channel]), expressionPalette(10, alpha=1/4)
	# 		)))
	


	# leg = layout[1, end+1] = LLegend(scene, 
	# [line1, scat1, line2, scat2], ["True", "Measured", "True", "Measured"])

	#layout[1,2] = LColorbar(scene,sc,width=10,height=Relative(1))

	
# 	###################################################### bind input events	
# 	on(scene.events.keyboardbuttons) do button
		
# 		if ispressed(button, Keyboard._1)
# 			μ₁[] = mouseevents.obs.val.data
			
# 		elseif ispressed(button, Keyboard._2)
# 			μ₂[] = mouseevents.obs.val.data
# 		end
# 	end
	
    RecordEvents( scene, "output" )
	println("Rendering Done")
end

# ╔═╡ 7cdadea8-4715-11eb-220c-475f60a98543
begin
	JSServe.with_session() do session, request
		return DOM.div( style,
			
			########################## embedding scatterplot
			DOM.div(scene),
			
			######################################### population color interactions
			DOM.div( DOM.label("Populations"), DOM.br(),
			[ ( DOM.input(type="color",name=name,id=name,value=colors[name].val,
								
					onchange=js"update_obs(($colors)[name],this.value)"),
			    DOM.label(name), DOM.br() ) for name ∈ populations ]),
		
			
			######################################## channel intensity
			DOM.div( DOM.input(type="checkbox",
				onchange=js"update_obs($fluorescence,this.checked);"),
			
			DOM.label("Fluorescence Intensity"), DOM.select(DOM.option.(channels),
			onchange=js"update_obs($channel,this.options[this.selectedIndex].text)"),
			)
		)
	end
end

# ╔═╡ 1973cde6-4523-11eb-3739-f9a9f97f197e
md"""
#### Cluster Exploration
In progress....
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
# ╠═48d505a4-469b-11eb-35b5-397a9a5498c2
# ╟─20596a96-4708-11eb-0b61-539eba64e3fd
# ╟─0b758ce0-4528-11eb-1df4-e97707ec4f1c
# ╠═aae5b128-436a-11eb-092b-0fc350961437
# ╠═7cdadea8-4715-11eb-220c-475f60a98543
# ╟─1973cde6-4523-11eb-3739-f9a9f97f197e
