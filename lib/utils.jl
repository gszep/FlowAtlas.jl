using GigaSOM,EzXML,DataFrames
using LightGraphs,MetaGraphs
using PolygonOps,StaticArrays

using DataFrames: CategoricalValue
import Base: show # display categorical types neatly
show(io::IO, ::MIME"text/html", x::CategoricalValue) = print(io, get(x))

function load(path::String; gating=nothing, cofactor=250, kwargs...)

    params, data = loadFCS(path; kwargs...)
    params = getMetaData(params)

    data = DataFrame(data, "S" ∈ names(params) ? params.S : params.N)
    @. data = asinh(data/cofactor)

    if isnothing(gating)
        return data

    else 
		return data, gatingGraph(path,gating; params= "S" ∈ names(params) ? params : nothing,
			cofactor=cofactor)
    end
end

function gatingGraph(path::String,gating::String; params=nothing,cofactor=250)

	channelmap = isnothing(params) ? nothing : Dict([ N=>S for (N,S) ∈ zip(params.N,params.S) ])
    workspace = root(readxml(gating))

	############################################## population names
	populations = findall("//DataSet[contains(@uri,'$path')]/..//Population",workspace)
	@assert( length(populations)>0, "gating not found in workspace $gating for sample\n$path")

	graph = MetaDiGraph{Int64,Bool}() ############# store strategy in graph
	set_props!(graph,Dict(:sample=>path))

	for population ∈ populations
		Gate,name = findfirst("Gate",population), population["name"]
		
		id = Gate["gating:id"]
		parent_id = haskey(Gate,"gating:parent_id") ? Gate["gating:parent_id"] : nothing
		
		#################### iterate through disconnected polygons
		for gate ∈ eachelement(Gate)
			
			channels = map( dimension-> isnothing(params) ? dimension["data-type:name"] : channelmap[dimension["data-type:name"]],
				findall("gating:dimension/data-type:fcs-dimension",gate) )
			
			################################# add polygons to graph
			if ~occursin("threshold",name)
				@assert(gate.name=="PolygonGate", "$(gate.name) not supported for population label $name in sample\n$path")
				
				vertices = map( coordinate-> parse(Float32,coordinate["data-type:value"]),
					findall("gating:vertex/gating:coordinate",gate) )
				
				polygon = map( (x,y)->SVector(asinh(x/cofactor),asinh(y/cofactor)), @view(vertices[1:2:end]), @view(vertices[2:2:end]) )
				push!(polygon,first(polygon))
				
				add_vertex!(graph) ################## store gate as vertex
				set_indexing_prop!( graph,nv(graph), :id,id )
				set_props!( graph,nv(graph), Dict(:channels=>channels,:polygon=>polygon,:name=>name) )
                
                ################# connect parent gates to children
				~isnothing(parent_id) && add_edge!(graph,graph[parent_id,:id],nv(graph))
			end
		end
    end
    
    return graph
end


function gate(graph::MetaDiGraph,idx::Integer;prefix="__gate__")
	return get_prop(graph,idx,:channels) =>
		ByRow( (x,y)->inpolygon( SVector(x,y), get_prop(graph,idx,:polygon);

			in=true, on=false, out=false)
	) => prefix*get_prop(graph,idx,:name)
end


function gate!(data::DataFrame,gating::MetaDiGraph,idx::Integer;prefix="__gate__")
	name = prefix*get_prop(gating,idx,:name)
	transform!( data, gate(gating,idx;prefix=prefix) )

	for parent ∈ inneighbors(gating,idx) ######### todo @gszep does not support multiple parents
		parent = prefix*get_prop(gating,parent,:name)
		transform!(data, [name,parent] => ByRow((x,y)->x&y) => name)
	end
	
	for child ∈ outneighbors(gating,idx)
		gate!(data,gating,child;prefix=prefix)
	end
end


function gate(data::DataFrame,gating::MetaDiGraph;prefix="__gate__")

	roots = [ idx for idx ∈ 1:nv(gating) if indegree(gating,idx)==0 ]
	names = [ prefix*get_prop(gating,idx,:name) for idx ∈ 1:nv(gating)]

	for idx ∈ roots
		gate!(data,gating,idx;prefix=prefix)
	end
	
	labels = select(data,names)
	select!(data,Not(names))

	rename!(labels,Dict([ name=>replace(name,prefix=>"") for name ∈ names ]))
	return labels
end