using GigaSOM,EzXML,DataFrames,TSne
using Serialization: serialize,deserialize
using Glob: GlobMatch

using DataFrames: CategoricalValue
import Base: show # display categorical types neatly
show(io::IO, ::MIME"text/html", x::CategoricalValue) = print(io, get(x))

function load(path::String; workspace::String="", cofactor::Number=250, channelMap::Dict=Dict(), kwargs...)

	#######################################
	path = replace(path,"\\"=>"/")
    params, data = loadFCS(path; kwargs...)
	params = getMetaData(params)
	
	###################################### data with human readable channel names
	data = DataFrame(data,params.N)
	function updateMap(laser::String,marker::String)

		if laser ∈ keys(channelMap)
			return channelMap[laser]

		elseif marker ∈ keys(channelMap)
			return channelMap[marker]

		else 
			return marker
		end
	end

	channelNames = Dict([ laser => updateMap(laser,marker)
		for (laser,marker) ∈ zip(params.N,"S" ∈ names(params) ? params.S : params.N) ])
	rename!(data,channelNames) 

	###################################### biexponential transformation
	@. data = asinh(data/cofactor)
	
	if isempty(workspace)
        return data,nothing,nothing,nothing

	else ################################ load metadata from workspace
		workspace = root(readxml(workspace))

		sampleID = findfirst("//DataSet[contains(@uri,'$path')]",workspace)["sampleID"]
		groups = DataFrame([ (group["name"]=>fill(true,size(data,1))) for group ∈ 
			findall("//SampleRefs/SampleRef[contains(@sampleID,'$sampleID')]/../..",workspace)
				if  group["name"] ≠ "All Samples" ])

		gating = gatingGraph(path,workspace;channelMap=channelNames,cofactor=cofactor)
		labels = gate(data,gating)

		return data,labels,groups,gating
    end
end


function load(pattern::GlobMatch; workspace::String="", cofactor::Number=250, channelMap::Dict=Dict(), kwargs...)

	data,labels,groups = DataFrame(),DataFrame(),DataFrame()
	gatings = Dict()
	
	for path ∈ readdir(pattern)
		fcs,label,group,gating = load(path;workspace=workspace,cofactor=cofactor,channelMap=channelMap,kwargs...)
		
		append!(data,fcs)
		gatings[path] = gating

		append!(labels,label,cols=:union)
		append!(groups,group,cols=:union)		
	end

	if isempty(workspace)
		return data,nothing,nothing,nothing
	end

	map( name->replace!(labels[!,name],missing=>false), names(labels) )
	map( name->replace!(groups[!,name],missing=>false), names(groups) )

	disallowmissing!(labels)
	disallowmissing!(groups)

	return data,labels, groups,gatings
end


function embed(data::DataFrame;path::String="",xdim::Int64=20,ydim::Int64=20,perplexity=10,maxIter=10000,eta=200.0)
	# Kratochvíl M, Koladiya A and Vondrášek J. Generalized EmbedSOM on quadtree-structured self-organizing maps. F1000Research 2020, 8:2120

	if isfile(path)
		som = deserialize(path)

	else
		som = initGigaSOM(data,xdim,ydim)
		som = trainGigaSOM(som,data)
		~isempty(path) && serialize(path,som)
	end

	########################### landmarks
	som.grid = tsne(som.codes,2,0,maxIter,perplexity;eta=eta)

	######################## extract clusters and embedding
	clusters = mapToGigaSOM(som,data)
	embedding = embedGigaSOM(som,data)

	rename!(clusters,"index"=>"cluster")
	return clusters,embedding
end


function encode( input; T::DataType=Int )
	v,encodedInt = 1,zero(T)
	iter = Iterators.reverse(input)

	for i ∈ iter
		encodedInt += v * i
		v <<= 1
	end
	return encodedInt
end


function toIndex(x::AbstractVector{<:Number},channelRange::AbstractVector{<:Number};nlevels::Int=10)

	min,max = extrema(channelRange)
	scaled = ( x .- min ) / ( max - min )
	
	@. scaled[scaled<0] = 0
	@. scaled[1<scaled] = 1

	return @. trunc(Int,(nlevels-1)*scaled) + 1
end