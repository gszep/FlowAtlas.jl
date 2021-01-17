using GigaSOM,EzXML,DataFrames
using Glob: GlobMatch

using DataFrames: CategoricalValue
import Base: show # display categorical types neatly
show(io::IO, ::MIME"text/html", x::CategoricalValue) = print(io, get(x))

function load(path::String; workspace::String="", cofactor::Number=250, channelMap::Dict=Dict(), kwargs...)

	#######################################
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
		append!(labels,label)

		append!(groups,group,cols=:union)
		gatings[path] = gating
		
	end

	if isempty(workspace)
		return data,nothing,nothing,nothing
	end

	for name ∈ names(groups)
		replace!(groups[!,name],missing=>false)
	end

	disallowmissing!(groups)
	return data,labels, groups,gatings
end