function selection!(context::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
        label = JSON.parse(String(request.body))
        t = @elapsed selection!(label)

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function selection!(label::Dict)
	global populationNames,conditionNames,groupNames

    if label["name"] ∈ selection.keys
	    selection[label["name"]] = label["selected"]

    elseif label["name"] == "populations"
        for name ∈ populationNames selection[name]=label["selected"] end

    elseif label["name"] == "conditions"
        for name ∈ conditionNames selection[name]=label["selected"] end

    elseif label["name"] == "groups"
        for name ∈ groupNames selection[name]=label["selected"] end

    elseif typeof(label["name"]) <: Dict
		populationNames = label["name"]["populations"]
		conditionNames = label["name"]["conditions"]
		groupNames = label["name"]["groups"]
    end

	populations = encode(map( name -> name ∈ populationNames, selection.keys ))
	conditions = encode(map( name -> name ∈ conditionNames, selection.keys ))
	groups = encode(map( name -> name ∈ groupNames, selection.keys ))

    selectionCodes = codes .& encode(selection.vals)
    @. selections = ( selectionCodes & populations ≠ 0 ) & ( selectionCodes & conditions ≠ 0 ) & ( selectionCodes & groups ≠ 0 )
	return selections
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

function decode( code::Int, codeMap::Vector{<:Any} )
	return codeMap[[ bit=='1' for bit ∈ bitstring(code)[end-length(codeMap)+1:end] ]]
end