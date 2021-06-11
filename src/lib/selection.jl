function selection!(context::NamedTuple,selections::NamedTuple,names::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
        label = JSON.parse(String(request.body))
        t = @elapsed selection!(label,selections,names)

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function selection!(label::Dict,selections::NamedTuple,names::NamedTuple)

    if label["name"] ∈ selections.names.keys
	    selections.names[label["name"]] = label["selected"]

    elseif label["name"] == "populations"
        for name ∈ names.populations selections.names[name]=label["selected"] end

    elseif label["name"] == "conditions"
        for name ∈ names.conditions selections.names[name]=label["selected"] end

    elseif label["name"] == "groups"
        for name ∈ names.groups selections.names[name]=label["selected"] end

    elseif typeof(label["name"]) <: Dict
		names.populations = label["name"]["populations"]
		names.conditions = label["name"]["conditions"]
		names.groups = label["name"]["groups"]
    end

	populations = encode(map( name -> name ∈ names.populations, selections.names.keys ))
	conditions = encode(map( name -> name ∈ names.conditions, selections.names.keys ))
	groups = encode(map( name -> name ∈ names.groups, selections.names.keys ))

    selectionCodes = selections.codes .& encode(selections.names.vals)
    @. selections.rows = ( selectionCodes & populations ≠ 0 ) & ( selectionCodes & conditions ≠ 0 ) & ( selectionCodes & groups ≠ 0 )
	return selections.rows
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