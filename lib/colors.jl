function colors!(colors::AbstractMatrix,context::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
        label = JSON.parse(String(request.body))
        t = @elapsed colors!(colors,label)

		@info """$(request.method) $(replace(request.target, r"&.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"&.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function colors!(colors::AbstractMatrix,label::Dict)
	selections = label["name"] ∈ names(labels) ? labels[!,label["name"]] : groups[!,label["name"]]

	colors[:,selections] .= channelview([parse(RGBA,label["color"])])
	return colors
end

function colors!(colors::AbstractMatrix)
    
    colors .= segments #fluorescence[] ? palette[ :, colorIndex[!,channel[]] ] : 

	selection!(selections)
    colors[ :, .~selections ] .= 0.0

    return colors
end

function selection!(selections::Vector{<:Bool})
    selectionCodes = codes .& encode(selection.vals)

	populations = encode(map( name -> name ∈ populationNames, selection.keys ))
	conditions = encode(map( name -> name ∈ conditionNames, selection.keys ))
	groups = encode(map( name -> name ∈ groupNames, selection.keys ))

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

function toIndex(x::AbstractVector{<:Number},channelRange::AbstractVector{<:Number};nlevels::Int=10)

	min,max = extrema(channelRange)
	scaled = ( x .- min ) / ( max - min )
	
	@. scaled[scaled<0] = 0
	@. scaled[1<scaled] = 1

	return @. trunc(Int,(nlevels-1)*scaled) + 1
end

