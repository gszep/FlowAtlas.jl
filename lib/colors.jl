function colors!(context::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
        label = JSON.parse(String(request.body))
        t = @elapsed colors!(label)

		@info """$(request.method) $(replace(request.target, r"&.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"&.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function colors!(label::Dict)
	colors[1:3, label["name"] ∈ names(labels) ? labels[!,label["name"]] : groups[!,label["name"]]] .= channelview([parse(RGB,label["color"])])
	return colors
end

function toIndex(x::AbstractVector{<:Number},channelRange::AbstractVector{<:Number};nlevels::Int=10)

	min,max = extrema(channelRange)
	scaled = ( x .- min ) / ( max - min )
	
	@. scaled[scaled<0] = 0
	@. scaled[1<scaled] = 1

	return @. trunc(Int,(nlevels-1)*scaled) + 1
end