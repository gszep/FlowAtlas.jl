function gate(context::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
        HTTP.Messages.setheader(response,"Content-type" => "text/csv")
        feature = JSON.parse(String(request.body))

        body = IOBuffer(response.body, write=true)
        t = @elapsed save( Stream(format"CSV",body), gate(feature["coordinates"]) )

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

gate(polygon::Vector{<:Any}) = gate( convert( Vector{SVector{2,Float64}}, polygon ) )
function gate(polygon::Vector{SVector{2,Float64}})

    inside(x::SVector{2,Float64}) = inpolygon(x,polygon;in=true,on=false,out=false)
    counts = combine(
    
        ############## for selected events in embedded polygon
        view( data, selections .& map( inside, embeddingCoordinates ), : ),
    
        ############## histogram per channel
        map( channel -> channel => ( x->fit(Histogram,x,channelRange).weights ) => channel, names(data) ))

    counts.channelValue = ( channelRange[1:end-1] + channelRange[2:end] ) / 2
    return counts
end