function tile(context::NamedTuple; extrema::Array{<:Tuple}=[(-2,2) (-2,2)])
    request = context.request
	try 
		response = HTTP.Response(200)
		HTTP.Messages.setheader(response,"Content-type" => "image/png")
		channel = match(r"(?<=channel=).*(?=&)",request.target)

		body = IOBuffer(response.body, write=true)
		t = @elapsed save( Stream(format"PNG",body), tile(request.target; extrema=extrema, channel=isnothing(channel) ? "" : channel.match) )

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function tile(url::AbstractString; extrema::Array{<:Tuple}=[(-2,2) (-2,2)], channel::AbstractString="")
    zxy = broadcast(parse, Int, match(r"(\d+)/(\d+)/(\d+)", url ).captures )
    return colorview(RGBA,tile(zxy...;extrema=extrema,channel=channel))
end

function tile( z::Int, x::Int, y::Int; extrema::Array{<:Tuple}=[(-2,2) (-2,2)], channel::AbstractString="", padIndex::Int=16)
	(xmin,xmax), (ymin,ymax) = extrema
	imageSize = max(xmax-xmin,ymax-ymin)

    tileSize = imageSize / 2^z
    padding = tileSize/padIndex

    x = xmin + x*tileSize
    y = ymin + y*tileSize

    return rasterKernelCircle( sqrt(z),
        rasterize( (256+2padIndex,256+2padIndex), embedding[:,selections],
			channel ∈ names(data) ? channelPalette[:,colorIndex[selections,channel]] : labelColors[:,selections],
		
            xlim=( x-padding, x+tileSize+padding),
            ylim=( y-padding, y+tileSize+padding)
		)
	)[:,padIndex+1:end-padIndex,padIndex+1:end-padIndex]
end