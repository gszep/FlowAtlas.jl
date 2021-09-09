function tile(context::NamedTuple, embedding::NamedTuple, selections::NamedTuple, colors::NamedTuple, names::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
		HTTP.Messages.setheader(response,"Content-type" => "image/png")
		channel = match(r"(?<=channel=)[^&]*(?=&)",request.target)
		scale = match(r"(?<=scale=)[^&]*(?=&)",request.target)

		body = IOBuffer(response.body, write=true)
		t = @elapsed save( Stream(format"PNG",body), tile(request.target,embedding,selections,colors,names; channel=isnothing(channel) ? "" : channel.match, scale=isnothing(scale) ? 1.0 : parse(Float32,scale.match)) )

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function tile(url::AbstractString, embedding::NamedTuple, selections::NamedTuple, colors::NamedTuple, names::NamedTuple; channel::AbstractString="", scale::Number=1.0)
    zxy = broadcast(parse, Int, match(r"(\d+)/(\d+)/(\d+)", url ).captures )
    return colorview(RGBA,tile(zxy...,embedding,selections,colors,names; channel=channel, scale=scale))
end

function tile( z::Int, x::Int, y::Int, embedding::NamedTuple, selections::NamedTuple, colors::NamedTuple, names::NamedTuple; channel::AbstractString="", padIndex::Int=16, scale::Number=1.0)
	(xmin,xmax), (ymin,ymax) = embedding.extrema
	imageSize = max(xmax-xmin,ymax-ymin)

    tileSize = imageSize / 2^z
    padding = tileSize/padIndex

    x = xmin + x*tileSize
    y = ymin + y*tileSize

    return rasterKernelCircle( scale*sqrt(z),
        rasterize( (256+2padIndex,256+2padIndex), embedding.array[:,selections.rows],
			channel ∈ names.channels ? colors.channels.colors[:,colors.channels.rows[selections.rows,channel]] : colors.labels.rows[:,selections.rows],
		
            xlim=( x-padding, x+tileSize+padding),
            ylim=( y-padding, y+tileSize+padding)
		)
	)[:,padIndex+1:end-padIndex,padIndex+1:end-padIndex]
end