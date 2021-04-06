function count(context::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
        HTTP.Messages.setheader(response,"Content-type" => "application/json")

        body = IOBuffer(response.body, write=true)
        t = @elapsed write( body, JSON.json(count()) )

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function count()
	return [ (names = decode(code, [names(labels);names(groups)]), count = count) for (code, count) ∈ countmap(codes[selections]) ]
end