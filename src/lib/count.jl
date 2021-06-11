function count(context::NamedTuple,selections::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
        HTTP.Messages.setheader(response,"Content-type" => "application/json")

        body = IOBuffer(response.body, write=true)
        t = @elapsed write( body, JSON.json(count(selections)) )

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function count(selections::NamedTuple)
	return [ (names = decode(code, selections.names.keys ), count = count) for (code, count) ∈ countmap(selections.codes[selections.rows]) ]
end