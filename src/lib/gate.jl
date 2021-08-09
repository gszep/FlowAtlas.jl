function gate(context::NamedTuple,data::DataFrame,embedding::NamedTuple,selections::NamedTuple;channelRange=range(-3,7,length=50))
    request = context.request
	try 
		response = HTTP.Response(200)
        HTTP.Messages.setheader(response,"Content-type" => "text/csv")
        feature = JSON.parse(String(request.body))

        body = IOBuffer(response.body, write=true)
        t = @elapsed save( Stream(format"CSV",body), gate(feature["coordinates"],data,embedding,selections;channelRange=channelRange) )

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

gate(polygon::Vector{<:Any},data::DataFrame,embedding::NamedTuple,selections::NamedTuple;kwargs...) = gate(convert(Vector{SVector{2,Float64}},polygon),data::DataFrame,embedding::NamedTuple,selections::NamedTuple;kwargs...)
function gate(polygon::Vector{SVector{2,Float64}},data::DataFrame,embedding::NamedTuple,selections::NamedTuple;channelRange=range(-3,7,length=50))

    inside(x::SVector{2,Float64}) = inpolygon(x,polygon;in=true,on=false,out=false)
    counts = combine(
    
        ############## for selected events in embedded polygon
        view( data, selections.rows .& map( inside, embedding.coordinates ), : ),
    
        ############## histogram per channel
        map( channel -> channel => ( x->fit(Histogram,filter(xi->~ismissing(xi),x),channelRange).weights ) => channel, Base.names(data) ))

    counts.channelValue = ( channelRange[1:end-1] + channelRange[2:end] ) / 2
    return counts
end