function colors!(context::NamedTuple,labels::DataFrame,groups::DataFrame,colors::NamedTuple)
    request = context.request
	try 
		response = HTTP.Response(200)
        label = JSON.parse(String(request.body))
        t = @elapsed colors!(label,labels,groups,colors)

		@info """$(request.method) $(replace(request.target, r"\?.+" => "")) | $t seconds"""
		return response

	catch exception
		printstyled("ERROR $(replace(request.target, r"\?.+" => "")) | $exception\n",color=:red)
		return HTTP.Response(500)
	end
end

function colors!(label::Dict,labels::DataFrame,groups::DataFrame,colors::NamedTuple)
	colors.labels.rows[1:3, label["name"] ∈ Base.names(labels) ? labels[!,label["name"]] : groups[!,label["name"]]] .= channelview([parse(RGB,label["color"])])
	return colors
end

function toIndex(x::AbstractVector{<:Union{Number,Missing}};nlevels::Int=10,p::Real=0.1)

	min,max = quantile(skipmissing(x),(p,1-p))
	scaled = ( x .- min ) / ( max - min )
	
	@. scaled[(scaled<0)&(~ismissing(scaled))] = 0
	@. scaled[(1<scaled)&(~ismissing(scaled))] = 1

	scaled = @. trunc(Union{Int,Missing},(nlevels-1)*scaled) + 1
	return replace!(scaled,missing=>nlevels+1)
end