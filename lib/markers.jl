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

