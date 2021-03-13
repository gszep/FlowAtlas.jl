struct ClosedVector{T} <: AbstractVector{T}
    x::AbstractVector{T}

    function ClosedVector(x::AbstractVector{T}) where T
        return new{T}(x)
    end
end

Base.size(A::ClosedVector) = (size(A.x,1)+1,)
Base.length(A::ClosedVector) = length(A.x)+1
Base.IndexStyle(::Type{ClosedVector}) = IndexCartesian() 

Base.getindex(A::ClosedVector, i::Int) = getindex(A.x, i == length(A) ? 1 : i)
Base.setindex!(A::ClosedVector, v, i::Int) = setindex!(A.x, v, i == length(A) ? 1 : i) 

struct Gate
    polygon::Observable{<:ClosedVector}
    selected::Observable{<:Bool}
    color::Observable{<:String}
end

Gate(x::AbstractVector{<:StaticVector}) = Gate( (Observable ∘ ClosedVector)(x), Observable(false), Observable("#EEEEEE") )
Gate(x::AbstractVector{<:StaticVector},select::Bool) = Gate( (Observable ∘ ClosedVector)(x), Observable(select), Observable("#EEEEEE") )
Gate(x::AbstractVector{<:StaticVector},color::String) = Gate( (Observable ∘ ClosedVector)(x), Observable(false), Observable(color) )

Gate(x::AbstractVector{<:StaticVector},select::Bool,color::String) = Gate( (Observable ∘ ClosedVector)(x), Observable(select), Observable(color) )
Gate(x::AbstractVector{<:StaticVector},color::String,select::Bool) = Gate( (Observable ∘ ClosedVector)(x), Observable(select), Observable(color) )

struct Group
    name::String
    selected::Observable{<:Bool}
    color::Observable{<:String}
end

Group(name::String) = Group(name,Observable(true),Observable("#EEEEEE"))
Group(name::String,select::Bool) = Group(name,Observable(select),Observable("#EEEEEE"))
Group(name::String,color::String) = Group(name,Observable(true),Observable(color))

Group(name::String,select::Bool,color::String) = Group(name,Observable(select),Observable(color))
Group(name::String,color::String,select::Bool) = Group(name,Observable(select),Observable(color))