"""
Methods to help handling of identical atoms in molecules
"""
module identical


export AbstractIdentical,
       areidentical,
       Identical


abstract type AbstractIdentical end


"""
    Identical <: AbstractIdentical

Struct to handle identical information for atoms or any object in general

# Fields
- `identical::Vector{Set}` :
"""
mutable struct Identical <: AbstractIdentical
    identical::Vector{Set}
    Identical() = new([])
    Identical(n::Integer) = new(Set.(1:n))
end

"""
    Base.push!(a::AbstractIdentical, x)

Adds identical information
"""
function Base.push!(a::AbstractIdentical, x)
    s = Set(x)
    if isempty(a.identical)
        push!(a.identical, s)
        return a
    end
    n = Set()
    d = []
    for i in eachindex(a.identical)
        if !isempty( intersect(s,a.identical[i]) )
            if isempty(n)
                n = a.identical[i]
                union!(n,s)
            else
                union!(n,a.identical[i])
                push!(d,i)
            end
        end
    end
    if isempty(d) && isempty(n)
        push!(a.identical,s)
    else
        deleteat!(a.identical, d)
    end
    return a
end

"""
    areidentical(a::AbstractIdentical, x)

Tests are objects x

# Arguments
- 'identical::AbstractIdentical' : object holding identical data
- 'x' : collectable holiding keys of object to be tested for identical information
"""
function areidentical(identical::AbstractIdentical, x)
    s = Set(x)
    for y in identical.identical
        if issubset(s,y)
            return true
        end
    end
    return false
end



end #module
