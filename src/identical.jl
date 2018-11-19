#import Base.push!
module identical


export AbstractIdentical, Identical, areidentical


abstract type AbstractIdentical end

mutable struct Identical <: AbstractIdentical
    identical::Vector{Set}
    Identical() = new([])
    Identical(n::Integer) = new(Set.(1:n))
end


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

function areidentical(a::AbstractIdentical, x)
    s = Set(x)
    for x in a.identical
        if issubset(s,x)
            return true
        end
    end
    return false
end



end #module
