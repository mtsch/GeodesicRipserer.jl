"""
    ThickSimplex{D, T, I<:Integer} <: AbstractSimplex{D, T, I}

Like `Simplex`, but it also keeps track of its shortest edge.
"""
struct ThickSimplex{D, T, I<:Integer} <: AbstractSimplex{D, T, I}
    index::I
    longest::T
    shortest::T

    function ThickSimplex{D, T, I}(index::Integer, longest, shortest) where {D, T, I}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        return new{D, T, I}(I(index), T(longest), T(shortest))
    end
end

function ThickSimplex{D}(vertices, longest::T, shortest::T) where {D, T}
    idx = index(TupleTools.sort(Tuple(vertices), rev=true))
    return ThickSimplex{D, T, eltype(vertices)}(idx, longest, shortest)
end
function ThickSimplex{1}(index_or_vertices, birth)
    return ThickSimplex{1}(index_or_vertices, birth, birth)
end

birth(sx::ThickSimplex) = sx.longest
index(sx::ThickSimplex) = abs(sx.index)
Base.sign(sx::ThickSimplex) = sign(sx.index)
function Base.:-(sx::ThickSimplex{D, T, I}) where {D, T, I}
    return ThickSimplex{D, T, I}(-sx.index, sx.longest, sx.shortest)
end

"""
    thickness(sx::ThickSimplex)

"""
thickness(sx::ThickSimplex) = sx.shortest

function Base.isless(sx1::S, sx2::S) where S<:ThickSimplex
    return (birth(sx1), thickness(sx1), -index(sx1)) < (birth(sx2), thickness(sx2), -index(sx2))
end

function Base.show(io::IO, sx::ThickSimplex{D}) where D
    print(io, sign(sx) == 1 ? :+ : :-, nameof(typeof(sx)), "{", D, "}(",
          vertices(sx), ", ", birth(sx), ", ", thickness(sx), ")")
end

# Dense version.
@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration,
    simplex::ThickSimplex,
    cofacet_vertices,
    new_vertex,
    sign,
) where {I, T, D, S<:ThickSimplex{D, T, I}}
    diameter = birth(simplex)
    shortest = simplex.shortest
    adj = adjacency_matrix(rips)
    for v in cofacet_vertices
        v == new_vertex && continue
        e = adj[new_vertex, v]
        (iszero(e) || e > threshold(rips)) && return nothing
        diameter = ifelse(e > diameter, e, diameter)
        shortest = ifelse(e < shortest, e, shortest)
    end
    return S(I(sign) * index(cofacet_vertices), diameter, shortest)
end

# Sparse version.
@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration,
    simplex::ThickSimplex,
    cofacet_vertices,
    ::Any,
    sign,
    new_edges,
) where {I, T, D, S<:ThickSimplex{D, T, I}}
    diameter = birth(simplex)
    shortest = simplex.shortest
    for e in new_edges
        e > threshold(rips) && return nothing
        diameter = ifelse(e > diameter, e, diameter)
        shortest = ifelse(e < shortest, e, shortest)
    end
    return S(I(sign) * index(cofacet_vertices), diameter, shortest)
end

function unsafe_simplex(
    ::Type{S}, rips::AbstractRipsFiltration{I, T}, vertices, sign
) where {I, T, S<:ThickSimplex{<:Any, T, I}}
    if dim(S) == 0
        return S(I(sign) * vertices[1], births(rips)[vertices[1]], zero(T))
    else
        adj = adjacency_matrix(rips)
        shortest = typemax(T)
        diameter = typemin(T)
        n = length(vertices)
        for i in 1:n, j in i+1:n
            e = adj[vertices[i], vertices[j]]
            (iszero(e) || e > threshold(rips)) && return nothing
            diameter = ifelse(e > diameter, e, diameter)
            shortest = ifelse(e < shortest, e, shortest)
        end
        return S(I(sign) * index(vertices), diameter, shortest)
    end
end
