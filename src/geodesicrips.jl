# TODO: remove/replace with FiltrationWrapper
struct EqRips{I, T, R<:AbstractRipsFiltration{I, T}} <: AbstractRipsFiltration{I, T}
    filtration::R
end

adjacency_matrix(e::EqRips) = adjacency_matrix(e.filtration)
threshold(e::EqRips) = threshold(e.filtration)
simplex_type(::Type{<:EqRips{I, T}}, D) where {I, T} = ThickSimplex{D, T, I}

emergent_pairs(::EqRips) = false

struct GeodesicRips{I, T, G, A<:AbstractMatrix{T}} <: AbstractRipsFiltration{I, T}
    graph::G
    dist::A
    paths::Matrix{Int}
    threshold::T
end

function GeodesicRips{I}(graph, distmx=weights(graph); threshold=nothing) where I
    fw = floyd_warshall_shortest_paths(graph, distmx)
    if isnothing(threshold)
        threshold = radius(fw.dists)
    end
    T = eltype(fw.dists)
    return GeodesicRips{I, T, typeof(graph), typeof(fw.dists)}(
        graph, fw.dists, fw.parents, T(threshold)
    )
end
GeodesicRips(args...; kwargs...) = GeodesicRips{Int}(args...; kwargs...)

adjacency_matrix(gr::GeodesicRips) = gr.dist
threshold(gr::GeodesicRips) = gr.threshold
simplex_type(::Type{<:GeodesicRips{I, T}}, D) where {I, T} = ThickSimplex{D, T, I}

emergent_pairs(::GeodesicRips) = false

#=
function _edges(u, vs, adj, r)
    res = map(vs) do v
        adj[u, v]
    end
    if maximum(res) ≤ r && minimum(res) > 0
        return res
    else
        return nothing
    end
end

function _visit!(stack, vs, visited)
    for v in vs
        if !visited[v]
            push!(stack, v)
            visited[v] = true
        end
    end
end

function _min_cofacet(grips::GeodesicRips, σ::ThickSimplex)
    diam = birth(σ)
    vxs = vertices(σ)
    starting_vertex = vxs[1]
    τ = nothing

    # DFS from one of the vertices, only taking vertices that are less than birth(σ)
    # away from all the others.
    visited = falses(nv(grips))
    stack = Int[]
    _visit!(stack, neighbors(grips.graph, starting_vertex), visited)

    while !isempty(stack)
        v = pop!(stack)
        if v in vxs
            _visit!(stack, neighbors(grips.graph, v), visited)
        else
            weights = _edges(v, vxs, grips.dist, diam)
            if !isnothing(weights)
                _visit!(stack, neighbors(grips.graph, v), visited)
                #new_vxs = TupleTools.sort(TupleTools.insertafter(Tuple(vxs), 0, (v,)), rev=true)
                new_vxs, sign = _signed_insert(Tuple(vxs), v)
                # Ignore the sign from now, we will get it from the call to boundary.
                candidate = unsafe_cofacet(grips, σ, new_vxs, v, sign, weights)
                if isnothing(τ) || !isnothing(candidate) && τ > candidate
                    τ = candidate
                end
            end
        end
    end
    return τ
end

function find_apparent_pairs(grips::GeodesicRips{<:Any, T}, cols, progress) where T
    # TODO only works fine for dense grips
    if true || issparse(adjacency_matrix(grips))
        return cols, ()
    end
    S = eltype(cols)
    C = simplex_type(grips, dim(S) + 1)
    cols_left = S[]
    apparent = Tuple{S, C}[]

    if progress
        progbar = Progress(length(cols); desc="Finding apparent pairs... ")
    end
    for σ in cols
        τ = _min_cofacet(grips, σ)
        if isnothing(τ)
            push!(cols_left, σ)
        else
            if σ == maximum(boundary(grips, τ))
                push!(apparent, (σ, τ))
            else
                push!(cols_left, σ)
            end
        end
        progress && next!(progbar)
    end
    progress && printstyled(stderr, "$(length(apparent)) apparent pairs.\n", color=:green)
    return cols_left, apparent
end
=#
function _signed_insert(vertices, vertex)
    n = length(vertices)
    sign = iseven(n) ? 1 : -1
    for i in 0:n-1
        if vertices[i+1] < vertex
            return TupleTools.insertafter(vertices, i, (vertex,)), sign
        end
        sign *= -1
    end
    return TupleTools.insertafter(vertices, n, (vertex,)), sign
end

function _min_cofacet(grips::GeodesicRips, sx::ThickSimplex{1})
    vxs = vertices(sx)
    u, v = vxs
    if has_edge(grips.graph, u, v)
    else
        w = typemax(eltype(sx))
        while true
            u′ = grips.paths[v, u]
            w = min(w, u′)
            u′ == v && break
            u = u′
        end
        new_vxs, sign = _signed_insert(Tuple(vxs), w)
        if w == v || w == typemax(eltype(sx))
            return nothing
        else
            return nothing
            return unsafe_cofacet(grips, sx, new_vxs, v, sign)
        end
    end
end


function find_apparent_pairs(
    grips::GeodesicRips, columns::AbstractVector{<:ThickSimplex{1}}, progress
)
    S = eltype(columns)
    C = simplex_type(grips, dim(S) + 1)
    cols_left = S[]
    apparent = Tuple{S, C}[]

    if progress
        progbar = Progress(length(columns); desc="Finding apparent pairs... ")
    end
    for column in columns
        pivot = _min_cofacet(grips, column)
        if !isnothing(pivot)
            if column == maximum(boundary(grips, pivot))
                push!(apparent, (column, pivot))
            else
                @show column pivot
                error()
                push!(cols_left, column)
            end
        else
            push!(cols_left, column)
        end
        progress && next!(progbar)
    end
    prog_println(progress, "$(length(apparent)) apparent pairs.")
    return cols_left, apparent
end

function filling(grips, element::AbstractChainElement)
    return filling(grips, simplex(element))
end

function filling(grips, sx::AbstractSimplex)
    if dim(sx) > 2
        throw(ArgumentError("currently only dims up to 2 are supported"))
    end
    vxs = vertices(sx)
    result = simplex_type(grips, 1)[]
    for (u, v) in IterTools.subsets(vxs, Val(2))
        while u ≠ v
            u′ = grips.paths[v, u]
            push!(result, simplex(grips, Val(1), (u, u′)))
            u = u′
        end
    end
    return result
end
