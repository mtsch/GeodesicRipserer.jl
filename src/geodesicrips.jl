# TODO: remove/replace with FiltrationWrapper
struct EqRips{I, T, R<:AbstractRipsFiltration{I, T}} <: AbstractRipsFiltration{I, T}
    filtration::R
end

adjacency_matrix(e::EqRips) = adjacency_matrix(e.filtration)
threshold(e::EqRips) = threshold(e.filtration)
simplex_type(::Type{<:EqRips{I, T}}, D) where {I, T} = ThickSimplex{D, T, I}

emergent_pairs(::EqRips) = false

struct GeodesicRips{I, T, G, A<:AbstractMatrix{T}, P} <: AbstractRipsFiltration{I, T}
    graph::G
    dist::A
    paths::Matrix{Int}
    points::P
    threshold::T
end

function GeodesicRips{I}(
    graph::AbstractGraph, distmx=weights(graph), points=nothing; threshold=nothing
) where I
    fw = floyd_warshall_shortest_paths(graph, distmx)
    if isnothing(threshold)
        threshold = radius(fw.dists)
        if isinf(threshold)
            error("TODO: space is disconnected.")
        end
    end
    T = eltype(fw.dists)
    return GeodesicRips{I, T, typeof(graph), typeof(fw.dists), typeof(points)}(
        graph, fw.dists, fw.parents, points, T(threshold)
    )
end
GeodesicRips(args...; kwargs...) = GeodesicRips{Int}(args...; kwargs...)

function GeodesicRips{I}(
    points::AbstractVector, r1, r2=nothing; metric=Euclidean(), kwargs...
) where I
    !isnothing(r2) && r1 > r2 && @warn "you probably want r1 ≤ r2"

    if isnothing(r2)
        r2 = r1
    else
        points = shuffle(points)
        tree = KDTree(SVector.(points), metric)
        # Create a subsample of points r1 apart.
        visited = falses(length(points))
        new_points = eltype(points)[]
        for i in eachindex(points)
            visited[i] && continue
            push!(new_points, points[i])
            visited[inrange(tree, SVector(points[i]), r1)] .= true
        end
        points = new_points
    end

    is = Int[]
    js = Int[]
    vs = typeof(metric(SVector(points[1]), SVector(points[2])))[]
    edgelist = Edge{Int}[]
    tree = KDTree(SVector.(points), metric)
    for u in eachindex(points)
        for v in inrange(tree, SVector(points[u]), 2r2)
            u == v && continue
            d = metric(SVector(points[u]), SVector(points[v]))
            append!(is, (u, v))
            append!(js, (v, u))
            append!(vs, (d, d))
            push!(edgelist, Edge(u, v))
        end
    end
    weights = sparse(is, js, vs, length(points), length(points), min)
    graph = Graph(edgelist)

    return GeodesicRips{I}(graph, weights, points; kwargs...)
end

adjacency_matrix(gr::GeodesicRips) = gr.dist
threshold(gr::GeodesicRips) = gr.threshold
simplex_type(::Type{<:GeodesicRips{I, T}}, D) where {I, T} = ThickSimplex{D, T, I}

emergent_pairs(::GeodesicRips) = false

function _edges(u, vxs, adj, r)
    res = map(vxs) do v
        adj[u, v]
    end
    if maximum(res) ≤ r && minimum(res) > 0
        return res
    else
        return nothing
    end
end

function _signed_insert(vertices, vertex)
    n = length(vertices)
    sign = iseven(n) ? 1 : -1
    for i in 0:n-1
        if vertices[i+1] < vertex
            return TupleTools.insertafter(Tuple(vertices), i, (vertex,)), sign
        end
        sign *= -1
    end
    return TupleTools.insertafter(Tuple(vertices), n, (vertex,)), sign
end

function _min_cofacet(grips, sx)
    vxs = vertices(sx)
    adj = adjacency_matrix(grips)
    graph = grips.graph
    result = nothing
    for u in vxs
        for v in neighbors(grips.graph, u)
            es = _edges(v, vxs, adj, birth(sx))
            isnothing(es) && continue
            new_vxs, sign = _signed_insert(vxs, v)
            candidate = unsafe_cofacet(grips, sx, new_vxs, v, sign, es)
            if isnothing(result) || candidate < result
                result = candidate
            end
        end
    end
    if !isnothing(result) && thickness(result) < thickness(sx)
        return result
    else
        return nothing
    end
end

function find_apparent_pairs(grips::GeodesicRips, columns, progress)
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
                push!(cols_left, column)
            end
        else
            push!(cols_left, column)
        end
        progress && next!(
            progbar, showvalues=((:left, length(cols_left)), (:apparent, length(apparent)))
        )
    end
    prog_println(progress, fmt_number(length(apparent)), " apparent pairs.")
    return cols_left, apparent
end

function postprocess_diagram(grips::GeodesicRips, diagram)
    intervals = map(diagram) do interval
        if isnothing(interval.death_simplex)
            death_filling = (;death_filling=simplex_type(grips, 1)[])
        else
            death_filling = (;death_filling=filling(grips, interval.death_simplex))
        end
        birth_filling = (;birth_filling=filling(grips, interval.birth_simplex))

        meta = (;interval.meta..., birth_filling..., death_filling...)
        return PersistenceInterval(interval.birth, interval.death, meta)
    end
    return PersistenceDiagram(intervals, diagram.meta)
end

function filling(grips::GeodesicRips, element::AbstractChainElement)
    return filling(grips, simplex(element))
end

function filling(grips::GeodesicRips, sx::AbstractSimplex)
    if dim(sx) == 0
        return [sx]
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

function circumference(grips::GeodesicRips, sx::AbstractSimplex)
    sum(birth(e) for e in filling(grips, sx))
end
