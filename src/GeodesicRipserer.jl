module GeodesicRipserer

using Base: @propagate_inbounds
using Compat
using IterTools
using LightGraphs
using NearestNeighbors
using PersistenceDiagrams
using ProgressMeter
using Random
using Ripserer
using SparseArrays
using StaticArrays
using TupleTools

# Filtration overloads
import Ripserer:
    adjacency_matrix, threshold, simplex_type, emergent_pairs, find_apparent_pairs,
    postprocess_diagram
# Simplex overloads
import Ripserer:
    birth, index, unsafe_cofacet, unsafe_simplex

using Ripserer:
    AbstractRipsFiltration, AbstractSimplex, AbstractChainElement,
    boundary, births, fmt_number, prog_print, prog_println, radius

export GeodesicRips, ThickSimplex, circumference, filling, thickness

include("thicksimplex.jl")
include("geodesicrips.jl")

end
