module GeodesicRipserer

using Base: @propagate_inbounds
using IterTools
using LightGraphs
using Ripserer
using TupleTools

# Filtration overloads
import Ripserer:
    adjacency_matrix, threshold, simplex_type, emergent_pairs, find_apparent_pairs
# Simplex overloads
import Ripserer:
    birth, index, unsafe_cofacet, unsafe_simplex

using Ripserer:
    AbstractRipsFiltration, AbstractSimplex, AbstractChainElement,
    births, prog_print, prog_println, radius

export GeodesicRips, ThickSimplex, thickness

include("thicksimplex.jl")
include("geodesicrips.jl")

end
