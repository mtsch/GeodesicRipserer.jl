using LightGraphs
using Ripserer
using GeodesicRipserer
using Test

using Ripserer: adjacency_matrix

function torus_graph(n, m)
    adj = zeros(Int, n*m, n*m)
    cartesians = CartesianIndices((n, m))
    linears = LinearIndices((n, m))
    for vertex_cart in cartesians
        for δi in -1:1, δj in -1:1
            δi == 0 && δj == 0 && continue
            neighbor_cart = vertex_cart + CartesianIndex(δi, δj)
            neighbor_lin = linears[mod1(neighbor_cart[1], n), mod1(neighbor_cart[2], m)]
            vertex_lin = linears[vertex_cart]
            adj[vertex_lin, neighbor_lin] = 1
        end
    end
    return SimpleGraph(adj)
end

function projective_graph(n, m)
    adj = zeros(Int, n*m, n*m)
    cartesians = CartesianIndices((n, m))
    linears = LinearIndices((n, m))
    for vertex_cart in cartesians
        for δi in -1:1, δj in -1:1
            δi == 0 && δj == 0 && continue
            ni = (vertex_cart + CartesianIndex(δi, δj))[1]
            nj = (vertex_cart + CartesianIndex(δi, δj))[2]
            if 1 ≤ ni ≤ n
                ni = mod1(ni, n)
                nj = mod1(nj, m)
            else
                ni = mod1(ni, n)
                nj = mod1(-nj, m)
            end
            neighbor_lin = linears[ni, nj]
            vertex_lin = linears[vertex_cart]
            adj[vertex_lin, neighbor_lin] = 1
        end
    end
    return SimpleGraph(adj)
end

@testset "GeodesicRips" begin
    @testset "4-vertex cycle - make sure small holes don't get skipped" begin
        cycle = cycle_graph(4)
        @test ripserer(GeodesicRips(cycle))[2] == [(1, 2)]
    end

    @testset "Unlike Rips, finds equilateral triangle in 6-vertex cycle" begin
        cycle = cycle_graph(6)
        grips = GeodesicRips(cycle)
        vs = vertices(ripserer(grips)[2][1].death_simplex)
        @test vs[1] - vs[2] == 2
        @test vs[2] - vs[3] == 2
        @test vs[3] - vs[1] + 6  == 2
    end

    @testset "Ex-bug" begin
        # This used to fail.
        # Example was generated randomly. Full precision is needed to make it fail.
        cycle = cycle_graph(6)
        weights = [
            0.0     0.842   0.0     0.0     0.0     0.8096
            0.842   0.0     0.7734  0.0     0.0     0.0
            0.0     0.7734  0.0     0.8613  0.0     0.0
            0.0     0.0     0.8613  0.0     0.169   0.0
            0.0     0.0     0.0     0.169   0.0     0.5273
            0.8096  0.0     0.0     0.0     0.5273  0.0
        ]
        grips = GeodesicRips(cycle, Float16.(weights))
        @test ripserer(grips) == ripserer(adjacency_matrix(grips))
    end

    for m in (2, 17)
        @testset "18 cycle graph modulus=$m" begin
            grips = GeodesicRips(cycle_graph(18))
            @test ripserer(grips; dim_max=2, modulus=2) ==
                ripserer(adjacency_matrix(grips); dim_max=2, modulus=m)
        end
    end

    for m in (2, 5)
        @testset "12×9 torus graph modulus=$m" begin
            torus = torus_graph(12, 9)
            grips = GeodesicRips(torus)
            result = ripserer(grips, modulus=m)

            @test result == ripserer(adjacency_matrix(grips), modulus=m)
            @test result[1] == vcat(fill((0, 1), 107), [(0, Inf)])
            @test result[2][1] == (1, 3)
            @test result[2][2] == (1, 4)
            @test circumference(grips, result[2][1].death_simplex) == 9
            @test circumference(grips, result[2][2].death_simplex) == 12
        end
    end

    @testset "5×7 projective plane graph modulus=2" begin
        plane = projective_graph(5, 7)
        grips = GeodesicRips(plane)
        result = ripserer(GeodesicRips(plane), dim_max=3)

        @test result == ripserer(adjacency_matrix(grips), dim_max=3)
        @test length(result[2]) == 2
        @test length(result[3]) == 1
        @test circumference(grips, result[2][1].death_simplex) == 5
        @test circumference(grips, result[2][2].death_simplex) == 4
    end

    @testset "5×7 projective plane graph modulus=3" begin
        plane = projective_graph(5, 7)
        grips = GeodesicRips(plane)
        result = ripserer(grips, modulus=3, dim_max=3)

        @test result == ripserer(adjacency_matrix(grips), modulus=3, dim_max=3)
        @test length(result[2]) == 1
        @test circumference(grips, result[2][1].death_simplex) == 4
    end
end
