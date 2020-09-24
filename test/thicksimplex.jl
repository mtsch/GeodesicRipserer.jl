using LightGraphs
using Ripserer
using GeodesicRipserer
using Test

@testset "ThickSimplex" begin
    @testset "Constructors" begin
        @test ThickSimplex{1}([2, 1], 1) == ThickSimplex{1}(1, 1, 2)
        @test ThickSimplex{1}((1, 3), 0f0) == ThickSimplex{1}(2, 0f0, 1f0)
        @test ThickSimplex{2}([5, 2, 1], 1.0, 2.0) == ThickSimplex{2}(5, 1.0, 2.0)
    end

    @testset "AbstractSimplex interface" begin
        for D in (0, 2),
            T in (Float32, Int),
            I in (Int32, Int64)
            @testset "ThickSimplex{$D, $T, $I}" begin
                b = T(9)
                t = T(18)
                i = I(100)
                sx = ThickSimplex{D}(i, b, t)

                @testset "Type stuff" begin
                    @test sx ≡ ThickSimplex{D, T, I}(i, b, t)
                    @test typeof(sx) ≡ ThickSimplex{D, T, I}
                    D > 0 && @test_throws DomainError ThickSimplex{-D}(i, b, t)
                    @test eltype(sx) == I
                end

                @testset "Getters" begin
                    @test index(sx) == i
                    @test index(-sx) == i
                    @test birth(sx) ≡ b
                    @test thickness(sx) ≡ t
                end

                @testset "Equality, hashing" begin
                    @test sx == sx
                    @test sx != ThickSimplex{D + 1}(i, b, t)
                    @test sx == ThickSimplex{D}(i, b + 1, t)
                    @test isequal(sx, ThickSimplex{D}(i, b + 1, t))
                    @test hash(sx) == hash(ThickSimplex{D}(i, b, t + 1))
                    @test hash(sx) == hash(index(sx))
                end

                @testset "Signs" begin
                    @test +sx === sx
                    @test sx == -sx
                    @test sx ≡ -(-sx)
                    @test sign(sx) == 1
                    @test sign(-sx) == -1
                    @test abs(sx) ≡ sx
                    @test abs(-sx) ≡ sx
                    @test dim(sx) == D
                end

                @testset "Ordering" begin
                    @test sx < ThickSimplex{D}(I(i + 1), b + 1, t)
                    @test sx > ThickSimplex{D}(I(i - 1), b - 1, t)

                    @test sx < ThickSimplex{D}(I(i - 1), b, t)
                    @test sx > ThickSimplex{D}(I(i + 1), b, t)

                    @test sx < ThickSimplex{D}(i, b, t + 1)
                    @test sx > ThickSimplex{D}(i, b, t - 1)
                end

                @testset "Array interface, vertices" begin
                    verts = vertices(sx)

                    @test eltype(sx) == eltype(verts)
                    @test length(sx) == length(verts) == D + 1
                    @test size(sx) == (D + 1,)
                    @test firstindex(sx) == 1
                    @test lastindex(sx) == D + 1

                    @test ThickSimplex{D}(verts, birth(sx), thickness(sx)) ≡ sx

                    for (i, v) in enumerate(sx)
                        @test v == verts[i]
                    end

                    @test begin @inferred vertices(sx); true end
                end

                @testset "Printing" begin
                    @test sprint(show, sx) ==
                        "+ThickSimplex{$D}($(vertices(sx)), $(b), $(t))"
                    @test sprint(show, -sx) ==
                        "-ThickSimplex{$D}($(vertices(sx)), $(b), $(t))"
                    @test sprint((i, s) -> show(i, MIME"text/plain"(), s), sx) ==
                        "$D-dimensional ThickSimplex(index=$i, birth=$b):\n  +$(vertices(sx))"
                end
            end
        end
    end
end
