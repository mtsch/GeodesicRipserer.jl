using SafeTestsets

@safetestset "Aqua.jl" begin
    include("aqua.jl")
end
@safetestset "thicksimplex" begin
    include("thicksimplex.jl")
end
@safetestset "geodesicrips" begin
    include("geodesicrips.jl")
end
