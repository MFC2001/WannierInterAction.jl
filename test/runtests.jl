using Test
using WannierInterAction

@testset "WannierInterAction.jl" begin
    lattice = Lattice([1 0 0; 0 1 0; 0 0 1])
    Glist = [
        ReducedCoordinates{Int}([0, 0, 0]),
        ReducedCoordinates{Int}([1, 0, 0]),
        ReducedCoordinates{Int}([0, 1, 0]),
        ReducedCoordinates{Int}([0, 0, 1])
    ]
    Gidx = [[1, 2, 3, 4]]
    c1 = [0.5+0.5im, -0.5+0.5im, 0.0+0.0im, 0.0+0.0im] # norm = 1
    c2 = [0.5+0.0im, 0.5+0.0im, 0.5+0.0im, 0.5+0.0im]  # norm = 1
    components = [hcat(c1, c2)]

    kgrid = RedKgrid(MonkhorstPack([1, 1, 1]))
    wave = WannierPWKgrid(kgrid, components, Gidx, Glist, lattice;
        nspinor=1, fft_grid=[4, 4, 4], period=[1, 1, 0])

    pair = [
        (1, 1, [0, 0, 0]),
        (1, 2, [0, 0, 0]),
        (1, 2, [1, 0, 0]),
    ]
    ρk = Rhok(wave)
    mktempdir() do tmp_dir
        int(ρk; ρRbuffer=:low, pair, outfolder=tmp_dir)
        @test isfile(joinpath(tmp_dir, "U_bare.jld2"))
        @test isfile(joinpath(tmp_dir, "pair_bare.jld2"))
        int(ρk; ρRbuffer=:high, pair, outfolder=tmp_dir)
        int(ρk; ρRbuffer=:medium, pair, outfolder=tmp_dir)
    end

    wave = WannierPWKgrid(kgrid, components, Gidx, Glist, lattice;
        nspinor=1, fft_grid=[4, 4, 4], period=[1, 1, 1])
    ρk = Rhok(wave)
    mktempdir() do tmp_dir
        int(ρk; ρRbuffer=:low, pair, outfolder=tmp_dir)
        @test isfile(joinpath(tmp_dir, "U_bare.jld2"))
        @test isfile(joinpath(tmp_dir, "pair_bare.jld2"))
    end
end
