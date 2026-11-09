export MirrorCorrection

include("./_getU_out_of_rcut.jl")
include("./Gauss.jl")
include("./KeldyshGauss.jl")

#TODO: 或许需要带着相位因子展开。

"""
	MirrorCorrection(sym::Symbol, U::AbstractHR, kgrid, lattice, orblocat; outfile::Union{Nothing,AbstractString}=nothing,
		logfile::AbstractString=isnothing(outfile) ? "./MirrorCorrection.log" : joinpath(dirname(outfile), "MirrorCorrection.log"),
		kwargs...)

Apply a mirror correction to the direct interaction `U`.

In a periodic calculation, `U` contains the interaction between the Wannier functions
and those of all periodic images (the "mirror" systems). The function subtracts the
periodic sum of a model long-range potential from `U` and adds back the corresponding
real-space potential, so that the corrected interaction describes a single cell in an
infinite system. 

The model potential is selected by `sym`: `:gauss` and `:gauss_optim`
use a Gaussian-damped Coulomb potential for 3D systems, whereas `:keldysh_gauss` and
`:keldysh_gauss_optim` use a Keldysh potential with Gaussian damping along the
non-periodic direction for 2D (slab) systems.

# Arguments
- `sym::Symbol`: The correction method, one of `:gauss`, `:gauss_optim`,
  `:keldysh_gauss` and `:keldysh_gauss_optim`.
- `U::AbstractHR`: The direct interaction between Wannier functions in the periodic cell,
  as computed by [`int`](@ref) (e.g. the `interaction` field in `U_bare.jld2` in the
  examples). It is converted to a full `HRFull` internally and must extend far enough for
  its remote pairs to have reached the asymptotic long-range behavior.
- `kgrid`: The k-point grid of the Wannier calculation, given as a `MonkhorstPack`, a
  `RedKgrid`, a vector of three integers or a 3-tuple of integers. Only the grid size is
  used, and the mesh is rebuilt internally with zero shift, so it must match the size of grid on
  which `U` was computed. For 2D methods its size along the non-periodic direction must
  be one.
- `lattice`: The real-space lattice of the periodic cell, as a `Lattice` or a 3×3 matrix
  whose columns are the lattice vectors. Lengths are in Å and energies in eV.
- `orblocat`: The Wannier-center positions, one per orbital and in the order of the
  orbital indices of `U`. Elements of type `<:CartesianCoordinates` are taken as
  Cartesian coordinates, other vectors as fractional (reduced) coordinates.

# Keyword arguments
- `rcut::Real`: The short-range cutoff in Å, shared by all `_optim` variants and
  required there. The long-range parameters are fitted only to the plane-wave-basis
  interaction values of the pairs beyond `rcut`, so that they are as reasonable as
  possible; every orbital pair must have at least one interaction beyond `rcut`.
- `outfile::Union{Nothing,AbstractString}=nothing`: Path of a JLD2 file in which the
  returned `NamedTuple` is saved; nothing is written if `nothing`.
- `logfile::AbstractString`: Path of the log file, defaulting to `MirrorCorrection.log`
  next to `outfile` (or in the current directory if `outfile` is `nothing`),
  it can be inspected during the run to track which tasks have completed and to estimate the time remaining.

!!! note
	Other method-specific keyword arguments are not documented here; see `/examples` for
	detailed usage.

Returns a `NamedTuple` with `Uwithlr` (the corrected interaction as a `DirectInterAction`,
to be used in later calculations), `Umirror` (the corrected real-space interaction) and
`Ushort` (`U` with the periodic tail removed only); the `_optim` variants additionally
return the fitted parameters.
"""
function MirrorCorrection(sym::Symbol, U::AbstractHR,
    kgrid::Union{MonkhorstPack,RedKgrid,AbstractVector{<:Integer},NTuple{3,<:Integer}},
    lattice::Union{Lattice,AbstractMatrix{<:Real}}, orblocat::AbstractVector{V};
    outfile::Union{Nothing,AbstractString}=nothing,
    logfile::AbstractString=isnothing(outfile) ? "./MirrorCorrection.log" : joinpath(dirname(outfile), "MirrorCorrection.log"),
    kwargs...,
) where {V}

    logio = get_logio(logfile)
    original_stdout = stdout
    redirect_stdout(logio)
    _exit_flush(logio)

    # preprocessing
    U = HRFull(U)
    if kgrid isa AbstractVector || kgrid isa NTuple
        kgrid = RedKgrid(MonkhorstPack(kgrid))
    else
        kgrid = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
    end
    lattice = Lattice(parent(lattice))
    orblocat = collect(orblocat)
    if V <: CartesianCoordinates
        orblocat = map(orblocat) do r
            lattice \ r
        end
    elseif V <: AbstractVector
        orblocat = map(orblocat) do r
            ReducedCoordinates(r)
        end
    else
        error("orblocat should be a vector of vector!")
    end

    # calculate
    result = MirrorCorrection(Val(sym), U, kgrid, lattice, orblocat; kwargs...)

    redirect_stdout(original_stdout)
    close(logio)

    if !isnothing(outfile)
        mkpath(dirname(outfile))
        jldsave(outfile; result...)
    end

    return result
end

# 虽然增加高斯展宽后的库伦势在实空间会很快趋于点电荷的结果，但是转换至有限大周期体系中的形式后，高斯展宽参数会体现出较为明显的影响，这为依据远距离的样本点进行参数优化提供了理论依据。
# 2D的要计算远一点。wannier基在非周期方向存在位移的情况下，KeldyshGauss也可以描述，但是需要更远的距离。
# 对比MoS2的结果，wannier基局域性更好后，Mirrorcorrection更难收敛。
