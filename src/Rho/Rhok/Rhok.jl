export Rhok
"""
	abstract type AbstractRhok{RhoT <: Number} end
	(ρk::AbstractRhok{RhoT})(k::Int, i::Int, j::Int, q::Int)::Vector{RhoT}
	(ρk::AbstractRhok{RhoT})(ρk_G::AbstractVector{RhoT}, k::Int, i::Int, j::Int, q::Int)::Vector{RhoT}
	(ρk::AbstractRhok{RhoT})(wavek3D::AbstractArray{RhoT, 3}, wavekq3D::AbstractArray{RhoT, 3}, k::Int, i::Int, j::Int, q::Int)::Vector{RhoT}
"""
abstract type AbstractRhok{RhoT<:Number} end
@inline function (ρk::AbstractRhok{RhoT})(k::Integer, i::Integer, j::Integer, q::Integer)::Vector{RhoT} where {RhoT}
    return ρk(Int(k), Int(i), Int(j), Int(q))
end
@inline function (ρk::AbstractRhok{RhoT})(k::Int, i::Int, j::Int, q::Int)::Vector{RhoT} where {RhoT}
    nG_q = ρk.nG[q]
    ρk_G = Vector{RhoT}(undef, nG_q)
    wavek3D = Array{RhoT,3}(undef, ρk.fft_grid)
    wavekq3D = Array{RhoT,3}(undef, ρk.fft_grid)
    return ρk(ρk_G, wavek3D, wavekq3D, k, i, j, q)
end
@inline function (ρk::AbstractRhok{RhoT})(ρk_G::AbstractVector{RhoT}, k::Int, i::Int, j::Int, q::Int)::Vector{RhoT} where {RhoT}
    wavek3D = Array{RhoT,3}(undef, ρk.fft_grid)
    wavekq3D = Array{RhoT,3}(undef, ρk.fft_grid)
    return ρk(ρk_G, wavek3D, wavekq3D, k, i, j, q)
end
@inline function (ρk::AbstractRhok{RhoT})(wavek3D::AbstractArray{RhoT,3}, wavekq3D::AbstractArray{RhoT,3},
    k::Int, i::Int, j::Int, q::Int)::Vector{RhoT} where {RhoT}
    nG_q = ρk.nG[q]
    ρk_G = Vector{RhoT}(undef, nG_q)
    return ρk(ρk_G, wavek3D, wavekq3D, k, i, j, q)
end

include("./wospin.jl")
include("./collinear.jl")
include("./noncollinear.jl")

"""
    Rhok(wave::WannierPWKgrid; ecutrho::Real = -1, fft_grid = wave.fft_grid,
         spin_atol::Real = 1e-6, upindex::AbstractVector{<:Integer} = Int[],
         dnindex::AbstractVector{<:Integer} = Int[])
	Rhok(eps::DielectricPWBasis, wave::WannierPWKgrid; ecutrho::Real = -1, fft_grid = wave.fft_grid,
         spin_atol::Real = 1e-6, upindex::AbstractVector{<:Integer} = Int[],
         dnindex::AbstractVector{<:Integer} = Int[])

Compute the density in reciprocal space (`ρ_{G}`) for a set of Wannier functions.

This function supports both spinless and noncollinear calculations.

# Arguments
- `eps`: The dielectric function in plane-wave basis.
- `wave`: The Wannier functions represented in plane-wave basis.
- `ecutrho`: Energy cutoff for the density. If not provided (or set to a negative value),
  it defaults to four times the truncated energy of the wave functions.
- `fft_grid`: FFT grid used during the computation. Defaults to `wave.fft_grid`.
- `spin_atol`: Absolute tolerance for automatic detection of spin polarization of Wannier functions.
- `upindex`: Indices of Wannier functions with spin-up character. If empty, they are determined automatically.
- `dnindex`: Indices of Wannier functions with spin-down character. If empty, they are determined automatically.

# Notes
- If `eps` is provided, the ordering of the low-energy basis set of the density is rearranged to match that of the dielectric function.
- When the Wannier functions have a clear spin orientation, specifying the `upindex`/`dnindex` explicitly can reduce computational overhead.
- The function automatically checks the spin polarization of the Wannier functions based on `spin_atol` even though the indices are manually provided.
"""
function Rhok(wave::WannierPWKgrid;
    ecutrho::Real=-1, fft_grid=wave.fft_grid,
    spin_atol::Real=1e-6,
    upindex::AbstractVector{<:Integer}=Int[],
    dnindex::AbstractVector{<:Integer}=Int[],
)

    normcutrho = ecutrho <= 0 ? -1 : sqrt(2 * mₑ * qₑ * ecutrho / ħ^2) * 0.1

    rldot = wave.rldot

    if normcutrho <= 0
        normcutwave = sqrt(maximum(G -> G ⋅ (rldot * G), wave.Glist))
        normcutrho = ceil(2 * normcutwave; digits=1)
    end

    qgrid = RedKgrid(MonkhorstPack(wave.kgrid.kgrid_size))
    Glist_rho, Gidx_rho = build_Glist_rho(qgrid, fft_grid, rldot; normcutrho)

    @info "Rhok calculation parameters" normcutrho=normcutrho nGlist_rho=length(Glist_rho) nG_max_rho=maximum(length.(Gidx_rho))

    if wave.nspinor == 1
        return Rhok_wospin(wave, qgrid, Gidx_rho, Glist_rho; fft_grid)
    elseif wave.nspinor == 2
        @assert all(iw->iw ∈ 1:wave.nw, upindex) "Wrong upindex, the index should be in 1:nw."
        @assert all(iw->iw ∈ 1:wave.nw, dnindex) "Wrong upindex, the index should be in 1:nw."
        @assert isdisjoint(upindex, dnindex) "Wrong upindex or dnindex, they can only have different index."
        spinindex = wave(:spin_index; atol=spin_atol)
        upindex = union(spinindex[1], upindex)
        dnindex = union(spinindex[2], dnindex)
        @assert isdisjoint(upindex, dnindex) "May be you need to exchange `upindex` and `dnindex`."
        sort!(upindex)
        sort!(dnindex)
        return Rhok_noncollinear(wave, qgrid, Gidx_rho, Glist_rho; fft_grid, upindex, dnindex)
    end
end
function Rhok(eps::DielectricPWBasis, wave::WannierPWKgrid;
    ecutrho::Real=-1, fft_grid=wave.fft_grid,
    spin_atol::Real=1e-6,
    upindex::AbstractVector{<:Integer}=Int[],
    dnindex::AbstractVector{<:Integer}=Int[],
)

    normcutrho = ecutrho <= 0 ? -1 : sqrt(2 * mₑ * qₑ * ecutrho / ħ^2) * 0.1

    rldot = wave.rldot

    if normcutrho <= 0
        normcutwave = sqrt(maximum(G -> G ⋅ (rldot * G), wave.Glist))
        normcutrho = ceil(2 * normcutwave; digits=1)
    end

    Glist_rho, Gidx_rho = build_Glist_rho(eps, fft_grid, rldot; normcutrho)

    @info "Rhok calculation parameters" normcutrho=normcutrho nGlist_rho=length(Glist_rho) nG_max_rho=maximum(length.(Gidx_rho))

    if wave.nspinor == 1
        return Rhok_wospin(wave, eps.qgrid, Gidx_rho, Glist_rho; fft_grid)
    elseif wave.nspinor == 2
        @assert all(iw->iw ∈ 1:wave.nw, upindex) "Wrong upindex, the index should be in 1:nw."
        @assert all(iw->iw ∈ 1:wave.nw, dnindex) "Wrong upindex, the index should be in 1:nw."
        @assert isdisjoint(upindex, dnindex) "Wrong upindex or dnindex, they can only have different index."
        spinindex = wave(:spin_index; atol=spin_atol)
        upindex = union(spinindex[1], upindex)
        dnindex = union(spinindex[2], dnindex)
        @assert isdisjoint(upindex, dnindex) "May be you need to exchange `upindex` and `dnindex`."
        sort!(upindex)
        sort!(dnindex)
        return Rhok_noncollinear(wave, eps.qgrid, Gidx_rho, Glist_rho; fft_grid, upindex, dnindex)
    end
end
"""
    Rhok(wave_up::WannierPWKgrid, wave_dn::WannierPWKgrid;
         ecutrho::Real = -1, fft_grid = wave_up.fft_grid)
	Rhok(eps::DielectricPWBasis, wave_up::WannierPWKgrid, wave_dn::WannierPWKgrid;
         ecutrho::Real = -1, fft_grid = wave_up.fft_grid)

Compute the density in reciprocal space (`ρ_{G}`) for a set of Wannier functions.

This function supports collinear calculations.

# Arguments
- `eps`: The dielectric function in plane-wave basis.
- `wave_up`: The spin-up Wannier functions represented in plane-wave basis.
- `wave_dn`: The spin-down Wannier functions represented in plane-wave basis.
- `ecutrho`: Energy cutoff for the density. If not provided (or set to a negative value),
  it defaults to four times the truncated energy of the wave functions.
- `fft_grid`: FFT grid used during the computation. Defaults to `wave.fft_grid`.

# Notes
- If `eps` is provided, the ordering of the low-energy basis set of the density is rearranged to match that of the dielectric function.
"""
function Rhok(wave_up::WannierPWKgrid, wave_dn::WannierPWKgrid;
    ecutrho::Real=-1, fft_grid=wave_up.fft_grid)

    if wave_up.nspinor ≠ 1 || wave_dn.nspinor ≠ 1
        error("Wrong wave, nspinor should be 1 for Rhok_collinear!")
    end

    normcutrho = ecutrho <= 0 ? -1 : sqrt(2 * mₑ * qₑ * ecutrho / ħ^2) * 0.1

    rldot = wave_up.rldot

    if normcutrho <= 0
        normcutwave = sqrt(maximum(G -> G ⋅ (rldot * G), wave_up.Glist))
        normcutrho = ceil(2 * normcutwave; digits=1)
    end

    qgrid = RedKgrid(MonkhorstPack(wave_up.kgrid.kgrid_size))
    Glist_rho, Gidx_rho = build_Glist_rho(qgrid, fft_grid, rldot; normcutrho)

    @info "Rhok calculation parameters" normcutrho=normcutrho nGlist_rho=length(Glist_rho) nG_max_rho=maximum(length.(Gidx_rho))

    return Rhok_collinear(wave_up, wave_dn, qgrid, Gidx_rho, Glist_rho; fft_grid)
end
function Rhok(eps::DielectricPWBasis, wave_up::WannierPWKgrid, wave_dn::WannierPWKgrid;
    ecutrho::Real=-1, fft_grid=wave_up.fft_grid)

    if wave_up.nspinor ≠ 1 || wave_dn.nspinor ≠ 1
        error("Wrong wave, nspinor should be 1 for Rhok_collinear!")
    end

    normcutrho = ecutrho <= 0 ? -1 : sqrt(2 * mₑ * qₑ * ecutrho / ħ^2) * 0.1

    rldot = wave_up.rldot

    if normcutrho <= 0
        normcutwave = sqrt(maximum(G -> G ⋅ (rldot * G), wave_up.Glist))
        normcutrho = ceil(2 * normcutwave; digits=1)
    end

    Glist_rho, Gidx_rho = build_Glist_rho(eps, fft_grid, rldot; normcutrho)

    @info "Rhok calculation parameters" normcutrho=normcutrho nGlist_rho=length(Glist_rho) nG_max_rho=maximum(length.(Gidx_rho))

    return Rhok_collinear(wave_up, wave_dn, eps.qgrid, Gidx_rho, Glist_rho; fft_grid)
end


