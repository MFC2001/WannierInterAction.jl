"""
	Rhok_wospin{RhoT} <: AbstractRhok{RhoT}

Abstract type for the spinless (without spin) density matrix in reciprocal space,
`ρ_k(q, G) = ⟨u_{k}| e^{-i(q+G)·r} |u_{k+q}⟩`.

This is the scalar (non-spinor) version used as a building block for both
collinear and noncollinear spin cases. It performs FFT-based computation of
the pair density via:

1. Expand wave functions ``u_{nk}`` and ``u_{mk+q}`` onto a 3D FFT grid.
2. Inverse FFT → real-space product → forward FFT.
3. Extract the relevant G-components weighted by ``1/(|q+G|·N_k)``.

See also: [`Rhok_collinear`](@ref), [`Rhok_noncollinear`](@ref).
"""
abstract type Rhok_wospin{RhoT} <: AbstractRhok{RhoT} end
struct _Rhok_wospin{WaveT <: Number, RhoT, ftT <: AbstractFFTs.Plan{RhoT}, iftT <: AbstractFFTs.Plan{RhoT}} <: Rhok_wospin{RhoT}
	#public
	nw::Int
	nk::Int
	nq::Int
	kgrid::RedKgrid
	qgrid::RedKgrid
	nG::Vector{Int} # length.(Gidx)
	nG_max::Int # greater than nGidx_negq0
	Gidx::Vector{Vector{Int}} # [iq][iG]
	q0_idx::Int
	G0idx_q0::Int # used for head term correction
	Glist::Vector{ReducedCoordinates{Int}}
	fft_grid::NTuple{3, Int}
	period::Vec3{Bool}
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	ldot::Mat3{Float64}
	rldot::Mat3{Float64}
	Ω::Float64
	rΩ::Float64
	#private
	wavek::Vector{Matrix{WaveT}}
	nG_wave::Vector{Int} # length.(Gidx_wave)
	nG_wave_max::Int
	Gidx_wave::Vector{Vector{Int}} # [iq][iG].
	Glist_wave::Vector{ReducedCoordinates{Int}}
	wavek_3DGmap::Vector{Int} # linear map of Glist_wave #线性索引远优于笛卡尔索引。
	wavek_3DGmap_ik::Vector{Vector{Int}} # linear map of Glist_wave of ik
	rhok_3DGmap::Vector{Int} # linear map of Glist_rhok
	ft_wave::ftT
	ift_wave::iftT
	kq_addmap::Matrix{Int} # [ik, iq], map to k
	kq_Gwind::Matrix{ReducedCoordinates{Int}} # k[ik]+q[iq] = k[kq_addmap[ik,iq]] + kq_Gwind[ik,iq]
	kq_Gwind_idx::Matrix{Int} # Gwind_Grho_addmap[kq_Gwind_idx[ik, iq], Gidx[iq][iG]]
	Grho_Gwind_addmap::Matrix{Int} # [iGρ, iGwind], map to Grho, missing G is 0.
	fft_factor::Float64
end
function Base.show(io::IO, ρk::Rhok_wospin)
	println(io, "ρk_wospin object:")
	println(io, "  nw = $(ρk.nw), nk = $(ρk.nk), nq = $(ρk.nq)")
	println(io, "  nG_max = $(ρk.nG_max)")
	println(io, "  fft_grid = $(ρk.fft_grid)")
	return nothing
end

@inbounds function (ρk::Rhok_wospin{RhoT})(
	ρk_G::AbstractVector{RhoT},
	wavek3D::AbstractArray{RhoT, 3}, wavekq3D::AbstractArray{RhoT, 3},
	k::Int, i::Int, j::Int, q::Int,
)::Vector{RhoT} where {RhoT}

	wavek = view(ρk.wavek[k], :, i)
	wavek3Dmap = ρk.wavek_3DGmap_ik[k]
	_get_wave_3D!(wavek3D, wavek, wavek3Dmap)
	kq = ρk.kq_addmap[k, q]
	wavekq = view(ρk.wavek[kq], :, j)
	wavekq3Dmap = ρk.wavek_3DGmap_ik[kq]
	_get_wave_3D!(wavekq3D, wavekq, wavekq3Dmap)
	Gidx_q = ρk.Gidx[q]
	Gwind_idx = ρk.kq_Gwind_idx[k, q]
	GGidx = view(ρk.Grho_Gwind_addmap, Gidx_q, Gwind_idx)

	ρk.ift_wave * wavek3D
	ρk.ift_wave * wavekq3D
	wavekq3D .*= conj!(wavek3D)
	ρk.ft_wave * wavekq3D

	for (iG, GG) in enumerate(GGidx)
		ρk_G[iG] = iszero(GG) ? zero(RhoT) : wavekq3D[ρk.rhok_3DGmap[GG]] * ρk.fft_factor
	end

	return ρk_G
end
#kgrid and qgrid have the same size, but qgrid must be Γ-centred.
function Rhok_wospin(wave::WannierPWKgrid{WaveT},
	qgrid::RedKgrid, Gidx::Vector{Vector{Int}}, Glist::Vector{ReducedCoordinates{Int}};
	fft_grid = wave.fft_grid,
) where {WaveT}

	nq = length(qgrid)
	nk = wave.nk
	q0_idx = findfirst(iszero, qgrid)

	fft_grid = (fft_grid[1], fft_grid[2], fft_grid[3])
	empty_fftgrid = Array{complex(WaveT), 3}(undef, fft_grid)
	ift_wave = plan_bfft!(empty_fftgrid; flags = FFTW.ESTIMATE)
	ft_wave = plan_fft!(empty_fftgrid; flags = FFTW.ESTIMATE)

	pgi = PeriodicGridIndex(fft_grid)
	wavek_3DGmap = pgi.(Val(:linear), wave.Glist)
	rhok_3DGmap = pgi.(Val(:linear), Glist)
	wavek_3DGmap_ik = map(wave.Gidx) do Gidx_ik
		wavek_3DGmap[Gidx_ik]
	end

	G0idx_q0 = findfirst(iG -> iszero(Glist[iG]), Gidx[q0_idx])

	(kq_addmap, kq_Gwind) = _kgrid_map(wave.kgrid.kdirect, wave.kgrid.kdirect, qgrid.kdirect, +) #[nk, nq]
	any(iszero, kq_addmap) && error("Mismatched kgrid of wave and qgrid of epsmat.")

	kq_Gwind_unique = unique(kq_Gwind)
	kq_Gwind_unique_Dict = Dict{ReducedCoordinates{Int}, Int}(G => i for (i, G) in enumerate(kq_Gwind_unique))
	kq_Gwind_idx = map(G -> kq_Gwind_unique_Dict[G], kq_Gwind)

	Grho_Gwind_addmap = _Glist_map(Glist, Glist, kq_Gwind_unique, +) #[nGρ, nGwind]

	nG = length.(Gidx)
	nG_max = maximum(nG)

	fft_factor = 1.0 / prod(fft_grid)

	RhoT = complex(WaveT)
	return _Rhok_wospin{WaveT, RhoT, typeof(ft_wave), typeof(ift_wave)}(
		wave.nw, nk, nq,
		wave.kgrid, qgrid,
		nG, nG_max, Gidx,
		q0_idx, G0idx_q0,
		Glist,
		fft_grid, wave.period, wave.lattice, wave.rlattice, wave.ldot, wave.rldot, wave.Ω, wave.rΩ,
		wave.components, wave.nG, wave.nG_max, wave.Gidx, wave.Glist,
		wavek_3DGmap, wavek_3DGmap_ik, rhok_3DGmap,
		ft_wave, ift_wave,
		kq_addmap, kq_Gwind, kq_Gwind_idx,
		Grho_Gwind_addmap, fft_factor,
	)
end


