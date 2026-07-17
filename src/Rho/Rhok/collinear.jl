abstract type Rhok_collinear{RhoT} <: AbstractRhok{RhoT} end
struct _Rhok_collinear{WaveT <: Number, RhoT,
	ftT_up <: AbstractFFTs.Plan{RhoT}, iftT_up <: AbstractFFTs.Plan{RhoT},
	ftT_dn <: AbstractFFTs.Plan{RhoT}, iftT_dn <: AbstractFFTs.Plan{RhoT},
} <: Rhok_collinear{RhoT}
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
	nw_up::Int
	nw_dn::Int
	upindex::Vector{Int}
	dnindex::Vector{Int}
	ρk_up::_Rhok_wospin{WaveT, RhoT, ftT_up, iftT_up}
	ρk_dn::_Rhok_wospin{WaveT, RhoT, ftT_dn, iftT_dn}
end
function Base.show(io::IO, ρk::Rhok_collinear)
	println(io, "ρk_collinear object:")
	println(io, "  nw = $(ρk.nw), nk = $(ρk.nk), nq = $(ρk.nq)")
	println(io, "  nG_max = $(ρk.nG_max)")
	println(io, "  fft_grid = $(ρk.fft_grid)")
	return nothing
end
@inbounds function (ρk::Rhok_collinear{RhoT})(
	ρk_G::AbstractVector{RhoT},
	wavek3D::Array{RhoT, 3}, wavekq3D::Array{RhoT, 3},
	k::Int, i::Int, j::Int, q::Int,
)::Vector{RhoT} where {RhoT}
	ispin = i <= ρk.nw_up ? :up : :down
	jspin = j <= ρk.nw_up ? :up : :down
	if ispin ≠ jspin
		return fill!(ρk_G, zero(RhoT))
	elseif ispin == :up
		return ρk.ρk_up(ρk_G, wavek3D, wavekq3D, k, i, j, q)
	else
		return ρk.ρk_dn(ρk_G, wavek3D, wavekq3D, k, i - ρk.nw_up, j - ρk.nw_up, q)
	end
end
function Rhok_collinear(wave_up::WannierPWKgrid{WaveT}, wave_dn::WannierPWKgrid{WaveT},
	qgrid::RedKgrid, Gidx::Vector{Vector{Int}}, Glist::Vector{ReducedCoordinates{Int}};
	fft_grid = wave_up.fft_grid,
) where {WaveT}

	nq = length(qgrid)
	q0_idx = findfirst(iszero, qgrid)

	kgrid_up = wave_up.kgrid.kgrid_size
	kgrid_dn = wave_dn.kgrid.kgrid_size
	if kgrid_up ≠ kgrid_dn
		error("Mismatched kgrid of wave.")
	end
	nk = length(wave_up.kgrid)

	nw_up = wave_up.nw
	nw_dn = wave_dn.nw
	nw = nw_up + nw_dn
	upindex = collect(1:nw_up)
	dnindex = collect((nw_up+1):nw)

	fft_grid = (fft_grid[1], fft_grid[2], fft_grid[3])

	ρk_up = Rhok_wospin(wave_up, qgrid, Gidx, Glist; fft_grid)
	ρk_dn = Rhok_wospin(wave_dn, qgrid, Gidx, Glist; fft_grid)

	RhoT = complex(WaveT)
	return _Rhok_collinear{WaveT, RhoT,
		typeof(ρk_up.ft_wave), typeof(ρk_up.ift_wave),
		typeof(ρk_dn.ft_wave), typeof(ρk_dn.ift_wave)}(
		nw, nk, nq,
		ρk_up.kgrid, qgrid,
		ρk_up.nG, ρk_up.nG_max, Gidx,
		q0_idx, ρk_up.G0idx_q0,
		Glist,
		fft_grid,
		wave_up.period, wave_up.lattice, wave_up.rlattice, wave_up.ldot, wave_up.rldot, wave_up.Ω, wave_up.rΩ,
		nw_up, nw_dn, upindex, dnindex, ρk_up, ρk_dn,
	)
end


