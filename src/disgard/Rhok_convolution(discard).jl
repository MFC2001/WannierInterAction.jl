struct ρk_wospin{waveT <: Number} <: Abstractρk
	#public
	nw::Int
	nk::Int
	nq::Int
	nGidx::Vector{Int} # length.(Gidx)
	nGidx_max::Int
	Gidx::Vector{Vector{Int}} # [iq][iG].
	kgrid::RedKgrid
	qgrid::RedKgrid
	smallq::ReducedCoordinates{Float64}
	smallq_idx::Int
	G0idx_smallq::Int # used for head term correction
	Glist::Vector{ReducedCoordinates{Int}}
	Glist_wave::Vector{ReducedCoordinates{Int}}
	qG_norm::Vector{Vector{Float64}} # [iq][iG] |q+G|
	lattice::Lattice{Float64}
	rlattice::ReciprocalLattice{Float64}
	ldot::Mat3{Float64}
	rldot::Mat3{Float64}
	Ω::Float64
	rΩ::Float64
	#private
	wavek::Array{waveT}
	wavek_smallq::Array{waveT} #For small q.
	sumGidx_wave::Base.OneTo{Int}
	kq_addmap::Matrix{Int} # [ik, iq], map to k
	kq_Gwind::Matrix{ReducedCoordinates{Int}} # k[ik]+q[iq] = k[kq_addmap[ik,iq]] + kq_Gwind[ik,iq]
	kq_Gwind_idx::Matrix{Int} # Gwind_Grho_addmap[kq_Gwind_idx[ik, iq], Gidx[iq][iG]]
	Gwave_Grho_addmap::Matrix{Int} # [iGwave, iGρ], map to Gwave, missing G is 0.
	Gwind_Grho_addmap::Matrix{Int} # [iGwind, iGρ], map to Grho, missing G is 0.
	Gidx_mat::Matrix{Int} # Matrix form of Gidx for faster access [nG, nq], missing G is 0.
	qGN_norm_inv::Matrix{Float64} # 1/|q+G|nk, missing G is 0.
end
@inline function (ρk::ρk_wospin{waveT})(k::Int, G::Int, i::Int, j::Int, q::Int)::waveT where {waveT}
	@views @inbounds begin
		Gwind_idx = ρk.kq_Gwind_idx[k, q]
		Grho_idx = ρk.Gidx_mat[G, q]
		iszero(Grho_idx) && return zero(waveT)
		Gidx = ρk.Gwind_Grho_addmap[Gwind_idx, Grho_idx]
		iszero(Gidx) && return zero(waveT)

		Gmap = ρk.Gwave_Grho_addmap[:, Gidx]

		mask = .!iszero.(Gmap)
		any(mask) || return zero(waveT)

		wavek = ρk.wavek[:, i, k]
		if q == ρk.smallq_idx
			wavekq = ρk.wavek_smallq[:, j, k]
		else
			wavekq = ρk.wavek[:, j, ρk.kq_addmap[k, q]]
		end

		wavek_vals = wavek[mask]
		wavekq_vals = wavekq[Gmap[mask]]

		return (wavek_vals ⋅ wavekq_vals) * ρk.qGN_norm_inv[G, q]
	end
end
@inline function (ρk::ρk_wospin{waveT})(k::Int, i::Int, j::Int, q::Int)::Vector{waveT} where {waveT}
	@views @inbounds begin
		wavek = ρk.wavek[:, i, k]
		if q == ρk.smallq_idx
			wavekq = ρk.wavek_smallq[:, j, k]
		else
			wavekq = ρk.wavek[:, j, ρk.kq_addmap[k, q]]
		end

		ρk_G = Vector{waveT}(undef, ρk.nGidx[q])

		Gwind_idx = ρk.kq_Gwind_idx[k, q]
		for G in Base.OneTo(ρk.nGidx[q])
			Grho_idx = ρk.Gidx_mat[G, q]
			iszero(Grho_idx) && begin
				ρk_G[G] = zero(waveT)
				continue
			end
			Gidx = ρk.Gwind_Grho_addmap[Gwind_idx, Grho_idx]
			iszero(Gidx) && begin
				ρk_G[G] = zero(waveT)
				continue
			end
			Gmap = ρk.Gwave_Grho_addmap[:, Gidx]
			mask = .!iszero.(Gmap)
			any(mask) || begin
				ρk_G[G] = zero(waveT)
				continue
			end
			wavek_vals = wavek[mask]
			wavekq_vals = wavekq[Gmap[mask]]

			ρk_G[G] = (wavek_vals ⋅ wavekq_vals) * ρk.qGN_norm_inv[G, q]
		end

		return ρk_G
	end
end
#kgrid and qgrid have the same size, but qgrid must be Γ-centred.
function ρk_wospin(wave::WannierWaveKgrid,
	qgrid::RedKgrid, Gidx::Vector{Vector{Int}}, Glist::Vector{ReducedCoordinates{Int}};
	smallq = ReducedCoordinates(0.0, 0.0, 0.0), smallq_idx = findfirst(iszero, qgrid),
	wavek_smallq::Union{Array{<:Number, 3}, Nothing} = nothing,
	qG_norm::Union{Vector{Vector{Float64}}, Nothing} = nothing,
	wavek_smallq_tol::Real = 1e-8,
)
	G0idx_smallq = findfirst(iG -> iszero(Glist[iG]), Gidx[smallq_idx])

	(kq_addmap, kq_Gwind) = _kgrid_map(wave.kgrid.kdirect, wave.kgrid.kdirect, qgrid.kdirect, +) #[nk, nq]
	any(iszero, kq_addmap) && error("Mismatched kgrid of wave and qgrid of epsmat.")

	Gwave_Grho_addmap = _Glist_map(wave.Glist, wave.Glist, Glist, +) #[nGwave, nGρ]

	kq_Gwind_unique = unique(kq_Gwind)
	kq_Gwind_unique_Dict = Dict{ReducedCoordinates{Int}, Int}(G => i for (i, G) in enumerate(kq_Gwind_unique))
	kq_Gwind_idx = getindex.(Ref(kq_Gwind_unique_Dict), kq_Gwind)

	Gwind_Grho_addmap = _Glist_map(Glist, kq_Gwind_unique, Glist, +) #[nGwind, nGρ]

	# wavek_smallq
	println("calculate wave_smallq")
	if isnothing(wavek_smallq)
		if iszero(smallq)
			wavek_smallq = wave.wavek
		else
			wavek_smallq = make_wave_smallq(wave.wavek, wave.kgrid, wave.Glist, smallq; tol = wavek_smallq_tol)
		end
	end
	println("calculate wave_smallq end")

	nq = length(qgrid)
	(nG, nw, nk) = size(wave.wavek)
	nGidx = length.(Gidx)
	nGidx_max = maximum(nGidx)

	Gidx_mat = Matrix{Int}(undef, nGidx_max, nq)
	for iq in Base.OneTo(nq)
		@inbounds begin
			nGidx_iq = nGidx[iq]
			Gidx_mat[1:nGidx_iq, iq] .= Gidx[iq]
			Gidx_mat[nGidx_iq+1:nGidx_max, iq] .= 0
		end
	end

	# qG_norm
	rldot = wave.rldot
	if isnothing(qG_norm)
		qG_norm = Vector{Vector{Float64}}(undef, nq)
		for iq in Base.OneTo(nq)
			if iq == smallq_idx
				q = smallq
			else
				q = qgrid[iq]
			end
			qG_norm[iq] = map(Gidx[iq]) do iG
				qG = q + Glist[iG]
				sqrt(qG ⋅ (rldot * qG))
			end
		end
	end

	qGN_norm_inv = Matrix{Float64}(undef, nGidx_max, nq)
	Threads.@threads :static for iq in Base.OneTo(nq)
		@inbounds begin
			qGN_norm_inv[1:nGidx[iq], iq] .= 1.0 ./ (qG_norm[iq] * nk)
			qGN_norm_inv[nGidx[iq]+1:nGidx_max, iq] .= 0.0
		end
	end

	return ρk_wospin(
		nw, nk, nq, nGidx, nGidx_max, Gidx,
		wave.kgrid, qgrid,
		smallq, smallq_idx, G0idx_smallq,
		Glist, wave.Glist, qG_norm,
		wave.lattice, wave.rlattice, wave.ldot, wave.rldot, wave.Ω, wave.rΩ,
		wave.wavek, wavek_smallq,
		Base.OneTo(nG),
		kq_addmap, kq_Gwind, kq_Gwind_idx,
		Gwave_Grho_addmap,
		Gwind_Grho_addmap,
		Gidx_mat, qGN_norm_inv,
	)
end


