
abstract type WannierWaveKgridShiftInterpolation{WaveT} end

function waveik_interpolation_pre(wave::WannierPWKgrid{WaveT}) where {WaveT}
	nk = wave.nk
	Gidx = wave.Gidx
	Glist = wave.Glist
	Glists = Vector{Vector{ReducedCoordinates{Int}}}(undef, nk)
	for ik in 1:nk
		Glists[ik] = Glist[Gidx[ik]]
	end

	kgrid = deepcopy(wave.kgrid)
	wavek_normalize!(Glists, kgrid)
	(Glist, Gidx) = merge_Glist(Glists...; rldot = wave.rldot, Gidx = true)

	return WannierPWKgrid{WaveT}(
		nk, wave.nw, wave.nspin, wave.nG, wave.nG_max, kgrid, Gidx, Glist, wave.components, wave.fft_grid, wave.Gmap,
		wave.period, wave.lattice, wave.rlattice, wave.ldot, wave.rldot, wave.Ω, wave.rΩ,
	)
end
function WannierWaveKgridShiftInterpolation(wave::WannierPWKgrid{WaveT}) where {WaveT}
	if wave.nspin == 1 || wave.nspin == 2
		return WannierWaveKgridShiftInterpolationScalar(wave)
	elseif wave.nspin == 4
		return WannierWaveKgridShiftInterpolationSpinor(wave)
	else
		error("Wrong wave!")
	end
end
struct WannierWaveKgridShiftInterpolationScalar{WaveT, etpT} <: WannierWaveKgridShiftInterpolation{WaveT}
	nk::Int
	nw::Int
	nG::Vector{Int}
	kgrid::RedKgrid
	Gidx::Vector{Vector{Int}}
	Glist::Vector{ReducedCoordinates{Int}}
	mapk::Vector{NTuple{3, Int}}
	mapG::Vector{NTuple{3, Int}}
	kshift_interpolation::NTuple{3, Float64}
	WT::Type{WaveT}
	extrapolation::Vector{etpT}
end
function (etp::WannierWaveKgridShiftInterpolationScalar{WaveT, etpT})(q::ReducedCoordinates{Float64}) where {WaveT, etpT}
	(nkx, nky, nkz) = etp.kgrid.kgrid_size
	(kx_shift, ky_shift, kz_shift) = etp.kshift_interpolation
	mapkq = Vector{NTuple{3, Float64}}(undef, etp.nk)
	for ik in Base.OneTo(etp.nk)
		kq = etp.kgrid[ik] + q
		kqx_idx = (kq[1] + kx_shift) * nkx + 1
		kqy_idx = (kq[2] + ky_shift) * nky + 1
		kqz_idx = (kq[3] + kz_shift) * nkz + 1
		mapkq[ik] = (kqx_idx, kqy_idx, kqz_idx)
	end
	wavekq = [Matrix{WaveT}(undef, etp.nG[ik], etp.nw) for ik in 1:etp.nk]
	Threads.@threads for I in CartesianIndices((etp.nw, etp.nk))
		(iw, ik) = Tuple(I)
		(kqx_idx, kqy_idx, kqz_idx) = mapkq[ik]
		wavekq_iw_ik = view(wavekq[ik], :, iw)
		for (iG, iG_all) in enumerate(etp.Gidx[ik])
			(Gx_idx, Gy_idx, Gz_idx) = etp.mapG[iG_all]
			wavekq_iw_ik[iG] = etp.extrapolation[iw](
				Gx_idx + kqx_idx,
				Gy_idx + kqy_idx,
				Gz_idx + kqz_idx,
			)
		end
		wave_norm = norm(wavekq_iw_ik)
		if !isapprox(wave_norm, 1; atol = 1e-6)
			wavekq_iw_ik ./= wave_norm
		end
	end
	return wavekq
end
function (etp::WannierWaveKgridShiftInterpolationScalar{WaveT, etpT})(wavekq::AbstractVector{<:AbstractMatrix{WaveT}}, q::ReducedCoordinates{Float64}) where {WaveT, etpT}
	(nkx, nky, nkz) = etp.kgrid.kgrid_size
	(kx_shift, ky_shift, kz_shift) = etp.kshift_interpolation
	mapkq = Vector{NTuple{3, Float64}}(undef, etp.nk)
	for ik in Base.OneTo(etp.nk)
		kq = etp.kgrid[ik] + q
		kqx_idx = (kq[1] + kx_shift) * nkx + 1
		kqy_idx = (kq[2] + ky_shift) * nky + 1
		kqz_idx = (kq[3] + kz_shift) * nkz + 1
		mapkq[ik] = (kqx_idx, kqy_idx, kqz_idx)
	end
	Threads.@threads for I in CartesianIndices((etp.nw, etp.nk))
		(iw, ik) = Tuple(I)
		(kqx_idx, kqy_idx, kqz_idx) = mapkq[ik]
		wavekq_iw_ik = view(wavekq[ik], :, iw)
		for (iG, iG_all) in enumerate(etp.Gidx[ik])
			(Gx_idx, Gy_idx, Gz_idx) = etp.mapG[iG_all]
			wavekq_iw_ik[iG] = etp.extrapolation[iw](
				Gx_idx + kqx_idx,
				Gy_idx + kqy_idx,
				Gz_idx + kqz_idx,
			)
		end
		wave_norm = norm(wavekq_iw_ik)
		if !isapprox(wave_norm, 1; atol = 1e-6)
			wavekq_iw_ik ./= wave_norm
		end
	end
	return wavekq
end
function WannierWaveKgridShiftInterpolationScalar(wave::WannierPWKgrid{WaveT}) where {WaveT}

	nk = wave.nk
	nw = wave.nw
	nG = wave.nG
	kgrid = wave.kgrid
	Gidx = wave.Gidx
	Glist = wave.Glist

	(minGx, maxGx) = extrema(G -> G[1], Glist)
	(minGy, maxGy) = extrema(G -> G[2], Glist)
	(minGz, maxGz) = extrema(G -> G[3], Glist)
	nGx = maxGx - minGx + 1
	nGy = maxGy - minGy + 1
	nGz = maxGz - minGz + 1
	(nkx, nky, nkz) = kgrid.kgrid_size
	if iszero(kgrid.kshift[1])
		kx_shift = isodd(nkx) ? floor(Int, nkx / 2) // nkx : 1 // 2
	else
		kx_shift = isodd(nkx) ? 1 // 2 : (nkx - 1) // (2 * nkx)
	end
	if iszero(kgrid.kshift[2])
		ky_shift = isodd(nky) ? floor(Int, nky / 2) // nky : 1 // 2
	else
		ky_shift = isodd(nky) ? 1 // 2 : (nky - 1) // (2 * nky)
	end
	if iszero(kgrid.kshift[3])
		kz_shift = isodd(nkz) ? floor(Int, nkz / 2) // nkz : 1 // 2
	else
		kz_shift = isodd(nkz) ? 1 // 2 : (nkz - 1) // (2 * nkz)
	end

	mapk = Vector{NTuple{3, Int}}(undef, nk)
	for ik in Base.OneTo(nk)
		k = kgrid[ik]
		kx_idx = Int((k[1] + kx_shift) * nkx) + 1
		ky_idx = Int((k[2] + ky_shift) * nky) + 1
		kz_idx = Int((k[3] + kz_shift) * nkz) + 1
		mapk[ik] = (kx_idx, ky_idx, kz_idx)
	end
	nGlist = length(Glist)
	mapG = Vector{NTuple{3, Int}}(undef, nGlist)
	for iG in Base.OneTo(nGlist)
		G = Glist[iG]
		Gx_idx = (G[1] - minGx) * nkx
		Gy_idx = (G[2] - minGy) * nky
		Gz_idx = (G[3] - minGz) * nkz
		mapG[iG] = (Gx_idx, Gy_idx, Gz_idx)
	end

	nkGx = nGx * nkx
	nkGy = nGy * nky
	nkGz = nGz * nkz

	C = Array{WaveT}(undef, nkGx, nkGy, nkGz)
	etps = map(1:nw) do iw
		fill!(C, zero(WaveT))
		Threads.@threads for ik in Base.OneTo(nk)
			(kx_idx, ky_idx, kz_idx) = mapk[ik]
			wave_ik = view(wave.components[ik], :, iw)
			for (iG, iG_all) in enumerate(Gidx[ik])
				(Gx_idx, Gy_idx, Gz_idx) = mapG[iG_all]
				C[
					Gx_idx+kx_idx,
					Gy_idx+ky_idx,
					Gz_idx+kz_idx,
				] = wave_ik[iG]
			end
		end
		C_itp = interpolate(C, BSpline(Cubic(Flat(OnGrid()))))
		C_etp = extrapolate(C_itp, Flat())
		return C_etp
	end

	return WannierWaveKgridShiftInterpolationScalar(
		nk, nw, nG, kgrid, Gidx, Glist, mapk, mapG,
		Float64.((kx_shift, ky_shift, kz_shift)), WaveT, etps,
	)
end
struct WannierWaveKgridShiftInterpolationSpinor{WaveT, etpT} <: WannierWaveKgridShiftInterpolation{WaveT}
	nk::Int
	nw::Int
	nG::Vector{Int}
	kgrid::RedKgrid
	Gidx::Vector{Vector{Int}}
	Glist::Vector{ReducedCoordinates{Int}}
	mapk::Vector{NTuple{3, Int}}
	mapG::Vector{NTuple{3, Int}}
	kshift_interpolation::NTuple{3, Float64}
	WT::Type{WaveT}
	extrapolation_up::Vector{etpT}
	extrapolation_dn::Vector{etpT}
end
function (etp::WannierWaveKgridShiftInterpolationSpinor{WaveT, etpT})(q::ReducedCoordinates{Float64}) where {WaveT, etpT}
	(nkx, nky, nkz) = etp.kgrid.kgrid_size
	(kx_shift, ky_shift, kz_shift) = etp.kshift_interpolation
	mapkq = Vector{NTuple{3, Float64}}(undef, etp.nk)
	for ik in Base.OneTo(etp.nk)
		kq = etp.kgrid[ik] + q
		kqx_idx = (kq[1] + kx_shift) * nkx + 1
		kqy_idx = (kq[2] + ky_shift) * nky + 1
		kqz_idx = (kq[3] + kz_shift) * nkz + 1
		mapkq[ik] = (kqx_idx, kqy_idx, kqz_idx)
	end
	wavekq = [Matrix{WaveT}(undef, etp.nG[ik] * 2, etp.nw) for ik in 1:etp.nk]
	Threads.@threads for I in CartesianIndices((etp.nw, etp.nk))
		(iw, ik) = Tuple(I)
		(kqx_idx, kqy_idx, kqz_idx) = mapkq[ik]
		wavekq_iw_ik_up = view(wavekq[ik], 1:etp.nG[ik], iw)
		wavekq_iw_ik_dn = view(wavekq[ik], (etp.nG[ik]+1):(2*etp.nG[ik]), iw)
		for (iG, iG_all) in enumerate(etp.Gidx[ik])
			(Gx_idx, Gy_idx, Gz_idx) = etp.mapG[iG_all]
			wavekq_iw_ik_up[iG] = etp.extrapolation_up[iw](
				Gx_idx + kqx_idx,
				Gy_idx + kqy_idx,
				Gz_idx + kqz_idx,
			)
			wavekq_iw_ik_dn[iG] = etp.extrapolation_dn[iw](
				Gx_idx + kqx_idx,
				Gy_idx + kqy_idx,
				Gz_idx + kqz_idx,
			)
		end
		wavekq_iw_ik = view(wavekq[ik], :, iw)
		wave_norm = norm(wavekq_iw_ik)
		if !isapprox(wave_norm, 1; atol = 1e-6)
			wavekq_iw_ik ./= wave_norm
		end
	end
	return wavekq
end
function (etp::WannierWaveKgridShiftInterpolationSpinor{WaveT, etpT})(wavekq::AbstractVector{<:AbstractMatrix{WaveT}}, q::ReducedCoordinates{Float64}) where {WaveT, etpT}
	(nkx, nky, nkz) = etp.kgrid.kgrid_size
	(kx_shift, ky_shift, kz_shift) = etp.kshift_interpolation
	mapkq = Vector{NTuple{3, Float64}}(undef, etp.nk)
	for ik in Base.OneTo(etp.nk)
		kq = etp.kgrid[ik] + q
		kqx_idx = (kq[1] + kx_shift) * nkx + 1
		kqy_idx = (kq[2] + ky_shift) * nky + 1
		kqz_idx = (kq[3] + kz_shift) * nkz + 1
		mapkq[ik] = (kqx_idx, kqy_idx, kqz_idx)
	end
	Threads.@threads for I in CartesianIndices((etp.nw, etp.nk))
		(iw, ik) = Tuple(I)
		(kqx_idx, kqy_idx, kqz_idx) = mapkq[ik]
		wavekq_iw_ik_up = view(wavekq[ik], 1:etp.nG[ik], iw)
		wavekq_iw_ik_dn = view(wavekq[ik], (etp.nG[ik]+1):(2*etp.nG[ik]), iw)
		for (iG, iG_all) in enumerate(etp.Gidx[ik])
			(Gx_idx, Gy_idx, Gz_idx) = etp.mapG[iG_all]
			wavekq_iw_ik_up[iG] = etp.extrapolation_up[iw](
				Gx_idx + kqx_idx,
				Gy_idx + kqy_idx,
				Gz_idx + kqz_idx,
			)
			wavekq_iw_ik_dn[iG] = etp.extrapolation_dn[iw](
				Gx_idx + kqx_idx,
				Gy_idx + kqy_idx,
				Gz_idx + kqz_idx,
			)
		end
		wavekq_iw_ik = view(wavekq[ik], :, iw)
		wave_norm = norm(wavekq_iw_ik)
		if !isapprox(wave_norm, 1; atol = 1e-6)
			wavekq_iw_ik ./= wave_norm
		end
	end
	return wavekq
end
function WannierWaveKgridShiftInterpolationSpinor(wave::WannierPWKgrid{WaveT}) where {WaveT}

	nk = wave.nk
	nw = wave.nw
	nG = wave.nG
	kgrid = wave.kgrid
	Gidx = wave.Gidx
	Glist = wave.Glist

	(minGx, maxGx) = extrema(G -> G[1], Glist)
	(minGy, maxGy) = extrema(G -> G[2], Glist)
	(minGz, maxGz) = extrema(G -> G[3], Glist)
	nGx = maxGx - minGx + 1
	nGy = maxGy - minGy + 1
	nGz = maxGz - minGz + 1
	(nkx, nky, nkz) = kgrid.kgrid_size
	if iszero(kgrid.kshift[1])
		kx_shift = isodd(nkx) ? floor(Int, nkx / 2) // nkx : 1 // 2
	else
		kx_shift = isodd(nkx) ? 1 // 2 : (nkx - 1) // (2 * nkx)
	end
	if iszero(kgrid.kshift[2])
		ky_shift = isodd(nky) ? floor(Int, nky / 2) // nky : 1 // 2
	else
		ky_shift = isodd(nky) ? 1 // 2 : (nky - 1) // (2 * nky)
	end
	if iszero(kgrid.kshift[3])
		kz_shift = isodd(nkz) ? floor(Int, nkz / 2) // nkz : 1 // 2
	else
		kz_shift = isodd(nkz) ? 1 // 2 : (nkz - 1) // (2 * nkz)
	end

	mapk = Vector{NTuple{3, Int}}(undef, nk)
	for ik in Base.OneTo(nk)
		k = kgrid[ik]
		kx_idx = Int((k[1] + kx_shift) * nkx) + 1
		ky_idx = Int((k[2] + ky_shift) * nky) + 1
		kz_idx = Int((k[3] + kz_shift) * nkz) + 1
		mapk[ik] = (kx_idx, ky_idx, kz_idx)
	end
	nGlist = length(Glist)
	mapG = Vector{NTuple{3, Int}}(undef, nGlist)
	for iG in Base.OneTo(nGlist)
		G = Glist[iG]
		Gx_idx = (G[1] - minGx) * nkx
		Gy_idx = (G[2] - minGy) * nky
		Gz_idx = (G[3] - minGz) * nkz
		mapG[iG] = (Gx_idx, Gy_idx, Gz_idx)
	end

	nkGx = nGx * nkx
	nkGy = nGy * nky
	nkGz = nGz * nkz
	C = Array{WT}(undef, nkGx, nkGy, nkGz)

	etps_up = map(1:nw) do iw
		fill!(C, zero(WaveT))
		Threads.@threads for ik in Base.OneTo(nk)
			(kx_idx, ky_idx, kz_idx) = mapk[ik]
			wave_ik = view(wave.components[ik], 1:nG[ik], iw)
			for (iG, iG_all) in enumerate(Gidx[ik])
				(Gx_idx, Gy_idx, Gz_idx) = mapG[iG_all]
				C[
					Gx_idx+kx_idx,
					Gy_idx+ky_idx,
					Gz_idx+kz_idx,
				] = wave_ik[iG]
			end
		end
		C_itp = interpolate(C, BSpline(Cubic(Flat(OnGrid()))))
		C_etp = extrapolate(C_itp, Flat())
		return C_etp
	end
	etps_dn = map(1:nw) do iw
		fill!(C, zero(WaveT))
		Threads.@threads for ik in Base.OneTo(nk)
			(kx_idx, ky_idx, kz_idx) = mapk[ik]
			wave_ik = view(wave.components[ik], (nG[ik]+1):(2*nG[ik]), iw)
			for (iG, iG_all) in enumerate(Gidx[ik])
				(Gx_idx, Gy_idx, Gz_idx) = mapG[iG_all]
				C[
					Gx_idx+kx_idx,
					Gy_idx+ky_idx,
					Gz_idx+kz_idx,
				] = wave_ik[iG]
			end
		end
		C_itp = interpolate(C, BSpline(Cubic(Flat(OnGrid()))))
		C_etp = extrapolate(C_itp, Flat())
		return C_etp
	end

	return WannierWaveKgridShiftInterpolationSpinor(
		nk, nw, nG, kgrid, Gidx, Glist, mapk, mapG,
		Float64.((kx_shift, ky_shift, kz_shift)), WaveT, etps_up, etps_dn,
	)
end

function wavek_normalize!(Glist::Vector{Vector{ReducedCoordinates{Int}}}, kgrid::RedKgrid)
	Threads.@threads for ik in eachindex(Glist)
		k = kgrid[ik]
		k_regular = normalize_kdirect(k)
		ΔG = ReducedCoordinates{Int}(k - k_regular)
		kgrid.kdirect[ik] = k_regular
		Glist[ik] .+= [ΔG]
	end
	return nothing
end
