function realspace(wave::WannierPWKgrid{WT}; kwargs...) where {WT <: Number}
	if wave.nspin == 4
		return realspace_spinor(wave; kwargs...)
	else
		return realspace_scalar(wave; kwargs...)
	end
end
function realspace_scalar(wave::WannierPWKgrid{WT}; range = :default, fft_grid = wave.fft_grid,
	cell::Union{Cell, Nothing} = nothing) where {WT <: Number}
	@assert wave.nspin ≠ 4 "Only work for scalar wave function!"

	nk = wave.nk
	nw = wave.nw
	kgrid = wave.kgrid

	if range == :default
		range = [1, 1, 1]
		range[wave.period] .= 3
		range = (range[1], range[2], range[3])
	end
	aimR = gridindex(range)
	nR = length(aimR)

	fft_grid = (Int(fft_grid[1]), Int(fft_grid[2]), Int(fft_grid[3]))
	if fft_grid == wave.fft_grid
		Gmap3D = wave.Gmap
	else
		pgi = PeriodicGridIndex(fft_grid)
		Gmap3D = pgi.(Val(:linear), wave.Glist)
	end

	ift_wave = plan_bfft!(Array{WT, 3}(undef, fft_grid); flags = FFTW.ESTIMATE)
	real_grid = [ReducedCoordinates((ix - 1) // fft_grid[1], (iy - 1) // fft_grid[2], (iz - 1) // fft_grid[3])
				 for ix in 1:fft_grid[1], iy in 1:fft_grid[2], iz in 1:fft_grid[3]]

	wannier_real_grid = Vector{Array{WT, 4}}(undef, nw)
	for iw in 1:nw
		wannier_real_grid[iw] = zeros(WT, fft_grid..., nR)
	end
	lk = ReentrantLock()
	Threads.@threads for ik in Base.OneTo(nk)
		wave3D = Array{WT, 3}(undef, fft_grid...)
		wave3D_R = Array{WT, 3}(undef, fft_grid...)
		Gmap3D_ik = Gmap3D[wave.Gidx[ik]]
		k = kgrid[ik]
		ekr = map(r -> cispi(2 * (k ⋅ r)), real_grid)
		ekR = [cispi(2 * (k ⋅ R)) for R in aimR]
		for iw in 1:nw
			wave_iw = view(wave.components[ik], :, iw)
			fill!(wave3D, zero(WT))
			for (Gmap, value) in zip(Gmap3D_ik, wave_iw)
				wave3D[Gmap] = value
			end
			ift_wave * wave3D
			wave3D .*= ekr
			for iR in eachindex(aimR)
				wave3D_R .= wave3D .* ekR[iR]
				lock(lk) do
					wannier_real_grid[iw][:, :, :, iR] .+= wave3D_R
				end
			end
		end
	end

	for iw in 1:nw
		wannier_real_grid_iw = wannier_real_grid[iw]
		(norm2_max, idx) = findmax(abs2, wannier_real_grid_iw)
		negative_phase = conj(wannier_real_grid_iw[idx]) / sqrt(norm2_max)
		negative_phase_inv_nk = negative_phase / nk
		wannier_real_grid_iw .*= negative_phase_inv_nk
	end

	(minRx, maxRx) = extrema(R -> R[1], aimR)
	(minRy, maxRy) = extrema(R -> R[2], aimR)
	(minRz, maxRz) = extrema(R -> R[3], aimR)
	nRx = maxRx - minRx + 1
	nRy = maxRy - minRy + 1
	nRz = maxRz - minRz + 1
	R_map = Vector{NTuple{3, Int}}(undef, length(aimR))
	for (iR, R) in enumerate(aimR)
		iRx = R[1] - minRx + 1
		iRy = R[2] - minRy + 1
		iRz = R[3] - minRz + 1
		R_map[iR] = (iRx, iRy, iRz)
	end

	wannier_real_grid_reshape = Vector{Array{WT, 3}}(undef, nw)
	for iw in Base.OneTo(nw)
		wannier_real_grid_reshape_iw = Array{WT, 3}(undef, fft_grid[1] * nRx, fft_grid[2] * nRy, fft_grid[3] * nRz)
		for iR in eachindex(aimR)
			(iRx, iRy, iRz) = R_map[iR]
			wannier_real_grid_reshape_iw[
				(fft_grid[1]*(iRx-1)+1):(fft_grid[1]*iRx),
				(fft_grid[2]*(iRy-1)+1):(fft_grid[2]*iRy),
				(fft_grid[3]*(iRz-1)+1):(fft_grid[3]*iRz),
			] .= wannier_real_grid[iw][:, :, :, iR]
		end
		wannier_real_grid_reshape[iw] = wannier_real_grid_reshape_iw
	end

	if isnothing(cell)
		return wannier_real_grid_reshape
	else
		return map(iw -> wannier90_xsf(wannier_real_grid_reshape[iw], range, cell), Base.OneTo(nw))
	end
end
function realspace_spinor(wave::WannierPWKgrid{WT}; range = :default, fft_grid = wave.fft_grid,
	cell::Union{Cell, Nothing} = nothing) where {WT <: Number}
	@assert wave.nspin == 4 "Only work for spinor wave function!"

	nk = wave.nk
	nw = wave.nw
	kgrid = wave.kgrid

	if range == :default
		range = [1, 1, 1]
		range[wave.period] .= 3
		range = (range[1], range[2], range[3])
	end
	aimR = gridindex(range)
	nR = length(aimR)

	fft_grid = (Int(fft_grid[1]), Int(fft_grid[2]), Int(fft_grid[3]))
	if fft_grid == wave.fft_grid
		Gmap3D = wave.Gmap
	else
		pgi = PeriodicGridIndex(fft_grid)
		Gmap3D = pgi.(Val(:linear), wave.Glist)
	end

	ift_wave = plan_bfft!(Array{WT, 3}(undef, fft_grid); flags = FFTW.ESTIMATE)
	real_grid = [ReducedCoordinates((ix - 1) // fft_grid[1], (iy - 1) // fft_grid[2], (iz - 1) // fft_grid[3])
				 for ix in 1:fft_grid[1], iy in 1:fft_grid[2], iz in 1:fft_grid[3]]

	wannier_real_grid = Matrix{Array{WT, 4}}(undef, 2, nw)
	for iw in 1:nw, ispin in 1:2
		wannier_real_grid[ispin, iw] = zeros(WT, fft_grid..., nR)
	end
	lk = ReentrantLock()
	Threads.@threads for ik in Base.OneTo(nk)
		wave3D = Array{WT, 3}(undef, fft_grid...)
		wave3D_R = Array{WT, 3}(undef, fft_grid...)
		nG_ik = wave.nG[ik]
		Gidx = [1:nG_ik, (nG_ik+1):(2*nG_ik)]
		Gmap3D_ik = Gmap3D[wave.Gidx[ik]]
		k = kgrid[ik]
		ekr = map(r -> cispi(2 * (k ⋅ r)), real_grid)
		ekR = [cispi(2 * (k ⋅ R)) for R in aimR]
		for iw in 1:nw, ispin in 1:2
			wave_iw = view(wave.components[ik], Gidx[ispin], iw)
			fill!(wave3D, zero(WT))
			for (Gmap, value) in zip(Gmap3D_ik, wave_iw)
				wave3D[Gmap] = value
			end
			ift_wave * wave3D
			wave3D .*= ekr
			for iR in eachindex(aimR)
				wave3D_R .= wave3D .* ekR[iR]
				lock(lk) do
					wannier_real_grid[ispin, iw][:, :, :, iR] .+= wave3D_R
				end
			end
		end
	end

	for iw in 1:nw
		wannier_real_grid_iw_up = wannier_real_grid[1, iw]
		wannier_real_grid_iw_dn = wannier_real_grid[2, iw]
		(norm2_max_up, idx_up) = findmax(abs2, wannier_real_grid_iw_up)
		(norm2_max_dn, idx_dn) = findmax(abs2, wannier_real_grid_iw_dn)
		if norm2_max_up > norm2_max_dn
			negative_phase = conj(wannier_real_grid_iw_up[idx_up]) / sqrt(norm2_max_up)
		else
			negative_phase = conj(wannier_real_grid_iw_dn[idx_dn]) / sqrt(norm2_max_dn)
		end
		negative_phase_inv_nk = negative_phase / nk
		wannier_real_grid_iw_up .*= negative_phase_inv_nk
		wannier_real_grid_iw_dn .*= negative_phase_inv_nk
	end

	(minRx, maxRx) = extrema(R -> R[1], aimR)
	(minRy, maxRy) = extrema(R -> R[2], aimR)
	(minRz, maxRz) = extrema(R -> R[3], aimR)
	nRx = maxRx - minRx + 1
	nRy = maxRy - minRy + 1
	nRz = maxRz - minRz + 1
	R_map = Vector{NTuple{3, Int}}(undef, length(aimR))
	for (iR, R) in enumerate(aimR)
		iRx = R[1] - minRx + 1
		iRy = R[2] - minRy + 1
		iRz = R[3] - minRz + 1
		R_map[iR] = (iRx, iRy, iRz)
	end

	wannier_real_grid_reshape = Matrix{Array{WT, 3}}(undef, 2, nw)
	for iw in Base.OneTo(nw), ispin in 1:2
		wannier_real_grid_reshape_iw = Array{WT, 3}(undef, fft_grid[1] * nRx, fft_grid[2] * nRy, fft_grid[3] * nRz)
		for iR in eachindex(aimR)
			(iRx, iRy, iRz) = R_map[iR]
			wannier_real_grid_reshape_iw[
				(fft_grid[1]*(iRx-1)+1):(fft_grid[1]*iRx),
				(fft_grid[2]*(iRy-1)+1):(fft_grid[2]*iRy),
				(fft_grid[3]*(iRz-1)+1):(fft_grid[3]*iRz),
			] .= wannier_real_grid[ispin, iw][:, :, :, iR]
		end
		wannier_real_grid_reshape[ispin, iw] = wannier_real_grid_reshape_iw
	end

	if isnothing(cell)
		return wannier_real_grid_reshape
	else
		return [wannier90_xsf(wannier_real_grid_reshape[ispin, iw], range, cell) for ispin in 1:2, iw in 1:nw]
	end
end
