"""
    realspace(wave::WannierPWKgrid{WT}, Rlist::AbstractVector; fft_grid = wave.fft_grid)

Compute the real-space representation of the Wannier functions contained in `wave` via an inverse FFT, 
evaluated at the lattice cells specified by `Rlist`.

# Arguments
- `wave::WannierPWKgrid{WT}`: The object containing the Wannier functions to be transformed into real space.
- `Rlist::AbstractVector`: The lattice cells over which the Wannier functions are evaluated.
  Each element must be a three-integer vector `[n1, n2, n3]` specifying the cell reached by the corresponding lattice vector.

# Keywords
- `fft_grid = wave.fft_grid`: The grid density, i.e. the number of grid points per unit cell used for the FFT.

# Returns
- For scalar Wannier functions: a vector of length `nw` (the number of Wannier functions).
  Each element is an array of size `[nfx, nfy, nfz, nR]`, where the last dimension is indexed consistently with `Rlist`.
- For spinor Wannier functions: a `[2, nw]` matrix, where the first dimension is the spin degree of freedom. 
  Each element is an array of size `[nfx, nfy, nfz, nR]`, with the last dimension indexed consistently with `Rlist`.
"""
function WannierInterActionBase.realspace(wave::WannierPWKgrid{WT}, Rlist::AbstractVector; kwargs...) where {WT<:Number}
    Rlist = collect(ReducedCoordinates{Int}, Rlist)
    if wave.nspinor == 2
        return realspace_spinor(wave, Rlist; kwargs...)
    else
        return realspace_scalar(wave, Rlist; kwargs...)
    end
end
function realspace_scalar(wave::WannierPWKgrid{WT}, Rlist::Vector{ReducedCoordinates{Int}};
    fft_grid=wave.fft_grid) where {WT<:Number}

    nk = wave.nk
    nw = wave.nw
    kgrid = wave.kgrid
    nR = length(Rlist)

    fft_grid = (Int(fft_grid[1]), Int(fft_grid[2]), Int(fft_grid[3]))
    if fft_grid == wave.fft_grid
        Gmap3D = wave.Gmap
    else
        pgi = PeriodicGridIndex(fft_grid)
        Gmap3D = pgi.(Val(:linear), wave.Glist)
    end

    ift_wave = plan_bfft!(Array{WT,3}(undef, fft_grid); flags=FFTW.ESTIMATE)
    real_grid = [ReducedCoordinates((ix - 1) // fft_grid[1], (iy - 1) // fft_grid[2], (iz - 1) // fft_grid[3])
                 for ix in 1:fft_grid[1], iy in 1:fft_grid[2], iz in 1:fft_grid[3]]

    wannier_real_grid = Vector{Array{WT,4}}(undef, nw)
    for iw in 1:nw
        wannier_real_grid[iw] = zeros(WT, fft_grid..., nR)
    end
    lks = [ReentrantLock() for _ in 1:nR, _ in 1:nw]
    Threads.@threads for ik in Base.OneTo(nk)
        wave3D = Array{WT,3}(undef, fft_grid...)
        wave3D_R = Array{WT,3}(undef, fft_grid...)
        Gmap3D_ik = Gmap3D[wave.Gidx[ik]]
        k = kgrid[ik]
        ekr = map(r -> cispi(2 * (k ⋅ r)), real_grid)
        ekR = [cispi(2 * (k ⋅ R)) for R in Rlist]
        for iw in 1:nw
            wave_iw = view(wave.components[ik], :, iw)
            fill!(wave3D, zero(WT))
            for (Gmap, value) in zip(Gmap3D_ik, wave_iw)
                wave3D[Gmap] = value
            end
            ift_wave * wave3D
            wave3D .*= ekr
            for iR in eachindex(Rlist)
                wave3D_R .= wave3D .* ekR[iR]
                lock(lks[iR, iw]) do
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

    return wannier_real_grid
end
function realspace_spinor(wave::WannierPWKgrid{WT}, Rlist::Vector{ReducedCoordinates{Int}};
    fft_grid=wave.fft_grid) where {WT<:Number}

    nk = wave.nk
    nw = wave.nw
    kgrid = wave.kgrid
    nR = length(Rlist)

    fft_grid = (Int(fft_grid[1]), Int(fft_grid[2]), Int(fft_grid[3]))
    if fft_grid == wave.fft_grid
        Gmap3D = wave.Gmap
    else
        pgi = PeriodicGridIndex(fft_grid)
        Gmap3D = pgi.(Val(:linear), wave.Glist)
    end

    ift_wave = plan_bfft!(Array{WT,3}(undef, fft_grid); flags=FFTW.ESTIMATE)
    real_grid = [ReducedCoordinates((ix - 1) // fft_grid[1], (iy - 1) // fft_grid[2], (iz - 1) // fft_grid[3])
                 for ix in 1:fft_grid[1], iy in 1:fft_grid[2], iz in 1:fft_grid[3]]

    wannier_real_grid = Matrix{Array{WT,4}}(undef, 2, nw)
    for iw in 1:nw, ispin in 1:2
        wannier_real_grid[ispin, iw] = zeros(WT, fft_grid..., nR)
    end
    lks = [ReentrantLock() for _ in 1:nR, _ in 1:2, _ in 1:nw]
    Threads.@threads for ik in Base.OneTo(nk)
        wave3D = Array{WT,3}(undef, fft_grid...)
        wave3D_R = Array{WT,3}(undef, fft_grid...)
        nG_ik = wave.nG[ik]
        Gidx = [1:nG_ik, (nG_ik+1):(2*nG_ik)]
        Gmap3D_ik = Gmap3D[wave.Gidx[ik]]
        k = kgrid[ik]
        ekr = map(r -> cispi(2 * (k ⋅ r)), real_grid)
        ekR = [cispi(2 * (k ⋅ R)) for R in Rlist]
        for iw in 1:nw, ispin in 1:2
            wave_iw = view(wave.components[ik], Gidx[ispin], iw)
            fill!(wave3D, zero(WT))
            for (Gmap, value) in zip(Gmap3D_ik, wave_iw)
                wave3D[Gmap] = value
            end
            ift_wave * wave3D
            wave3D .*= ekr
            for iR in eachindex(Rlist)
                wave3D_R .= wave3D .* ekR[iR]
                lock(lks[iR, ispin, iw]) do
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

    return wannier_real_grid
end
