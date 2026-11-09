"""
    realspace(wave::WannierPWKgrid{WT}; range = :default, fft_grid = wave.fft_grid, cell::Union{Cell, Nothing} = nothing)

Compute the real-space representation of the Wannier functions contained in `wave` via an inverse FFT.

# Arguments
- `wave::WannierPWKgrid{WT}`: The object containing the Wannier functions to be transformed into real space.

# Keywords
- `range = :default`: The real-space range over which the Wannier functions are evaluated.
  A tuple `(n1, n2, n3)` corresponds to the ranges along the three directions; 
  for example, `(3, 3, 1)` maps to `(-1:1, -1:1, 0)`. With the default
  `:default`, the range follows the periodicity specified by `wave.period`: periodic
  directions use `-1:1`, while non-periodic directions are fixed to `0`.
- `fft_grid = wave.fft_grid`: The grid density, i.e. the number of grid points per unit cell used for the FFT.
- `cell::Union{Cell, Nothing} = nothing`: 
  If provided (as a [`Cell`](@ref)), the function returns a vector of [`wannier90_xsf`](@ref) objects; 
  otherwise, it returns a vector of arrays.

# Returns
A vector whose length equals the number of Wannier functions. 
Each element is an array, or a [`wannier90_xsf`](@ref) object when `cell` is provided. 
If `wave` is a spinor, the result is instead a matrix whose first dimension has size 2 (the spin degree of freedom).
"""
function WannierInterActionBase.realspace(wave::WannierPWKgrid{WT}; kwargs...) where {WT<:Number}
    if wave.nspinor == 2
        return realspace_spinor(wave; kwargs...)
    else
        return realspace_scalar(wave; kwargs...)
    end
end
function realspace_scalar(wave::WannierPWKgrid{WT}; range=:default, fft_grid=wave.fft_grid,
    cell::Union{Cell,Nothing}=nothing) where {WT<:Number}

    if range == :default
        range = [1, 1, 1]
        range[wave.period] .= 3
        range = (range[1], range[2], range[3])
    end
    Rlist = gridindex(range)
    wannier_real_grid = realspace_scalar(wave, Rlist; fft_grid)

    (minRx, maxRx) = extrema(R -> R[1], Rlist)
    (minRy, maxRy) = extrema(R -> R[2], Rlist)
    (minRz, maxRz) = extrema(R -> R[3], Rlist)
    nRx = maxRx - minRx + 1
    nRy = maxRy - minRy + 1
    nRz = maxRz - minRz + 1
    R_map = Vector{NTuple{3,Int}}(undef, length(Rlist))
    for (iR, R) in enumerate(Rlist)
        iRx = R[1] - minRx + 1
        iRy = R[2] - minRy + 1
        iRz = R[3] - minRz + 1
        R_map[iR] = (iRx, iRy, iRz)
    end

    nw = wave.nw
    wannier_real_grid_reshape = Vector{Array{WT,3}}(undef, nw)
    for iw in Base.OneTo(nw)
        wannier_real_grid_reshape_iw = Array{WT,3}(undef, fft_grid[1] * nRx, fft_grid[2] * nRy, fft_grid[3] * nRz)
        for iR in eachindex(Rlist)
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
function realspace_spinor(wave::WannierPWKgrid{WT}; range=:default, fft_grid=wave.fft_grid,
    cell::Union{Cell,Nothing}=nothing) where {WT<:Number}

    if range == :default
        range = [1, 1, 1]
        range[wave.period] .= 3
        range = (range[1], range[2], range[3])
    end
    Rlist = gridindex(range)
    wannier_real_grid = realspace_spinor(wave, Rlist; fft_grid)

    (minRx, maxRx) = extrema(R -> R[1], Rlist)
    (minRy, maxRy) = extrema(R -> R[2], Rlist)
    (minRz, maxRz) = extrema(R -> R[3], Rlist)
    nRx = maxRx - minRx + 1
    nRy = maxRy - minRy + 1
    nRz = maxRz - minRz + 1
    R_map = Vector{NTuple{3,Int}}(undef, length(Rlist))
    for (iR, R) in enumerate(Rlist)
        iRx = R[1] - minRx + 1
        iRy = R[2] - minRy + 1
        iRz = R[3] - minRz + 1
        R_map[iR] = (iRx, iRy, iRz)
    end

    nw = wave.nw
    wannier_real_grid_reshape = Matrix{Array{WT,3}}(undef, 2, nw)
    for iw in Base.OneTo(nw), ispin in 1:2
        wannier_real_grid_reshape_iw = Array{WT,3}(undef, fft_grid[1] * nRx, fft_grid[2] * nRy, fft_grid[3] * nRz)
        for iR in eachindex(Rlist)
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
