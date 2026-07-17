function WannierInterActionBase.PlaneWaveCoulomb(::Val{:bulk}, ρk::AbstractRhok; kwargs...)
	@assert count(ρk.period) == 3 "Wrong period, bulk is used in 3D system."
	return PlaneWaveCoulomb(Val(:bulk), ρk.lattice, ρk.Glist; kwargs...)
end
function WannierInterActionBase.PlaneWaveCoulomb(::Val{:slab}, ρk::AbstractRhok; kwargs...)
	@assert count(ρk.period) == 2 "Wrong period, slab truncation is used in 2D system."
	return PlaneWaveCoulomb(Val(:slab), ρk.lattice, ρk.Glist; period = ρk.period, kwargs...)
end
# function WannierInterActionBase.PlaneWaveCoulomb(::Val{:wire}, ρk::AbstractRhok; fft_factor::Integer = 2, kwargs...)
# 	@assert count(ρk.period) == 1 "Wrong period, wire truncation is used in 1D system."
# 	fft_grid = ρk.fft_grid .* fft_factor
# 	return PlaneWaveCoulomb(Val(:wire), ρk.lattice, ρk.Glist; period = ρk.period, fft_grid, kwargs...)
# end
# function WannierInterActionBase.PlaneWaveCoulomb(::Val{:cylinder}, ρk::AbstractRhok; kwargs...)
# 	@assert count(ρk.period) == 1 "Wrong period, cylinder truncation is used in 1D system."
# 	# 还有关键词参数柱半径R.
# 	return PlaneWaveCoulomb(Val(:cylinder), ρk.lattice, ρk.Glist; period = ρk.period, kwargs...)
# end
# function WannierInterActionBase.PlaneWaveCoulomb(::Val{:box}, ρk::AbstractRhok; fft_factor::Integer = 1, kwargs...)
# 	@assert count(ρk.period) == 0 "Wrong period, box truncation is used in 0D system."
# 	# REDFT默认在一半区间上变换，所以fft_factor=1便代表每个方向都将fft_grid扩大二倍。
# 	fft_grid = ρk.fft_grid .* fft_factor
# 	# 还有关键词参数fft_nt。
# 	return PlaneWaveCoulomb(Val(:box), ρk.lattice, ρk.Glist; fft_grid, kwargs...)
# end
# function WannierInterActionBase.PlaneWaveCoulomb(::Val{:sphere}, ρk::AbstractRhok; kwargs...)
# 	@assert count(ρk.period) == 0 "Wrong period, sphere truncation is used in 0D system."
# 	# 还有关键词参数球半径R.
# 	return PlaneWaveCoulomb(Val(:sphere), ρk.lattice, ρk.Glist; kwargs...)
# end
