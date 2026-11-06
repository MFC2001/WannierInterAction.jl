"""
	PlaneWaveCoulomb(sym::Symbol, ρk::AbstractRhok; kwargs...)

Construct a truncated Coulomb potential in plane-wave basis for a given system defined by `ρk`.

	PlaneWaveCoulomb(:bulk, ρk::AbstractRhok)
	PlaneWaveCoulomb(:slab, ρk::AbstractRhok; α::Real = 0)
	PlaneWaveCoulomb(:wire, ρk::AbstractRhok; fft_factor::Integer = 2, α::Real=0)
	PlaneWaveCoulomb(:cylinder, ρk::AbstractRhok; R::Union{Real,Nothing}=nothing, α::Real=0)
	PlaneWaveCoulomb(:box, ρk::AbstractRhok; fft_factor::Integer = 1, fft_nt::Integer=Threads.nthreads())
	PlaneWaveCoulomb(:sphere, ρk::AbstractRhok; R::Union{Real,Nothing}=nothing)
"""
function WannierInterActionBase.PlaneWaveCoulomb(::Val{:bulk}, ρk::AbstractRhok)
	@assert count(ρk.period) == 3 "Wrong period, bulk is used in 3D system."
	return PlaneWaveCoulomb(Val(:bulk), ρk.lattice, ρk.Glist)
end
function WannierInterActionBase.PlaneWaveCoulomb(::Val{:slab}, ρk::AbstractRhok; α::Real = 0)
	@assert count(ρk.period) == 2 "Wrong period, slab truncation is used in 2D system."
	# 关键词参数α，详见PhysRevB.73.233103。参考BerkeleyGW，默认设置为0。
	return PlaneWaveCoulomb(Val(:slab), ρk.lattice, ρk.Glist; period = ρk.period, α)
end
function WannierInterActionBase.PlaneWaveCoulomb(::Val{:wire}, ρk::AbstractRhok; fft_factor::Integer = 2, α::Real=0)
	@assert count(ρk.period) == 1 "Wrong period, wire truncation is used in 1D system."
	fft_grid = ρk.fft_grid .* fft_factor
	return PlaneWaveCoulomb(Val(:wire), ρk.lattice, ρk.Glist; period = ρk.period, fft_grid, α)
end
function WannierInterActionBase.PlaneWaveCoulomb(::Val{:cylinder}, ρk::AbstractRhok; R::Union{Real,Nothing}=nothing, α::Real=0)
	@assert count(ρk.period) == 1 "Wrong period, cylinder truncation is used in 1D system."
	# 关键词参数柱半径R.
	return PlaneWaveCoulomb(Val(:cylinder), ρk.lattice, ρk.Glist; period = ρk.period, R, α)
end
function WannierInterActionBase.PlaneWaveCoulomb(::Val{:box}, ρk::AbstractRhok; fft_factor::Integer = 1, fft_nt::Integer=Threads.nthreads())
	@assert count(ρk.period) == 0 "Wrong period, box truncation is used in 0D system."
	# REDFT默认在一半区间上变换，所以fft_factor=1便代表每个方向都将fft_grid扩大二倍。
	fft_grid = ρk.fft_grid .* fft_factor
	# 关键词参数fft_nt。
	return PlaneWaveCoulomb(Val(:box), ρk.lattice, ρk.Glist; fft_grid, fft_nt)
end
function WannierInterActionBase.PlaneWaveCoulomb(::Val{:sphere}, ρk::AbstractRhok; R::Union{Real,Nothing}=nothing)
	@assert count(ρk.period) == 0 "Wrong period, sphere truncation is used in 0D system."
	# 关键词参数球半径R.
	return PlaneWaveCoulomb(Val(:sphere), ρk.lattice, ρk.Glist; R)
end
