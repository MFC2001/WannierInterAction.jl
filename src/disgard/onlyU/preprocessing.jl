function _preprocessing_int_onlyU(ρk::AbstractRhok{RhoT}, nω, nchunk) where {RhoT}

	nw = ρk.nw
	qgrid = ρk.qgrid
	lattice = ρk.lattice

	ijji_R = gridindex(Val(:WignerSeitz), lattice, qgrid.kgrid_size)
	ijji_idx = [(iw1, iw2) for iw1 in 1:nw, iw2 in 1:nw]
	ijji_value = iszero(nω) ? zeros(RhoT, nw, nw, length(ijji_R)) : zeros(RhoT, nω + 1, nw, nw, length(ijji_R))
	ijji_value_iq = iszero(nω) ? Array{RhoT}(undef, nw, nw) : Array{RhoT}(undef, nω + 1, nw, nw)

	# _idx
	ijji_idx = _split_tasks(ijji_idx, nchunk)

	return (nchunk = nchunk, ijji_value = ijji_value, ijji_R = ijji_R, ijji_value_iq = ijji_value_iq, ijji_idx = ijji_idx)
end
function _preprocessing_buffer_onlyU(::Val{3}, ρk::AbstractRhok{RhoT}, nchunk) where {RhoT}
	nw = ρk.nw
	nG_rho_max = ρk.nG_max

	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nw)
	println("ρR: $RhoT, dims=($nG_rho_max, $nw), total=$ρR_size")

	flush(stdout)

	ρR = Array{RhoT}(undef, nG_rho_max, nw)

	return (ρR = ρR,)
end
function _preprocessing_buffer_onlyU(::Union{Val{2}, Val{1}}, ρk::AbstractRhok{RhoT}, nchunk) where {RhoT}
	nw = ρk.nw
	nG_rho_max = ρk.nG_max

	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nw)
	println("ρR: $RhoT, dims=($nG_rho_max, $nw), total=$ρR_size")

	buffer_rho_size = format_memory(sizeof(RhoT) * nG_rho_max)
	total_size = format_memory(sizeof(RhoT) * nG_rho_max * nchunk)
	println("buffer_rho: $RhoT, len=$nG_rho_max, single=$buffer_rho_size, nchunk=$nchunk, total=$total_size")

	flush(stdout)

	ρR = Array{RhoT}(undef, nG_rho_max, nw)
	buffer_rho = [Vector{RhoT}(undef, nG_rho_max) for _ in 1:nchunk]

	return (ρR = ρR, rho = buffer_rho)
end
function _preprocessing_buffer_onlyU(::Val{3}, eps, ρk::AbstractRhok{RhoT}, nchunk) where {RhoT}
	nw = ρk.nw
	nG_eps_max = eps.nG_max
	nG_rho_max = ρk.nG_max

	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nw)
	println("ρR: $RhoT, dims=($nG_rho_max, $nw), total=$ρR_size")

	buffer_eps_size = format_memory(sizeof(RhoT) * nG_eps_max)
	total_size = format_memory(sizeof(RhoT) * nG_eps_max * nchunk * 2)
	println("buffer_eps and buffer_eps_mul: $RhoT, len=$nG_eps_max, single=$buffer_eps_size, nchunk=$nchunk, total=$total_size")

	flush(stdout)

	ρR = Array{RhoT}(undef, nG_rho_max, nw)
	buffer_eps_mul = [Vector{RhoT}(undef, nG_eps_max) for _ in 1:nchunk]
	buffer_eps = [Vector{RhoT}(undef, nG_eps_max) for _ in 1:nchunk]

	return (ρR = ρR, eps_mul = buffer_eps_mul, eps = buffer_eps)
end
function _preprocessing_buffer_onlyU(::Union{Val{2}, Val{1}}, eps, ρk::AbstractRhok{RhoT}, nchunk) where {RhoT}
	nw = ρk.nw
	nG_eps_max = eps.nG_max
	nG_rho_max = ρk.nG_max


	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nw)
	println("ρR: $RhoT, dims=($nG_rho_max, $nw), total=$ρR_size")

	buffer_eps_size = format_memory(sizeof(RhoT) * nG_eps_max)
	total_size = format_memory(sizeof(RhoT) * nG_eps_max * nchunk * 2)
	println("buffer_eps and buffer_eps_mul: $RhoT, len=$nG_eps_max, single=$buffer_eps_size, nchunk=$nchunk, total=$total_size")

	buffer_rho_size = format_memory(sizeof(RhoT) * nG_rho_max)
	total_size = format_memory(sizeof(RhoT) * nG_rho_max * nchunk)
	println("buffer_rho: $RhoT, len=$nG_rho_max, single=$buffer_rho_size, nchunk=$nchunk, total=$total_size")

	flush(stdout)

	ρR = Array{RhoT}(undef, nG_rho_max, nw)
	buffer_eps_mul = [Vector{RhoT}(undef, nG_eps_max) for _ in 1:nchunk]
	buffer_eps = [Vector{RhoT}(undef, nG_eps_max) for _ in 1:nchunk]
	buffer_rho = [Vector{RhoT}(undef, nG_rho_max) for _ in 1:nchunk]

	return (ρR = ρR, eps_mul = buffer_eps_mul, eps = buffer_eps, rho = buffer_rho)
end
