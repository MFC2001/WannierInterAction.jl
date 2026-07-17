struct RhoR_onlyU{RhoT <: Number, RhokT <: AbstractRhok{RhoT}} <: AbstractRhoR{RhoT}
	nk::Int
	nw::Int
	nq::Int
	kgrid::RedKgrid
	qgrid::RedKgrid
	nG::Vector{Int}
	nG_max::Int
	have_negq0::Bool
	nG_negq0::Int
	fft_grid::NTuple{3, Int}
	ρk::RhokT
	ik_chunks::Vector{Tuple{Int, UnitRange{Int}}}
	# buffer
	nchunk::Int
	idx_chunks::Vector{Vector{Tuple{Int, UnitRange{Int}, Int}}} # ik_start, ik_else, i
	ρR_buffer::Vector{Vector{RhoT}}
	ρk_buffer::Vector{Vector{RhoT}}
	wavek_buffer::Vector{Array{RhoT, 3}}
	wavekq_buffer::Vector{Array{RhoT, 3}}
	locks::Vector{ReentrantLock}
end
function (ρR::RhoR_onlyU{RhoT, RhokT})(q::Int) where {RhoT, RhokT}
	if q == -1
		ρR.have_negq0 || error("No negative q0 data!")
		nG_q = ρR.nG_negq0
	else
		nG_q = ρR.nG[q]
	end
	ρR_data = Matrix{RhoT}(undef, nG_q, ρR.nw)
	return ρR(ρR_data, q)
end
function (ρR::RhoR_onlyU{RhoT, RhokT})(ρR_data::AbstractArray{RhoT, 2}, q::Int) where {RhoT, RhokT}
	if q == -1
		ρR.have_negq0 || error("No negative q0 data!")
		nG_q = ρR.nG_negq0
	else
		nG_q = ρR.nG[q]
	end
	ρR_view = view(ρR_data, 1:nG_q, 1:ρR.nw)
	fill!(ρR_view, zero(RhoT))
	Threads.@threads for ichunk in Base.OneTo(ρR.nchunk)
		@inbounds begin
			ρR_ichunk = view(ρR.ρR_buffer[ichunk], 1:nG_q)
			ρk_ichunk = view(ρR.ρk_buffer[ichunk], 1:nG_q)
			wavek_ichunk = ρR.wavek_buffer[ichunk]
			wavekq_ichunk = ρR.wavekq_buffer[ichunk]
			for (ik_first, ik_else, i) in ρR.idx_chunks[ichunk]
				ρR.ρk(ρk_ichunk, wavek_ichunk, wavekq_ichunk, ik_first, i, i, q)
				ρR_ichunk .= ρk_ichunk
				for ik in ik_else
					ρR.ρk(ρk_ichunk, wavek_ichunk, wavekq_ichunk, ik, i, i, q)
					ρR_ichunk .+= ρk_ichunk
				end
				ρR_i = view(ρR_view, :, i)
				lock(ρR.locks[i])
				try
					ρR_i .+= ρR_ichunk
				finally
					unlock(ρR.locks[i])
				end
			end
		end
	end
	return ρR_data
end
function RhoR_onlyU(ρk::AbstractRhok{RhoT}; nchunk::Integer = Threads.nthreads()) where {RhoT}
	nk = ρk.nk
	nw = ρk.nw

	ik_chunks = _RhoR_split_tasks(nk, nchunk)
	idx_chunks = _RhoR_split_tasks_wistart(nk, nw, nchunk)

	nG_max = ρk.nG_max
	fft_grid = ρk.fft_grid

	# ρR buffer
	ρR_size = format_memory(sizeof(RhoT) * nG_max * nchunk)
	println("ρR: $RhoT, len=$nG_max, nchunk=$nchunk, total=$ρR_size")
	# ρk buffer
	ρk_size = format_memory(sizeof(RhoT) * nG_max * nchunk)
	println("ρk: $RhoT, len=$nG_max, nchunk=$nchunk, total=$ρk_size")
	# Wave functions (ψk & ψkq)
	wave_size = format_memory(sizeof(RhoT) * prod(fft_grid))
	total_size = format_memory(sizeof(RhoT) * prod(fft_grid) * nchunk * 2)
	println("ψ_{k} and ψ_{k+q}: $RhoT, FFT_grid=$fft_grid, single=$wave_size, nchunk=$nchunk, total=$total_size")

	flush(stdout)

	ρR_buffer = [Vector{RhoT}(undef, nG_max) for _ in 1:nchunk]
	ρk_buffer = [Vector{RhoT}(undef, nG_max) for _ in 1:nchunk]
	wavek_buffer = [Array{RhoT, 3}(undef, fft_grid) for _ in 1:nchunk]
	wavekq_buffer = [Array{RhoT, 3}(undef, fft_grid) for _ in 1:nchunk]

	locks = [ReentrantLock() for _ in 1:nw]

	return RhoR_onlyU{RhoT, typeof(ρk)}(
		nk, nw, ρk.nq,
		ρk.kgrid, ρk.qgrid,
		ρk.nG, nG_max,
		ρk.have_negq0, ρk.nG_negq0,
		fft_grid, ρk, ik_chunks,
		nchunk, idx_chunks,
		ρR_buffer, ρk_buffer, wavek_buffer, wavekq_buffer,
		locks,
	)
end
