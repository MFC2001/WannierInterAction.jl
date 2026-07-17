struct RhoR_nG{RhoT <: Number, RhokT <: AbstractRhok{RhoT}} <: AbstractRhoR{RhoT}
	nR::Int
	np_ij::Int
	ρRlist::Vector{ReducedCoordinates{Int}}
	pair_ij::Vector{NTuple{2, Int}}
	ekρR::Matrix{Complex{Float64}}
	ρk::RhokT
	# buffer
	nchunk::Int
	idx_chunks::Vector{Vector{Tuple{UnitRange{Int}, Int}}} # iks, ij
	ρR_buffer::Vector{Vector{RhoT}}
	ρk_buffer::Vector{Vector{RhoT}}
	wavek_buffer::Vector{Array{RhoT, 3}}
	wavekq_buffer::Vector{Array{RhoT, 3}}
	locks::Matrix{ReentrantLock}
end
function (ρR::RhoR_nG{RhoT, RhokT})(iq::Int) where {RhoT, RhokT}
	ρR_data = Vector{RhoT}(undef, ρR.nG[iq], ρR.nR, ρR.np_ij)
	return ρR(ρR_data, iq)
end
function (ρR::RhoR_nG{RhoT, RhokT})(ρR_data::AbstractArray{RhoT, 3}, iq::Int) where {RhoT, RhokT}
	nG_q = ρR.nG[iq]
	ρR_view = view(ρR_data, 1:nG_q, 1:ρR.nR, 1:ρR.np_ij)
	fill!(ρR_view, zero(RhoT))
	Threads.@threads for ichunk in Base.OneTo(ρR.nchunk)
		@inbounds begin
			ρR_ichunk = view(ρR.ρR_buffer[ichunk], 1:nG_q)
			ρk_ichunk = view(ρR.ρk_buffer[ichunk], 1:nG_q)
			wavek_ichunk = ρR.wavek_buffer[ichunk]
			wavekq_ichunk = ρR.wavekq_buffer[ichunk]
			for (iks, ij) in ρR.idx_chunks[ichunk]
				(i, j) = ρR.pair_ij[ij]
				for ik in iks
					ρR.ρk(ρk_ichunk, wavek_ichunk, wavekq_ichunk, ik, i, j, iq)
					for iR in Base.OneTo(ρR.nR)
						ρR_ichunk .= ρk_ichunk .* ρR.ekρR[iR, ik]
						ρR_iR_ij = view(ρR_view, :, iR, ij)
						lock(ρR.locks[iR, ij])
						try
							ρR_iR_ij .+= ρR_ichunk
						finally
							unlock(ρR.locks[iR, ij])
						end
					end
				end
			end
		end
	end
	return ρR_data
end
# 通常来说，最高并行数为nk * np_ij
function RhoR_nG(ρk::AbstractRhok{RhoT},
	ρRlist::Vector{ReducedCoordinates{Int}},
	pair_ij::Vector{NTuple{2, Int}};
	nchunk::Integer = Threads.nthreads(),
) where {RhoT}

	nk = ρk.nk
	np_ij = length(pair_ij)
	nR = length(ρRlist)

	idx_chunks = _RhoR_split_tasks_wostart(nk, np_ij, nchunk)

	ekρR = [cispi(2 * (k ⋅ R)) / nk for R in ρRlist, k in ρk.kgrid]

	nG_rho_max = ρk.nG_max
	fft_grid = ρk.fft_grid

	# ρR buffer
	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nchunk)
	println("ρR: $RhoT, len=$nG_rho_max, nchunk=$nchunk, total=$ρR_size")
	# ρk buffer
	ρk_size = format_memory(sizeof(RhoT) * nG_rho_max * nchunk)
	println("ρk: $RhoT, len=$nG_rho_max, nchunk=$nchunk, total=$ρk_size")
	# Wave functions (ψk & ψkq)
	wave_size = format_memory(sizeof(RhoT) * prod(fft_grid))
	total_size = format_memory(sizeof(RhoT) * prod(fft_grid) * nchunk * 2)
	println("ψ_{k} and ψ_{k+q}: $RhoT, FFT_grid=$fft_grid, single=$wave_size, nchunk=$nchunk, total=$total_size")

	flush(stdout)

	ρR_buffer = [Vector{RhoT}(undef, nG_rho_max) for _ in 1:nchunk]
	ρk_buffer = [Vector{RhoT}(undef, nG_rho_max) for _ in 1:nchunk]
	wavek_buffer = [Array{RhoT, 3}(undef, fft_grid) for _ in 1:nchunk]
	wavekq_buffer = [Array{RhoT, 3}(undef, fft_grid) for _ in 1:nchunk]

	locks = [ReentrantLock() for _ in 1:nR, _ in 1:np_ij]

	return RhoR_nG{RhoT, typeof(ρk)}(nR, np_ij, ρRlist, pair_ij, ekρR, ρk,
		nchunk, idx_chunks, ρR_buffer, ρk_buffer, wavek_buffer, wavekq_buffer, locks,
	)
end

