struct RhoR_nGnR{RhoT <: Number, RhokT <: AbstractRhok{RhoT}} <: AbstractRhoR{RhoT}
	nR::Int
	np_ij::Int
	ρRlist::Vector{ReducedCoordinates{Int}}
	pair_ij::Vector{NTuple{2, Int}}
	ekρR::Matrix{Complex{Float64}}
	ρk::RhokT
	# buffer
	nchunk::Int
	idx_chunks::Vector{Vector{Tuple{Int, UnitRange{Int}, Int}}} # ik_start, ik_else, ij
	ρR_buffer::Vector{Matrix{RhoT}}
	ρk_buffer::Vector{Vector{RhoT}}
	wavek_buffer::Vector{Array{RhoT, 3}}
	wavekq_buffer::Vector{Array{RhoT, 3}}
	locks::Vector{ReentrantLock}
end
function (ρR::RhoR_nGnR{RhoT, RhokT})(iq::Int) where {RhoT, RhokT}
	ρR_data = Array{RhoT}(undef, ρR.nG[iq], ρR.nR, ρR.np_ij)
	return ρR(ρR_data, iq)
end
function (ρR::RhoR_nGnR{RhoT, RhokT})(ρR_data::AbstractArray{RhoT, 3}, iq::Int) where {RhoT, RhokT}
	nG_q = ρR.nG[iq]
	ρR_view = view(ρR_data, 1:nG_q, 1:ρR.nR, 1:ρR.np_ij)
	fill!(ρR_view, zero(RhoT))
	Threads.@threads for ichunk in Base.OneTo(ρR.nchunk)
		@inbounds begin
			ρR_ichunk = view(ρR.ρR_buffer[ichunk], 1:nG_q, :)
			ρk_ichunk = view(ρR.ρk_buffer[ichunk], 1:nG_q)
			wavek_ichunk = ρR.wavek_buffer[ichunk]
			wavekq_ichunk = ρR.wavekq_buffer[ichunk]
			for (ik_first, ik_else, ij) in ρR.idx_chunks[ichunk]
				(i, j) = ρR.pair_ij[ij]
				ρR.ρk(ρk_ichunk, wavek_ichunk, wavekq_ichunk, ik_first, i, j, iq)
				ekρR_ik = transpose(view(ρR.ekρR, :, ik_first))
				ρR_ichunk .= ρk_ichunk .* ekρR_ik
				for ik in ik_else
					ρR.ρk(ρk_ichunk, wavek_ichunk, wavekq_ichunk, ik, i, j, iq)
					ekρR_ik = transpose(view(ρR.ekρR, :, ik))
					ρR_ichunk .+= ρk_ichunk .* ekρR_ik
				end
				ρR_ij = view(ρR_view,:,:,ij)
				lock(ρR.locks[ij])
				try
					ρR_ij .+= ρR_ichunk
				finally
					unlock(ρR.locks[ij])
				end
			end
		end
	end
	return ρR_data
end
# 通常来说，最高并行数为nk * np_ij
function RhoR_nGnR(ρk::AbstractRhok{RhoT},
	ρRlist::Vector{ReducedCoordinates{Int}},
	pair_ij::Vector{NTuple{2, Int}};
	nchunk::Integer = Threads.nthreads(),
) where {RhoT}

	nk = ρk.nk
	np_ij = length(pair_ij)
	nR = length(ρRlist)

	idx_chunks = _RhoR_split_tasks_wistart(nk, np_ij, nchunk)

	ekρR = [cispi(2 * (k ⋅ R)) / nk for R in ρRlist, k in ρk.kgrid]

	nG_rho_max = ρk.nG_max
	fft_grid = ρk.fft_grid

	# ρR buffer
	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nR * nchunk)
	println("ρR: $RhoT, dims=($nG_rho_max, $nR), nchunk=$nchunk, total=$ρR_size")
	# ρk buffer
	ρk_size = format_memory(sizeof(RhoT) * nG_rho_max * nchunk)
	println("ρk: $RhoT, len=$nG_rho_max, nchunk=$nchunk, total=$ρk_size")
	# Wave functions (ψk & ψkq)
	wave_size = format_memory(sizeof(RhoT) * prod(fft_grid))
	total_size = format_memory(sizeof(RhoT) * prod(fft_grid) * nchunk * 2)
	println("ψ_{k} and ψ_{k+q}: $RhoT, FFT_grid=$fft_grid, single=$wave_size, nchunk=$nchunk, total=$total_size")

	flush(stdout)

	ρR_buffer = [Matrix{RhoT}(undef, nG_rho_max, nR) for _ in 1:nchunk]
	ρk_buffer = [Vector{RhoT}(undef, nG_rho_max) for _ in 1:nchunk]
	wavek_buffer = [Array{RhoT, 3}(undef, fft_grid) for _ in 1:nchunk]
	wavekq_buffer = [Array{RhoT, 3}(undef, fft_grid) for _ in 1:nchunk]

	locks = [ReentrantLock() for _ in 1:np_ij]

	return RhoR_nGnR{RhoT, typeof(ρk)}(nR, np_ij, ρRlist, pair_ij, ekρR, ρk,
		nchunk, idx_chunks, ρR_buffer, ρk_buffer, wavek_buffer, wavekq_buffer, locks,
	)
end
