struct RhoR_nGnknij{RhoT <: Number, RhokT <: AbstractRhok{RhoT}} <: AbstractRhoR{RhoT}
	nR::Int
	np_ij::Int
	ρRlist::Vector{ReducedCoordinates{Int}}
	pair_ij::Vector{NTuple{2, Int}}
	ekρR::Matrix{Complex{Float64}}
	ρk::RhokT
	# buffer
	nchunk::Int
	idx_chunks::Vector{Vector{Tuple{UnitRange{Int}, Int}}} # iks, ij
	ρk_buffer::Array{RhoT, 3}
	wavek_buffer::Vector{Array{RhoT, 3}}
	wavekq_buffer::Vector{Array{RhoT, 3}}
end
function (ρR::RhoR_nGnknij{RhoT, RhokT})(iq::Int) where {RhoT, RhokT}
	ρR_data = Array{RhoT}(undef, ρR.nG[iq], ρR.nR, ρR.np_ij)
	return ρR(ρR_data, iq)
end
function (ρR::RhoR_nGnknij{RhoT, RhokT})(ρR_data::AbstractArray{RhoT, 3}, iq::Int) where {RhoT, RhokT}
	nG_q = ρR.nG[iq]
	Threads.@threads for ichunk in Base.OneTo(ρR.nchunk)
		@inbounds begin
			wavek_ichunk = ρR.wavek_buffer[ichunk]
			wavekq_ichunk = ρR.wavekq_buffer[ichunk]
			for (iks, ij) in ρR.idx_chunks[ichunk]
				(i, j) = ρR.pair_ij[ij]
				for ik in iks
					ρk_G = view(ρR.ρk_buffer, 1:nG_q, ik, ij)
					ρR.ρk(ρk_G, wavek_ichunk, wavekq_ichunk, ik, i, j, iq)
				end
			end
		end
	end
	Threads.@threads for ij in Base.OneTo(ρR.np_ij)
		ρk_G = view(ρR.ρk_buffer, 1:nG_q, :, ij)
		ρR_G = view(ρR_data, 1:nG_q, 1:ρR.nR, ij)
		mul!(ρR_G, ρk_G, ρR.ekρR)
	end
	return ρR_data
end
# 通常来说，最高并行数为nk * np_ij
function RhoR_nGnknij(ρk::AbstractRhok{RhoT},
	ρRlist::Vector{ReducedCoordinates{Int}},
	pair_ij::Vector{NTuple{2, Int}};
	nchunk::Integer = Threads.nthreads(),
) where {RhoT}

	nk = ρk.nk
	np_ij = length(pair_ij)
	nR = length(ρRlist)

	idx_chunks = _RhoR_split_tasks_wostart(nk, np_ij, nchunk)

	# ekρR = [cispi(2 * (k ⋅ R)) / nk for R in ρRlist, k in ρk.kgrid]
	ekρR = [cispi(2 * (k ⋅ R)) / nk for k in ρk.kgrid, R in ρRlist]

	nG_rho_max = ρk.nG_max
	fft_grid = ρk.fft_grid

	# ρk buffer
	ρk_size = format_memory(sizeof(RhoT) * nG_rho_max * nk * np_ij)
	println("ρk_buffer: $RhoT, dims=($nG_rho_max, $nk, $np_ij), total=$ρk_size")
	# Wave functions (ψk & ψkq)
	wave_size = format_memory(sizeof(RhoT) * prod(fft_grid))
	total_size = format_memory(sizeof(RhoT) * prod(fft_grid) * nchunk * 2)
	println("ψ_{k} and ψ_{k+q}: $RhoT, FFT_grid=$fft_grid, single=$wave_size, nchunk=$nchunk, total=$total_size")

	flush(stdout)

	ρk_buffer = Array{RhoT}(undef, nG_rho_max, nk, np_ij)
	wavek_buffer = [Array{RhoT, 3}(undef, fft_grid) for _ in 1:nchunk]
	wavekq_buffer = [Array{RhoT, 3}(undef, fft_grid) for _ in 1:nchunk]


	return RhoR_nGnknij{RhoT, typeof(ρk)}(nR, np_ij, ρRlist, pair_ij, ekρR, ρk,
		nchunk, idx_chunks, ρk_buffer, wavek_buffer, wavekq_buffer,
	)
end
