abstract type AbstractRhoR{RhoT} end
function Base.getproperty(ρR::AbstractRhoR{RhoT}, name::Symbol) where {RhoT}
	if name === :nw
		return ρR.ρk.nw
	elseif name === :nk
		return ρR.ρk.nk
	elseif name === :nq
		return ρR.ρk.nq
	elseif name === :kgrid
		return ρR.ρk.kgrid
	elseif name === :qgrid
		return ρR.ρk.qgrid
	elseif name === :nG
		return ρR.ρk.nG
	elseif name === :nG_max
		return ρR.ρk.nG_max
	elseif name === :fft_grid
		return ρR.ρk.fft_grid
	else
		return getfield(ρR, name)
	end
end
function (ρR::AbstractRhoR{RhoT})(R::AbstractVector{<:Integer}, i::Integer, j::Integer, iq::Integer; kwargs...) where {RhoT}
	return ρR(ReducedCoordinates{Int}(R), Int(i), Int(j), Int(iq); kwargs...)
end
function (ρR::AbstractRhoR{RhoT})(ρR_data::AbstractArray{RhoT, 1}, R::AbstractVector{<:Integer}, i::Integer, j::Integer, iq::Integer; kwargs...) where {RhoT}
	return ρR(ρR_data, ReducedCoordinates{Int}(R), Int(i), Int(j), Int(iq); kwargs...)
end
function (ρR::AbstractRhoR{RhoT})(R::ReducedCoordinates{Int}, i::Int, j::Int, iq::Int; kwargs...) where {RhoT}
	ρR_data = Vector{RhoT}(undef, ρR.nG[iq])
	return ρR(ρR_data, R, i, j, iq; kwargs...)
end
function (ρR::AbstractRhoR{RhoT})(ρR_data::AbstractArray{RhoT, 1}, R::ReducedCoordinates{Int}, i::Int, j::Int, iq::Int; nchunk::Integer = Threads.nthreads()) where {RhoT}
	ik_chunks = _RhoR_split_tasks(ρR.nk, nchunk)
	l = ReentrantLock()
	nG_q = ρR.nG[iq]
	ρR_view = view(ρR_data, 1:nG_q)
	fill!(ρR_view, zero(RhoT))
	Threads.@threads for ichunk in Base.OneTo(nchunk)
		(ik_first, ik_else) = ik_chunks[ichunk]
		@inbounds if !iszero(ik_first)
			ρR_ichunk = Vector{RhoT}(undef, nG_q)
			ρk_ichunk = Vector{RhoT}(undef, nG_q)
			wavek_ichunk = Array{RhoT, 3}(undef, ρR.fft_grid)
			wavekq_ichunk = Array{RhoT, 3}(undef, ρR.fft_grid)
			ρR.ρk(ρk_ichunk, wavek_ichunk, wavekq_ichunk, ik_first, i, j, iq)
			phase = cispi(2 * (ρR.kgrid[ik_first] ⋅ R)) / ρR.nk
			ρR_ichunk .= ρk_ichunk .* phase
			for ik in ik_else
				ρR.ρk(ρk_ichunk, wavek_ichunk, wavekq_ichunk, ik, i, j, iq)
				phase = cispi(2 * (ρR.kgrid[ik] ⋅ R)) / ρR.nk
				ρR_ichunk .+= ρk_ichunk .* phase
			end
			lock(l)
			try
				ρR_view .+= ρR_ichunk
			finally
				unlock(l)
			end
		end
	end
	return ρR_data
end

include("./RhoR_nG.jl")
include("./RhoR_nGnknij.jl")
include("./RhoR_nGnR.jl")
include("./RhoR_nobuffer.jl")

"""
	RhoR(ρk::AbstractRhok{RhoT}, pair_ij::AbstractVector{<:NTuple{2, <:Integer}}, ρRlist::AbstractVector{<:ReducedCoordinates{<:Integer}};
		nchunk::Integer = Threads.nthreads(), ρRbuffer::Symbol = :auto) where {RhoT}

Used for high performance computing.
"""
function RhoR(ρk::AbstractRhok{RhoT},
	ρRlist::AbstractVector{<:ReducedCoordinates{<:Integer}},
	pair_ij::AbstractVector{<:NTuple{2, <:Integer}};
	nchunk::Integer = Threads.nthreads(),
	ρRbuffer::Symbol = :auto,
) where {RhoT}

	if isempty(pair_ij)
		error("pair_ij shouldn't be empty!")
	end
	pair_ij = NTuple{2, Int}.(collect(pair_ij))
	if isempty(ρRlist)
		ρRlist = [ReducedCoordinates{Int}(0, 0, 0)]
	else
		ρRlist = ReducedCoordinates{Int}.(collect(ρRlist))
	end

	nk = ρk.nk
	np_ij = length(pair_ij)
	nR = length(ρRlist)
	nG_rho_max = ρk.nG_max
	fft_grid = ρk.fft_grid

	if nR == 1
		ρRbuffer = :low
	end

	if ρRbuffer == :auto
		high_mem = sizeof(RhoT) * (nG_rho_max * nk * np_ij + prod(fft_grid) * nchunk * 2)
		medium_mem = sizeof(RhoT) * (nG_rho_max * (nR + 1) * nchunk + prod(fft_grid) * nchunk * 2)
		low_mem = sizeof(RhoT) * (nG_rho_max * 2 * nchunk + prod(fft_grid) * nchunk * 2)
		free_mem = Sys.free_memory()
		# 0.5GB的余量
		if free_mem > high_mem + 1024^2 * 512
			ρRbuffer = :high
		elseif free_mem > medium_mem + 1024^2 * 512
			ρRbuffer = :medium
		elseif free_mem > low_mem + 1024^2 * 512
			ρRbuffer = :low
		else
			println("There is a possibility of memory overflow.")
			ρRbuffer = :low
		end
	end

	if ρRbuffer == :high
		return RhoR_nGnknij(ρk, ρRlist, pair_ij; nchunk)
	elseif ρRbuffer == :medium
		return RhoR_nGnR(ρk, ρRlist, pair_ij; nchunk)
	elseif ρRbuffer == :low
		return RhoR_nG(ρk, ρRlist, pair_ij; nchunk)
	end
end
function _RhoR_split_tasks(nk::Integer, nchunk::Integer)
	ik_chunks = _split_tasks(nk, nchunk)
	ik_chunks_tuple = Vector{Tuple{Int, UnitRange{Int}}}(undef, nchunk)
	for ichunk in 1:nchunk
		iks = ik_chunks[ichunk]
		if isempty(iks)
			ik_chunks_tuple[ichunk] = (0, 1:0)
		else
			ik_chunks_tuple[ichunk] = (first(iks), (first(iks)+1):last(iks))
		end
	end
	return ik_chunks_tuple
end
function _RhoR_split_tasks_wistart(nk::Integer, npij::Integer, nchunk::Integer)
	idx_chunks = _split_tasks(nk * npij, nchunk)
	all_ikiws = [(ik, iw) for ik in 1:nk, iw in 1:npij]
	idx_chunks_tuple = Vector{Vector{Tuple{Int, UnitRange{Int}, Int}}}(undef, nchunk)
	for ichunk in 1:nchunk
		idx_chunks_tuple_ichunk = Vector{Tuple{Int, UnitRange{Int}, Int}}(undef, 0)
		if !isempty(idx_chunks[ichunk])
			ikiws_ichunk = all_ikiws[idx_chunks[ichunk]]
			(ik_start, iw_present) = ikiws_ichunk[1]
			for (i, (ik, iw)) in enumerate(ikiws_ichunk)
				if i == lastindex(ikiws_ichunk)
					push!(idx_chunks_tuple_ichunk, (ik_start, (ik_start+1):ik, iw_present))
				elseif ikiws_ichunk[i+1][2] ≠ iw_present
					push!(idx_chunks_tuple_ichunk, (ik_start, (ik_start+1):ik, iw_present))
					(ik_start, iw_present) = ikiws_ichunk[i+1]
				end
			end
		end
		idx_chunks_tuple[ichunk] = idx_chunks_tuple_ichunk
	end
	return idx_chunks_tuple
end
function _RhoR_split_tasks_wostart(nk::Integer, npij::Integer, nchunk::Integer)
	idx_chunks = _split_tasks(nk * npij, nchunk)
	all_ikijs = [(ik, ij) for ik in 1:nk, ij in 1:npij]
	idx_chunks_tuple = Vector{Vector{Tuple{UnitRange{Int}, Int}}}(undef, nchunk)
	for ichunk in 1:nchunk
		idx_chunks_tuple_ichunk = Vector{Tuple{UnitRange{Int}, Int}}(undef, 0)
		if !isempty(idx_chunks[ichunk])
			ikijs_ichunk = all_ikijs[idx_chunks[ichunk]]
			(ik_start, ij_present) = ikijs_ichunk[1]
			for (i, (ik, ij)) in enumerate(ikijs_ichunk)
				if i == lastindex(ikijs_ichunk)
					push!(idx_chunks_tuple_ichunk, (ik_start:ik, ij_present))
				elseif ikijs_ichunk[i+1][2] ≠ ij_present
					push!(idx_chunks_tuple_ichunk, (ik_start:ik, ij_present))
					(ik_start, ij_present) = ikijs_ichunk[i+1]
				end
			end
		end
		idx_chunks_tuple[ichunk] = idx_chunks_tuple_ichunk
	end
	return idx_chunks_tuple
end
