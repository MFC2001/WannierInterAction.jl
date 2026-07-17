function _preprocessing_pair(pair, pair_range, nw, grid, lattice, uplo)
	if pair_range == :default
		grid = (grid .+ 1) .* 2 .÷ 3
		pair_R = gridindex(Val(:WignerSeitz), lattice, grid)
	elseif pair_range isa NTuple{3, <:Integer} || pair_range isa AbstractVector{<:Integer}
		pair_R = gridindex(pair_range)
	elseif pair_range isa AbstractVector
		pair_R = ReducedCoordinates{Int}.(pair_range)
	else
		error("Wrong pair_range!")
	end

	if pair == :default
		pair = vec([(iw1, iw2, ReducedCoordinates{Int}(0, 0, 0)) for iw1 in 1:nw, iw2 in 1:nw])
	elseif pair isa NTuple{3, <:Integer} || pair isa AbstractVector{<:Integer}
		pair_pR = gridindex(pair)
		pair = vec([(iw1, iw2, pR) for iw1 in 1:nw, iw2 in 1:nw, pR in pair_pR])
	elseif pair isa AbstractVector{<:Union{NTuple{5, <:Integer}, AbstractVector{<:Integer}}}
		pair = map(pair) do p
			(p[1], p[2], ReducedCoordinates{Int}(p[3], p[4], p[5]))
		end
	else
		error("Wrong pair!")
	end

	pair = completepair(pair)
	pair_idx = uplopairidx(pair, uplo)

	pair_pR = unique(map(p->p[3], pair))
	pair_pR_to_idx = Dict{ReducedCoordinates{Int}, Int}(R => i for (i, R) in enumerate(pair_pR))
	pair = map(p->(p[1], p[2], pair_pR_to_idx[p[3]]), pair)

	return pair, pair_R, pair_pR, pair_idx
end
function _preprocessing_pair(pair, pair_range, nw, grid, lattice, uplo, upindex, dnindex)
	if pair_range == :default
		grid = (grid .+ 1) .÷ 2
		pair_R = gridindex(Val(:WignerSeitz), lattice, grid)
	elseif pair_range isa NTuple{3, <:Integer} || pair_range isa AbstractVector{<:Integer}
		pair_R = gridindex(pair_range)
	elseif pair_range isa AbstractVector
		pair_R = map(pair_range) do r
			ReducedCoordinates{Int}(r...)
		end
	else
		error("Wrong pair_range!")
	end

	if pair == :default
		pair = vec([(iw1, iw2, ReducedCoordinates{Int}(0, 0, 0)) for iw1 in 1:nw, iw2 in 1:nw])
	elseif pair isa NTuple{3, <:Integer} || pair isa AbstractVector{<:Integer}
		pair_pR = gridindex(pair)
		pair = vec([(iw1, iw2, pR) for iw1 in 1:nw, iw2 in 1:nw, pR in pair_pR])
	elseif pair isa AbstractVector{<:Union{NTuple{5, <:Integer}, AbstractVector{<:Integer}}}
		pair = map(pair) do p
			(p[1], p[2], ReducedCoordinates{Int}(p[3], p[4], p[5]))
		end
	else
		error("Wrong pair!")
	end

	# delete pair which have different spin.
	pair = filter(pair) do p
		!((p[1] ∈ upindex && p[2] ∈ dnindex) || (p[1] ∈ dnindex && p[2] ∈ upindex))
	end

	pair = completepair(pair)
	pair_idx = uplopairidx(pair, uplo)

	pair_pR = unique(map(p->p[3], pair))
	pair_pR_to_idx = Dict{ReducedCoordinates{Int}, Int}(R => i for (i, R) in enumerate(pair_pR))
	pair = map(p->(p[1], p[2], pair_pR_to_idx[p[3]]), pair)

	return pair, pair_R, pair_pR, pair_idx
end
