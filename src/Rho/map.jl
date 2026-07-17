"""
	_kgrid_map(aimkdirects,	kdirects1, kdirects2, mode = +) -> Matrix{Int}

kdirects is Vector{ReducedCoordinates{Rational{Int}}}.
"""
function _kgrid_map(aimkdirects::Vector{ReducedCoordinates{Rational{Int}}},
	kdirects1::Vector{ReducedCoordinates{Rational{Int}}},
	kdirects2::Vector{ReducedCoordinates{Rational{Int}}},
	mode = +,
)
	# aimk_to_idx = Dict{ReducedCoordinates{Int}, Int}(k => i for (i, k) in enumerate(aimkdirects))

	nk1 = length(kdirects1)
	nk2 = length(kdirects2)

	mapmat = Matrix{Int}(undef, nk1, nk2)
	mapGwind = Matrix{ReducedCoordinates{Int}}(undef, nk1, nk2)
	Threads.@threads for I in CartesianIndices(mapmat)
		@inbounds begin
			(i, j) = Tuple(I)
			k1 = kdirects1[i]
			k2 = kdirects2[j]
			Tk = mode(k1, k2)

			idx = findfirst(aimkdirects) do k
				all(isinteger, k - Tk)
			end
			if isnothing(idx)
				mapmat[I] = 0
				mapGwind[I] = ReducedCoordinates{Int}(0, 0, 0)
			else
				mapmat[I] = idx
				mapGwind[I] = Tk - aimkdirects[idx]
			end
		end
	end
	return mapmat, mapGwind
end
"""
	_Glist_map(aimGlist, Glist1, Glist2, mode = +) -> Matrix{Int}

Glist is Vector{ReducedCoordinates{Int}}.
"""
function _Glist_map(aimGlist::Vector{ReducedCoordinates{Int}},
	Glist1::Vector{ReducedCoordinates{Int}},
	Glist2::Vector{ReducedCoordinates{Int}},
	mode = +,
)::Matrix{Int}
	aimG_to_idx = Dict{ReducedCoordinates{Int}, Int}(G => i for (i, G) in enumerate(aimGlist))

	nG1 = length(Glist1)
	nG2 = length(Glist2)
	mapmat = Matrix{Int}(undef, nG1, nG2)
	Threads.@threads for iG2 in 1:nG2
		@inbounds begin
			G2 = Glist2[iG2]
			for iG1 in 1:nG1
				G1 = Glist1[iG1]
				mapmat[iG1, iG2] = get(aimG_to_idx, mode(G1, G2), 0)
			end
		end
	end

	return mapmat
end
