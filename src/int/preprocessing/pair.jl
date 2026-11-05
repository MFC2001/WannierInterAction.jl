function _preprocessing_pair(pair, pair_range, nw, lattice, grid, grid_atol, uplo)
    if pair_range == :default
        pair_R, pair_ndegen = gridindex(Val(:WignerSeitz), lattice, grid; atol=grid_atol)
    elseif pair_range isa NTuple{3,<:Integer} || pair_range isa AbstractVector{<:Integer}
        pair_R, pair_ndegen = gridindex(Val(:WignerSeitz), lattice, pair_range; atol=grid_atol)
        pair_ndegen = ones(Int, length(pair_R))
    elseif pair_range isa AbstractVector
        pair_R = ReducedCoordinates{Int}.(pair_range)
        pair_ndegen = ones(Int, length(pair_R))
    else
        error("Wrong pair_range!")
    end

    if pair == :default
        pair = vec([(i, j, ReducedCoordinates{Int}(0, 0, 0)) for i in 1:nw, j in 1:nw])
    elseif pair isa NTuple{3,<:Integer} || pair isa AbstractVector{<:Integer}
        pair_pR = gridindex(pair)
        pair = vec([(i, j, RR) for i in 1:nw, j in 1:nw, RR in pair_pR])
    elseif pair isa AbstractVector{<:WanIntPairVecType}
        pair = map(pair) do p
            (p[1], p[2], ReducedCoordinates{Int}(p[3]))
        end
    else
        error("Wrong pair!")
    end

    pair_pR = unique(map(p->p[3], pair))
    pR_to_idx = Dict{ReducedCoordinates{Int},Int}(R => i for (i, R) in enumerate(pair_pR))
    pair = map(p->(p[1], p[2], pR_to_idx[p[3]]), pair)

    pair_isonsite = map(p -> WannierInterActionBase._isonsite_pair(p, pair_pR), pair)

    np = length(pair)
    pair_idx = vec([(ip, jp) for ip in Base.OneTo(np), jp in Base.OneTo(np)])
    # delete direct-type!
    filter!(pair_idx) do (ip, jp)
        !(pair_isonsite[ip] && pair_isonsite[jp])
    end
    if uplo == :U
        filter!(pair_idx) do (ip, jp)
            ip ≤ jp
        end
    elseif uplo == :L
        filter!(pair_idx) do (ip, jp)
            ip ≥ jp
        end
    elseif uplo == :all
    else
        error("Wrong pair_uplo!")
    end

    # make sure that negR is included.
    if uplo ≠ :all
        negR = [-RR for RR in pair_R]
        allR = unique([pair_R; negR])
        sort!(allR)
        allR2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(allR))
        R_to_allR = [allR2idx[RR] for RR in pair_R]
        negR_to_allR = [allR2idx[RR] for RR in negR]
        allnedegen = zeros(Int, length(allR))
        for (iR, iallR) in enumerate(R_to_allR)
            allnedegen[iallR] = pair_ndegen[iR]
        end
        for (iR, iallR) in enumerate(negR_to_allR)
            if iszero(allnedegen[iallR])
                allnedegen[iallR] = pair_ndegen[iR]
            end
        end
        pair_R = allR
        pair_ndegen = allnedegen
    end

    return pair, pair_pR, pair_isonsite, pair_R, pair_ndegen, pair_idx
end
function _preprocessing_pair(pair, pair_range, nw, lattice, grid, grid_atol, uplo, upindex, dnindex)
    if pair_range == :default
        pair_R, pair_ndegen = gridindex(Val(:WignerSeitz), lattice, grid; atol=grid_atol)
    elseif pair_range isa NTuple{3,<:Integer} || pair_range isa AbstractVector{<:Integer}
        pair_R, pair_ndegen = gridindex(Val(:WignerSeitz), lattice, pair_range; atol=grid_atol)
        pair_ndegen = ones(Int, length(pair_R))
    elseif pair_range isa AbstractVector
        pair_R = ReducedCoordinates{Int}.(pair_range)
        pair_ndegen = ones(Int, length(pair_R))
    else
        error("Wrong pair_range!")
    end

    if pair == :default
        pair = vec([(i, j, ReducedCoordinates{Int}(0, 0, 0)) for i in 1:nw, j in 1:nw])
    elseif pair isa NTuple{3,<:Integer} || pair isa AbstractVector{<:Integer}
        pair_pR = gridindex(pair)
        pair = vec([(i, j, pR) for i in 1:nw, j in 1:nw, pR in pair_pR])
    elseif pair isa AbstractVector{<:WanIntPairVecType}
        pair = map(pair) do p
            (p[1], p[2], ReducedCoordinates{Int}(p[3]))
        end
    else
        error("Wrong pair!")
    end

    # delete pair which have different spin.
    pair = filter(pair) do p
        !((p[1] ∈ upindex && p[2] ∈ dnindex) || (p[1] ∈ dnindex && p[2] ∈ upindex))
    end

    pair_pR = unique(map(p->p[3], pair))
    pR_to_idx = Dict{ReducedCoordinates{Int},Int}(R => i for (i, R) in enumerate(pair_pR))
    pair = map(p->(p[1], p[2], pR_to_idx[p[3]]), pair)

    pair_isonsite = map(p -> WannierInterActionBase._isonsite_pair(p, pair_pR), pair)

    np = length(pair)
    pair_idx = vec([(ip, jp) for ip in Base.OneTo(np), jp in Base.OneTo(np)])
    # delete direct-type!
    filter!(pair_idx) do (ip, jp)
        !(pair_isonsite[ip] && pair_isonsite[jp])
    end
    if uplo == :U
        filter!(pair_idx) do (ip, jp)
            ip ≤ jp
        end
    elseif uplo == :L
        filter!(pair_idx) do (ip, jp)
            ip ≥ jp
        end
    elseif uplo == :all
    else
        error("Wrong pair_uplo!")
    end

    # make sure that negR is included.
    if uplo ≠ :all
        negR = [-RR for RR in pair_R]
        allR = unique([pair_R; negR])
        sort!(allR)
        allR2idx = Dict{ReducedCoordinates{Int},Int}(RR=>i for (i, RR) in enumerate(allR))
        R_to_allR = [allR2idx[RR] for RR in pair_R]
        negR_to_allR = [allR2idx[RR] for RR in negR]
        allnedegen = zeros(Int, length(allR))
        for (iR, iallR) in enumerate(R_to_allR)
            allnedegen[iallR] = pair_ndegen[iR]
        end
        for (iR, iallR) in enumerate(negR_to_allR)
            if iszero(allnedegen[iallR])
                allnedegen[iallR] = pair_ndegen[iR]
            end
        end
        pair_R = allR
        pair_ndegen = allnedegen
    end

    return pair, pair_pR, pair_isonsite, pair_R, pair_ndegen, pair_idx
end
