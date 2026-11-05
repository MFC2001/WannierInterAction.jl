function _preprocessing_exchange(exchange_range, nw, period)
    if exchange_range isa NTuple{3,<:Integer} || exchange_range isa AbstractVector{<:Integer}
        exchange_range = [Int(x) for x in exchange_range]
        if !period[1]
            exchange_range[1] = 1
        elseif !period[2]
            exchange_range[2] = 1
        else
            exchange_range[3] = 1
        end
        exchange_R = gridindex(exchange_range)
    elseif exchange_range isa AbstractVector
        exchange_R = ReducedCoordinates{Int}.(exchange_range)
    else
        error("Wrong exchange_range!")
    end
    exchange_nR = length(exchange_R)
    exchange_R0idx = findfirst(R -> iszero(R), exchange_R)
    exchange_idx = [(iw1, iw2, iR) for iw1 in 1:nw, iw2 in 1:nw, iR in 1:exchange_nR]
    # delete direct-type!
    if !isnothing(exchange_R0idx)
        for iw in 1:nw
            exchange_idx[iw, iw, exchange_R0idx] = (0, 0, 0)
        end
        exchange_idx = vec(exchange_idx)
        filter!(idx -> !any(iszero, idx), exchange_idx)
    else
        exchange_idx = vec(exchange_idx)
    end
    return exchange_R, exchange_idx
end
