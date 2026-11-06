function _int_make_hr_nwnwnR(value::AbstractArray{T,3}, R, ndegen) where {T}
    (nw, _, nR) = size(value)
    N = nw * nw * nR

    is = Vector{Int}(undef, N)
    js = Vector{Int}(undef, N)
    iRs = Vector{Int}(undef, N)
    vals = Vector{T}(undef, N)

    n=0
    for iR in Base.OneTo(nR), j in Base.OneTo(nw), i in Base.OneTo(nw)
        n+=1
        is[n] = i
        js[n] = j
        iRs[n] = iR
        vals[n] = value[i, j, iR]
    end

    return HRFull(is, js, iRs, vals, nw, R; ndegen, inbound=true, dofilter=false)
end
function _int_make_hrs_nwnwnRnω(value::AbstractArray{T,4}, R, ndegen) where {T}
    (nw, _, nR, nω) = size(value)
    N = nw * nw * nR

    hrs = Vector{HRFull{T}}(undef, nω)
    for iω in 1:nω
        is = Vector{Int}(undef, N)
        js = Vector{Int}(undef, N)
        iRs = Vector{Int}(undef, N)
        vals = Vector{T}(undef, N)
        n = 0
        for iR in Base.OneTo(nR), j in Base.OneTo(nw), i in Base.OneTo(nw)
            n+=1
            is[n] = i
            js[n] = j
            iRs[n] = iR
            vals[n] = value[i, j, iR]
        end
        hrs[iω] = HRFull(is, js, iRs, vals, nw, R; ndegen, inbound=true, dofilter=false)
    end

    return hrs
end
