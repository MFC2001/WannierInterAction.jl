function _postprocessing_woeps_pair(::Val{:txt}, outfolder, predata)
    pairint = _postprocessing_pair(predata.pair_uplo, predata)
    # write(joinpath(outfolder, "pair_bare.dat"), pairint, WanIntijklR)
    jldsave(joinpath(outfolder, "pair_bare.jld2");
        interaction=pairint,
        head_mpc=predata.pair_value_bare_head_mpc)
end
function _postprocessing_wieps_pair(::Val{:txt}, outfolder, predata)
    nω = size(predata.pair_value, 4)
    for iω in Base.OneTo(nω)
        pairint = _postprocessing_pair(predata.pair_uplo, predata, iω)
        # write(joinpath(outfolder, "pair_$iω.dat"), pairint, WanIntijklR)
        jldsave(joinpath(outfolder, "pair_$iω.jld2");
            interaction=pairint,
            head_mpc=predata.pair_value_head_mpc[:, :, iω],
            wing_mpc=predata.pair_value_wing_mpc[:, :, iω])
    end
end
function _postprocessing_woeps_pair(::Val{:jld2}, outfolder, predata)
    pairint = _postprocessing_pair(predata.pair_uplo, predata)
    jldsave(joinpath(outfolder, "pair_bare.jld2");
        interaction=pairint,
        head_mpc=predata.pair_value_bare_head_mpc)
end
function _postprocessing_wieps_pair(::Val{:jld2}, outfolder, predata)
    nω = size(predata.pair_value, 4)
    for iω in Base.OneTo(nω)
        pairint = _postprocessing_pair(predata.pair_uplo, predata, iω)
        jldsave(joinpath(outfolder, "pair_$iω.jld2");
            interaction=pairint,
            head_mpc=predata.pair_value_head_mpc[:, :, iω],
            wing_mpc=predata.pair_value_wing_mpc[:, :, iω])
    end
end

_postprocessing_pair(uplo::Symbol, predata) =
    _postprocessing_pair(Val(uplo), predata.pair_value_bare, predata)
_postprocessing_pair(uplo::Symbol, predata, iω::Integer) =
    _postprocessing_pair(Val(uplo), view(predata.pair_value, :, :, iω), predata)
function _postprocessing_pair(::Val{:all}, value::AbstractMatrix{T}, predata) where {T}

    atol2 = predata.pair_atol * predata.pair_atol
    # count
    pair_nidx = length(predata.pair_idx)
    nR = length(predata.pair_R)
    N = 0
    @inbounds for iR in Base.OneTo(nR), ipjp in Base.OneTo(pair_nidx)
        if abs2(value[ipjp, iR]) > atol2
            N += 1
        end
    end

    ips = Vector{Int}(undef, N)
    jps = Vector{Int}(undef, N)
    iRs = Vector{Int}(undef, N)
    vals = Vector{T}(undef, N)

    n=0
    for iR in Base.OneTo(nR), ipjp in Base.OneTo(pair_nidx)
        v = value[ipjp, iR]
        if abs2(v) > atol2
            n+=1
            (ip, jp) = predata.pair_idx[ipjp]
            ips[n] = ip
            jps[n] = jp
            iRs[n] = iR
            vals[n] = v
        end
    end

    return PairIntHRFull(ips, jps, iRs, vals, predata.nw,
        predata.pair, predata.pair_pR, predata.pair_R; ndegen=predata.pair_ndegen,
        inbound=true, dofilter=true)
end
function _postprocessing_pair(::Val{:U}, value::AbstractMatrix{T}, predata) where {T}

    atol2 = predata.pair_atol * predata.pair_atol
    # count
    pair_nidx = length(predata.pair_idx)
    nR = length(predata.pair_R)
    N = 0
    @inbounds for iR in Base.OneTo(nR), ipjp in Base.OneTo(pair_nidx)
        if abs2(value[ipjp, iR]) > atol2
            N += 1
        end
    end

    ips = Vector{Int}(undef, N)
    jps = Vector{Int}(undef, N)
    iRs = Vector{Int}(undef, N)
    vals = Vector{T}(undef, N)

    n=0
    for iR in Base.OneTo(nR), ipjp in Base.OneTo(pair_nidx)
        v = value[ipjp, iR]
        if abs2(v) > atol2
            (ip, jp) = predata.pair_idx[ipjp]
            if ip == jp
                _isnegative_vector(predata.pair_R[iR]) && continue
                n+=1
                ips[n] = ip
                jps[n] = ip
                iRs[n] = iR
                vals[n] = v
            else
                n+=1
                ips[n] = ip
                jps[n] = jp
                iRs[n] = iR
                vals[n] = v
            end
        end
    end
    # ip==jp的项存在多余的，mask的计数偏大。
    resize!(ips, n)
    resize!(jps, n)
    resize!(iRs, n)
    resize!(vals, n)

    return PairIntHRHermitian(ips, jps, iRs, vals, predata.nw,
        predata.pair, predata.pair_pR, predata.pair_R; ndegen=predata.pair_ndegen,
        inbound=true, dofilter=true, have_negR=true, docheck=false)
end
function _postprocessing_pair(::Val{:L}, value::AbstractMatrix{T}, predata) where {T}

    R2idx = Dict(RR=>i for (i, RR) in enumerate(predata.pair_R))
    negR = collect(R2idx[-RR] for RR in predata.pair_R)

    atol2 = predata.pair_atol * predata.pair_atol
    # count
    pair_nidx = length(predata.pair_idx)
    nR = length(predata.pair_R)
    N = 0
    @inbounds for iR in Base.OneTo(nR), ipjp in Base.OneTo(pair_nidx)
        if abs2(value[ipjp, iR]) > atol2
            N += 1
        end
    end

    ips = Vector{Int}(undef, N)
    jps = Vector{Int}(undef, N)
    iRs = Vector{Int}(undef, N)
    vals = Vector{T}(undef, N)

    n=0
    for iR in Base.OneTo(nR), ipjp in Base.OneTo(pair_nidx)
        v = value[ipjp, iR]
        if abs2(v) > atol2
            (ip, jp) = predata.pair_idx[ipjp]
            if ip == jp
                _ispositive_vector(predata.pair_R[iR]) && continue
                n+=1
                ips[n] = ip
                jps[n] = ip
                iRs[n] = negR[iR]
                vals[n] = conj(v)
            else
                n+=1
                ips[n] = jp
                jps[n] = ip
                iRs[n] = negR[iR]
                vals[n] = conj(v)
            end
        end
    end
    # ip==jp的项存在多余的，mask的计数偏大。
    resize!(ips, n)
    resize!(jps, n)
    resize!(iRs, n)
    resize!(vals, n)

    return PairIntHRHermitian(ips, jps, iRs, vals, predata.nw,
        predata.pair, predata.pair_pR, predata.pair_R; ndegen=predata.pair_ndegen,
        inbound=true, dofilter=true, have_negR=true, docheck=false)
end
