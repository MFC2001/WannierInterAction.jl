# 对于bare相互作用参数，仅涉及向量点乘，
function _wotrunc_woeps_q0_direct(::Val{false}, predata, ρR, head)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (i, j, R31, ki, R24, jl) in predata.direct_idx_chunk[ichunk]
            ρR_left = view(ρR, :, R31, ki)
            ρR_right = view(ρR, :, R24, jl)
            predata.direct_value_bare_head[i, j] = head
            predata.direct_value_bare_iq[i, j] = ρR_left ⋅ ρR_right + predata.direct_value_bare_head[i, j]
        end
    end
    Threads.@threads for iR in eachindex(predata.direct_R)
        @inbounds begin
            direct_value_iR = view(predata.direct_value_bare,:,:,iR)
            direct_value_iR .+= predata.direct_value_bare_iq
        end
    end
end
@inbounds function _wotrunc_woeps_q0_direct(::Val{true}, predata, ρR, head)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (i, j, R31, ki, R24, jl) in predata.direct_idx_chunk[ichunk]
            ρR_left = view(ρR, :, R31, ki)
            ρR_right = view(ρR, :, R24, jl)
            AR_left = view(predata.AR, :, R31, ki)
            AR_right = view(predata.AR, :, R24, jl)
            CR_left = predata.CR[R31, ki]
            CR_right = predata.CR[R24, jl]
            predata.direct_value_bare_head[i, j] = head
            predata.direct_value_bare_head_mpc[i, j] = predata.head_mpc_coff * (AR_left ⋅ AR_right / 3 + (conj(CR_left) + CR_right) / 6)
            predata.direct_value_bare_iq[i, j] = ρR_left ⋅ ρR_right + predata.direct_value_bare_head[i, j] +
                                                 predata.direct_value_bare_head_mpc[i, j]
        end
    end
    Threads.@threads for iR in eachindex(predata.direct_R)
        @inbounds begin
            direct_value_iR = view(predata.direct_value_bare,:,:,iR)
            direct_value_iR .+= predata.direct_value_bare_iq
        end
    end
end
@inbounds function _woeps_not_q0_direct(predata, ρR, q)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (i, j, R31, ki, R24, jl) in predata.direct_idx_chunk[ichunk]
            ρR_left = view(ρR, :, R31, ki)
            ρR_right = view(ρR, :, R24, jl)
            predata.direct_value_bare_iq[i, j] = ρR_left ⋅ ρR_right
        end
    end
    Threads.@threads for iR in eachindex(predata.direct_R)
        @inbounds begin
            direct_value_iR = view(predata.direct_value_bare,:,:,iR)
            phase = cispi(-2 * (q ⋅ predata.direct_R[iR]))
            direct_value_iR .+= predata.direct_value_bare_iq .* phase
        end
    end
end
