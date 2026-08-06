function _witrunc_woeps_q0_direct(predata, ρR, head)
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
function _witrunc_wieps_q0_direct(predata, ρR_eps, ρR_eps_mul, iω, head, weight)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (i, j, R31, ki, R24, jl) in predata.direct_idx_chunk[ichunk]
            ρR_eps_left = view(ρR_eps, :, R31, ki)
            ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
            predata.direct_value_head[i, j, iω] += head * weight
            predata.direct_value_iq[i, j] = (predata.direct_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
                                             head - predata.direct_value_bare_head[i, j]) * weight
        end
    end
    Threads.@threads for iR in eachindex(predata.direct_R)
        @inbounds begin
            direct_value_iR = view(predata.direct_value,:,:,iR,iω)
            direct_value_iR .+= predata.direct_value_iq
        end
    end
end
