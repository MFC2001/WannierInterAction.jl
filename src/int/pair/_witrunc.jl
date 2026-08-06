function _witrunc_woeps_q0_pair(predata, ρR)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (ip1, ip2, R31, ki, R24, jl) in predata.pair_idx_chunk[ichunk]
            ρR_left = view(ρR, :, R31, ki)
            ρR_right = view(ρR, :, R24, jl)
            predata.pair_value_bare_iq[ip1, ip2] = ρR_left ⋅ ρR_right
        end
    end
    Threads.@threads for iR in eachindex(predata.pair_R)
        @inbounds begin
            pair_value_iR = view(predata.pair_value_bare,:,:,iR)
            pair_value_iR .+= predata.pair_value_bare_iq
        end
    end
end
function _witrunc_wieps_q0_pair(predata, ρR_eps, ρR_eps_mul, iω, weight)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (ip1, ip2, R31, ki, R24, jl) in predata.pair_idx_chunk[ichunk]
            ρR_eps_left = view(ρR_eps, :, R31, ki)
            ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
            predata.pair_value_iq[ip1, ip2] = (predata.pair_value_bare_iq[ip1, ip2] + ρR_eps_left ⋅ ρR_eps_right) * weight
        end
    end
    Threads.@threads for iR in eachindex(predata.pair_R)
        @inbounds begin
            pair_value_iR = view(predata.pair_value,:,:,iR,iω)
            pair_value_iR .+= predata.pair_value_iq
        end
    end
end
