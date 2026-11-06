function _wotrunc_wieps_q0_pair(::Val{false}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (ipjp, R31, ki, R24, jl) in predata.pair_idx_chunk[ichunk]
            ρR_eps_left = view(ρR_eps, :, R31, ki)
            ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
            predata.pair_value_iq[ipjp] = (predata.pair_value_bare_iq[ipjp] + ρR_eps_left ⋅ ρR_eps_right) * weight
        end
    end
    Threads.@threads for iR in eachindex(predata.pair_R)
        @inbounds begin
            pair_value_iR = view(predata.pair_value, :, iR, iω)
            pair_value_iR .+= predata.pair_value_iq
        end
    end
end
function _wotrunc_wieps_q0_pair(::Val{true}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
    invϵ = eps00 + 1
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (ipjp, R31, ki, R24, jl) in predata.pair_idx_chunk[ichunk]
            (ip, jp) = predata.pair_idx[ipjp]
            ρR_eps_left = view(ρR_eps, :, R31, ki)
            ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
            wing_left_sumG = ρR_eps_left ⋅ epsG0
            wing_right_sumG = ρR_eps_right[G0idx_q0]
            AR_left = view(predata.AR, :, R31, ki)
            AR_right = view(predata.AR, :, R24, jl)
            wing_mpc = _wing_mpc_cal(Val(predata.wingaverage), predata.wing_mpc_coff,
                AR_left, AR_right, wing_left_sumG, wing_right_sumG, q̂₀)
            predata.pair_value_head_mpc[ip, jp, iω] += predata.pair_value_bare_head_mpc[ip, jp] * invϵ * weight
            predata.pair_value_wing_mpc[ip, jp, iω] += wing_mpc * weight
            predata.pair_value_iq[ipjp] = (predata.pair_value_bare_iq[ipjp] + ρR_eps_left ⋅ ρR_eps_right +
                predata.pair_value_bare_head_mpc[ip, jp] * eps00 + wing_mpc) * weight
        end
    end
    Threads.@threads for iR in eachindex(predata.pair_R)
        @inbounds begin
            pair_value_iR = view(predata.pair_value, :, iR, iω)
            pair_value_iR .+= predata.pair_value_iq
        end
    end
end
function _wieps_not_q0_pair(predata, ρR_eps, ρR_eps_mul, q, iω)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (ipjp, R31, ki, R24, jl) in predata.pair_idx_chunk[ichunk]
            ρR_eps_left = view(ρR_eps, :, R31, ki)
            ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
            predata.pair_value_iq[ipjp] = predata.pair_value_bare_iq[ipjp] + ρR_eps_left ⋅ ρR_eps_right
        end
    end
    Threads.@threads for iR in eachindex(predata.pair_R)
        @inbounds begin
            pair_value_iR = view(predata.pair_value, :, iR, iω)
            phase = cispi(-2 * (q ⋅ predata.pair_R[iR]))
            pair_value_iR .+= predata.pair_value_iq .* phase
        end
    end
end
