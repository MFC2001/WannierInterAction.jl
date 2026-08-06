function _wotrunc_wieps_q0_direct(::Val{false}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
    invϵ = eps00 + 1.0
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (i, j, R31, ki, R24, jl) in predata.direct_idx_chunk[ichunk]
            ρR_eps_left = view(ρR_eps, :, R31, ki)
            ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
            predata.direct_value_head[i, j, iω] += predata.direct_value_bare_head[i, j] * invϵ * weight
            predata.direct_value_iq[i, j] = (predata.direct_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
                                             predata.direct_value_bare_head[i, j] * eps00) * weight
        end
    end
    Threads.@threads for iR in eachindex(predata.direct_R)
        @inbounds begin
            direct_value_iR = view(predata.direct_value,:,:,iR,iω)
            direct_value_iR .+= predata.direct_value_iq
        end
    end
end
function _wotrunc_wieps_q0_direct(::Val{true}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
    invϵ = eps00 + 1.0
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (i, j, R31, ki, R24, jl) in predata.direct_idx_chunk[ichunk]
            ρR_eps_left = view(ρR_eps, :, R31, ki)
            ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
            wing_left_sumG = ρR_eps_left ⋅ epsG0
            wing_right_sumG = ρR_eps_right[G0idx_q0]
            AR_left = view(predata.AR, :, R31, ki)
            AR_right = view(predata.AR, :, R24, jl)
            wing_mpc = _wing_mpc_cal(Val(predata.wingaverage), predata.wing_mpc_coff,
                AR_left, AR_right, wing_left_sumG, wing_right_sumG, q̂₀)
            predata.direct_value_head[i, j, iω] += predata.direct_value_bare_head[i, j] * invϵ * weight
            predata.direct_value_head_mpc[i, j, iω] += predata.direct_value_bare_head_mpc[i, j] * invϵ * weight
            predata.direct_value_wing_mpc[i, j, iω] += wing_mpc * weight
            predata.direct_value_iq[i, j] = (predata.direct_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
                                             predata.direct_value_bare_head[i, j] * eps00 +
                                             predata.direct_value_bare_head_mpc[i, j] * eps00 + wing_mpc) * weight
        end
    end
    Threads.@threads for iR in eachindex(predata.direct_R)
        @inbounds begin
            direct_value_iR = view(predata.direct_value,:,:,iR,iω)
            direct_value_iR .+= predata.direct_value_iq
        end
    end
end
function _wieps_not_q0_direct(predata, ρR_eps, ρR_eps_mul, q, iω)
    Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
        @inbounds for (i, j, R31, ki, R24, jl) in predata.direct_idx_chunk[ichunk]
            ρR_eps_left = view(ρR_eps, :, R31, ki)
            ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
            predata.direct_value_iq[i, j] = predata.direct_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right
        end
    end
    Threads.@threads for iR in eachindex(predata.direct_R)
        @inbounds begin
            direct_value_iR = view(predata.direct_value,:,:,iR,iω)
            phase = cispi(-2 * (q ⋅ predata.direct_R[iR]))
            direct_value_iR .+= predata.direct_value_iq .* phase
        end
    end
end
