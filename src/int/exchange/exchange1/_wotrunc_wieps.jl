@inbounds function _wotrunc_wieps_q0_exchange1(::Val{false}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
	invϵ = eps00 + 1
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, R31, ki, R24, jl) in predata.exchange1_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, R31, ki)
			ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
			predata.exchange1_value[i, j, iR, iω] += (predata.exchange1_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right) * weight
		end
	end
end
@inbounds function _wotrunc_wieps_q0_exchange1(::Val{true}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
	invϵ = eps00 + 1
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, R31, ki, R24, jl) in predata.exchange1_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, R31, ki)
			ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
			wing_left_sumG = ρR_eps_left ⋅ epsG0
			wing_right_sumG = ρR_eps_right[G0idx_q0]
			AR_left = view(predata.AR, :, R31, ki)
			AR_right = view(predata.AR, :, R24, jl)
			wing_mpc = _wing_mpc_cal(Val(predata.wingaverage), predata.wing_mpc_coff,
				AR_left, AR_right, wing_left_sumG, wing_right_sumG, q̂₀)
			predata.exchange1_value_head_mpc[i, j, iR, iω] += predata.exchange1_value_bare_head_mpc[i, j, iR] * invϵ * weight
			predata.exchange1_value_wing_mpc[i, j, iR, iω] += wing_mpc * weight
			predata.exchange1_value[i, j, iR, iω] += (predata.exchange1_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right +
												 predata.exchange1_value_bare_head_mpc[i, j, iR] * eps00 + wing_mpc) * weight
		end
	end
end
@inbounds function _wieps_not_q0_exchange1(predata, ρR_eps, ρR_eps_mul, q, iω)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, R31, ki, R24, jl) in predata.exchange1_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, R31, ki)
			ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
			predata.exchange1_value[i, j, iR, iω] += predata.exchange1_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right
		end
	end
end
