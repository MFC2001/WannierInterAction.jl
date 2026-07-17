@inbounds function _wotrunc_wieps_q0_ijij(::Val{false}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
	invϵ = eps00 + 1
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.ijij_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			predata.ijij_value[i, j, iR, iω] += (predata.ijij_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right) * weight
		end
	end
end
@inbounds function _wotrunc_wieps_q0_ijij(::Val{true}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
	invϵ = eps00 + 1
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.ijij_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			wing_left_sumG = ρR_eps_left ⋅ epsG0
			wing_right_sumG = ρR_eps_right[G0idx_q0]
			AR_left = view(predata.AR, :, Ridx_left, li)
			AR_right = view(predata.AR, :, Ridx_right, jk)
			wing_mpc = _wing_mpc_cal(Val(predata.wingaverage), predata.wing_mpc_coff,
				AR_left, AR_right, wing_left_sumG, wing_right_sumG, q̂₀)
			predata.ijij_value_head_mpc[i, j, iR, iω] += predata.ijij_value_bare_head_mpc[i, j, iR] * invϵ * weight
			predata.ijij_value_wing_mpc[i, j, iR, iω] += wing_mpc * weight
			predata.ijij_value[i, j, iR, iω] += (predata.ijij_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right +
												 predata.ijij_value_bare_head_mpc[i, j, iR] * eps00 + wing_mpc) * weight
		end
	end
end
@inbounds function _wieps_not_q0_ijij(predata, ρR_eps, ρR_eps_mul, q, iω)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.ijij_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			predata.ijij_value[i, j, iR, iω] += predata.ijij_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right
		end
	end
end
