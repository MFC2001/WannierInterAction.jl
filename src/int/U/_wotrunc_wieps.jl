@inbounds function _wotrunc_wieps_q0_ijji(::Val{false}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
	invϵ = eps00 + 1.0
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, Ridx_left, li, Ridx_right, jk) in predata.ijji_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			predata.ijji_value_head[i, j, iω] += predata.ijji_value_bare_head[i, j] * invϵ * weight
			predata.ijji_value_iq[i, j] = (predata.ijji_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
										   predata.ijji_value_bare_head[i, j] * eps00) * weight
		end
	end
	Threads.@threads for iR in eachindex(predata.ijji_R)
		ijji_value_iR = view(predata.ijji_value,:,:,iR,iω)
		ijji_value_iR .+= predata.ijji_value_iq
	end
end
@inbounds function _wotrunc_wieps_q0_ijji(::Val{true}, predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
	invϵ = eps00 + 1.0
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, Ridx_left, li, Ridx_right, jk) in predata.ijji_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			wing_left_sumG = ρR_eps_left ⋅ epsG0
			wing_right_sumG = ρR_eps_right[G0idx_q0]
			AR_left = view(predata.AR, :, Ridx_left, li)
			AR_right = view(predata.AR, :, Ridx_right, jk)
			wing_mpc = _wing_mpc_cal(Val(predata.wingaverage), predata.wing_mpc_coff,
				AR_left, AR_right, wing_left_sumG, wing_right_sumG, q̂₀)
			predata.ijji_value_head[i, j, iω] += predata.ijji_value_bare_head[i, j] * invϵ * weight
			predata.ijji_value_head_mpc[i, j, iω] += predata.ijji_value_bare_head_mpc[i, j] * invϵ * weight
			predata.ijji_value_wing_mpc[i, j, iω] += wing_mpc * weight
			predata.ijji_value_iq[i, j] = (predata.ijji_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
										   predata.ijji_value_bare_head[i, j] * eps00 +
										   predata.ijji_value_bare_head_mpc[i, j] * eps00 + wing_mpc) * weight
		end
	end
	Threads.@threads for iR in eachindex(predata.ijji_R)
		ijji_value_iR = view(predata.ijji_value,:,:,iR,iω)
		ijji_value_iR .+= predata.ijji_value_iq
	end
end
@inbounds function _wieps_not_q0_ijji(predata, ρR_eps, ρR_eps_mul, q, iω)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, Ridx_left, li, Ridx_right, jk) in predata.ijji_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			predata.ijji_value_iq[i, j] = predata.ijji_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right
		end
	end
	Threads.@threads for iR in eachindex(predata.ijji_R)
		ijji_value_iR = view(predata.ijji_value,:,:,iR,iω)
		phase = cispi(-2 * (q ⋅ predata.ijji_R[iR]))
		ijji_value_iR .+= predata.ijji_value_iq .* phase
	end
end
