@inbounds function _witrunc_woeps_q0_iijj(predata, ρR)
	# eqR = [cispi(-2 * (q ⋅ R)) for R in predata.iijj_R]
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.iijj_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.iijj_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right
		end
	end
	predata.iijj_value_bare .+= predata.iijj_value_bare_iq
end
@inbounds function _witrunc_wieps_q0_iijj(predata, ρR_eps, ρR_eps_mul, iω, weight)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.iijj_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			predata.iijj_value[i, j, iR, iω] += (predata.iijj_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right) * weight
		end
	end
end
