@inbounds function _wotrunc_woeps_q0_iijj(::Val{false}, predata, ρR)
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
@inbounds function _wotrunc_woeps_q0_iijj(::Val{true}, predata, ρR)
	# eqR = [cispi(-2 * (q ⋅ R)) for R in predata.iijj_R]
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.iijj_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			AR_left = view(predata.AR, :, Ridx_left, li)
			AR_right = view(predata.AR, :, Ridx_right, jk)
			predata.iijj_value_bare_head_mpc[i, j, iR] = predata.head_mpc_coff * (AR_left ⋅ AR_right / 3)
			predata.iijj_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right +
												   predata.iijj_value_bare_head_mpc[i, j, iR]
		end
	end
	predata.iijj_value_bare .+= predata.iijj_value_bare_iq
end
@inbounds function _woeps_not_q0_iijj(predata, ρR, q)
	eqR = [cispi(-2 * (q ⋅ R)) for R in predata.iijj_R]
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.iijj_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.iijj_value_bare_iq[i, j, iR] = eqR[iR] * (ρR_left ⋅ ρR_right)
		end
	end
	predata.iijj_value_bare .+= predata.iijj_value_bare_iq
end
