@inbounds function _wotrunc_woeps_q0_ijij(::Val{false}, predata, ρR)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.ijij_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.ijij_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right
		end
	end
	predata.ijij_value_bare .+= predata.ijij_value_bare_iq
end
@inbounds function _wotrunc_woeps_q0_ijij(::Val{true}, predata, ρR)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.ijij_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			AR_left = view(predata.AR, :, Ridx_left, li)
			AR_right = view(predata.AR, :, Ridx_right, jk)
			predata.ijij_value_bare_head_mpc[i, j, iR] = predata.head_mpc_coff * (AR_left ⋅ AR_right / 3)
			predata.ijij_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right +
												   predata.ijij_value_bare_head_mpc[i, j, iR]
		end
	end
	predata.ijij_value_bare .+= predata.ijij_value_bare_iq
end
@inbounds function _woeps_not_q0_ijij(predata, ρR, q)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.ijij_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.ijij_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right
		end
	end
	predata.ijij_value_bare .+= predata.ijij_value_bare_iq
end
