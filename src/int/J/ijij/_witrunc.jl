@inbounds function _witrunc_woeps_q0_ijij(predata, ρR)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.ijij_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.ijij_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right
		end
	end
	predata.ijij_value_bare .+= predata.ijij_value_bare_iq
end
@inbounds function _witrunc_wieps_q0_ijij(predata, ρR_eps, ρR_eps_mul, iω, weight)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, Ridx_left, li, Ridx_right, jk) in predata.ijij_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			predata.ijij_value[i, j, iR, iω] += (predata.ijij_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right) * weight
		end
	end
end
