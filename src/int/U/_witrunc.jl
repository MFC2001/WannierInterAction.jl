@inbounds function _witrunc_woeps_q0_ijji(predata, ρR, head)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, Ridx_left, li, Ridx_right, jk) in predata.ijji_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.ijji_value_bare_head[i, j] = head
			predata.ijji_value_bare_iq[i, j] = ρR_left ⋅ ρR_right + predata.ijji_value_bare_head[i, j]
		end
	end
	Threads.@threads for iR in eachindex(predata.ijji_R)
		ijji_value_iR = view(predata.ijji_value_bare,:,:,iR)
		ijji_value_iR .+= predata.ijji_value_bare_iq
	end
end
@inbounds function _witrunc_wieps_q0_ijji(predata, ρR_eps, ρR_eps_mul, iω, head, weight)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, Ridx_left, li, Ridx_right, jk) in predata.ijji_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			predata.ijji_value_head[i, j, iω] += head * weight
			predata.ijji_value_iq[i, j] = (predata.ijji_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
										   head - predata.ijji_value_bare_head[i, j]) * weight
		end
	end
	Threads.@threads for iR in eachindex(predata.ijji_R)
		ijji_value_iR = view(predata.ijji_value,:,:,iR,iω)
		ijji_value_iR .+= predata.ijji_value_iq
	end
end
