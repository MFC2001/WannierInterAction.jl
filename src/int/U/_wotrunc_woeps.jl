# 对于bare相互作用参数，仅涉及向量点乘，
@inbounds function _wotrunc_woeps_q0_ijji(::Val{false}, predata, ρR, head)
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
@inbounds function _wotrunc_woeps_q0_ijji(::Val{true}, predata, ρR, head)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, Ridx_left, li, Ridx_right, jk) in predata.ijji_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			AR_left = view(predata.AR, :, Ridx_left, li)
			AR_right = view(predata.AR, :, Ridx_right, jk)
			CR_left = predata.CR[Ridx_left, li]
			CR_right = predata.CR[Ridx_right, jk]
			predata.ijji_value_bare_head[i, j] = head
			predata.ijji_value_bare_head_mpc[i, j] = predata.head_mpc_coff * (AR_left ⋅ AR_right / 3 + (conj(CR_left) + CR_right) / 6)
			predata.ijji_value_bare_iq[i, j] = ρR_left ⋅ ρR_right + predata.ijji_value_bare_head[i, j] +
											   predata.ijji_value_bare_head_mpc[i, j]
		end
	end
	Threads.@threads for iR in eachindex(predata.ijji_R)
		ijji_value_iR = view(predata.ijji_value_bare,:,:,iR)
		ijji_value_iR .+= predata.ijji_value_bare_iq
	end
end
@inbounds function _woeps_not_q0_ijji(predata, ρR, q)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, Ridx_left, li, Ridx_right, jk) in predata.ijji_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.ijji_value_bare_iq[i, j] = ρR_left ⋅ ρR_right
		end
	end
	Threads.@threads for iR in eachindex(predata.ijji_R)
		ijji_value_iR = view(predata.ijji_value_bare,:,:,iR)
		phase = cispi(-2 * (q ⋅ predata.ijji_R[iR]))
		ijji_value_iR .+= predata.ijji_value_bare_iq .* phase
	end
end
