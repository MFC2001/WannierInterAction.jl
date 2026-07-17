@inbounds function _wotrunc_woeps_q0_pair(::Val{false}, predata, ρR)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (ip1, ip2, Ridx_left, li, Ridx_right, jk) in predata.pair_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.pair_value_bare_iq[ip1, ip2] = ρR_left ⋅ ρR_right
		end
	end
	Threads.@threads for iR in eachindex(predata.pair_R)
		pair_value_iR = view(predata.pair_value_bare,:,:,iR)
		pair_value_iR .+= predata.pair_value_bare_iq
	end
end
@inbounds function _wotrunc_woeps_q0_pair(::Val{true}, predata, ρR)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (ip1, ip2, Ridx_left, li, Ridx_right, jk) in predata.pair_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			AR_left = view(predata.AR, :, Ridx_left, li)
			AR_right = view(predata.AR, :, Ridx_right, jk)
			CR_left = predata.pair_isp0[ip2] ? predata.CR[Ridx_left, li] : zero(ComplexF64)
			CR_right = predata.pair_isp0[ip1] ? predata.CR[Ridx_right, jk] : zero(ComplexF64)
			predata.pair_value_bare_head_mpc[ip1, ip2] = predata.head_mpc_coff * (AR_left ⋅ AR_right / 3 + (conj(CR_left) + CR_right) / 6)
			predata.pair_value_bare_iq[ip1, ip2] = ρR_left ⋅ ρR_right +
												   predata.pair_value_bare_head_mpc[ip1, ip2]
		end
	end
	Threads.@threads for iR in eachindex(predata.pair_R)
		pair_value_iR = view(predata.pair_value_bare,:,:,iR)
		pair_value_iR .+= predata.pair_value_bare_iq
	end
end
@inbounds function _woeps_not_q0_pair(predata, ρR, q)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (ip1, ip2, Ridx_left, li, Ridx_right, jk) in predata.pair_idx_chunk[ichunk]
			ρR_left = view(ρR, :, Ridx_left, li)
			ρR_right = view(ρR, :, Ridx_right, jk)
			predata.pair_value_bare_iq[ip1, ip2] = ρR_left ⋅ ρR_right
		end
	end
	Threads.@threads for iR in eachindex(predata.pair_R)
		pair_value_iR = view(predata.pair_value_bare,:,:,iR)
		phase = cispi(-2 * (q ⋅ predata.pair_R[iR]))
		pair_value_iR .+= predata.pair_value_bare_iq .* phase
	end
end
