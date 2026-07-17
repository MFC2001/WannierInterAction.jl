@inbounds function _witrunc_woeps_q0_pair(predata, ρR)
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
@inbounds function _witrunc_wieps_q0_pair(predata, ρR_eps, ρR_eps_mul, iω, weight)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (ip1, ip2, Ridx_left, li, Ridx_right, jk) in predata.pair_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, Ridx_left, li)
			ρR_eps_right = view(ρR_eps_mul, :, Ridx_right, jk)
			predata.pair_value_iq[ip1, ip2] = (predata.pair_value_bare_iq[ip1, ip2] + ρR_eps_left ⋅ ρR_eps_right) * weight
		end
	end
	Threads.@threads for iR in eachindex(predata.pair_R)
		pair_value_iR = view(predata.pair_value,:,:,iR,iω)
		pair_value_iR .+= predata.pair_value_iq
	end
end
