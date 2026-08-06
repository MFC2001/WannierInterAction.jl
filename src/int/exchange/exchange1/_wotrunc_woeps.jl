@inbounds function _wotrunc_woeps_q0_exchange1(::Val{false}, predata, ρR)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, R31, ki, R24, jl) in predata.exchange1_idx_chunk[ichunk]
			ρR_left = view(ρR, :, R31, ki)
			ρR_right = view(ρR, :, R24, jl)
			predata.exchange1_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right
		end
	end
	predata.exchange1_value_bare .+= predata.exchange1_value_bare_iq
end
@inbounds function _wotrunc_woeps_q0_exchange1(::Val{true}, predata, ρR)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, R31, ki, R24, jl) in predata.exchange1_idx_chunk[ichunk]
			ρR_left = view(ρR, :, R31, ki)
			ρR_right = view(ρR, :, R24, jl)
			AR_left = view(predata.AR, :, R31, ki)
			AR_right = view(predata.AR, :, R24, jl)
			predata.exchange1_value_bare_head_mpc[i, j, iR] = predata.head_mpc_coff * (AR_left ⋅ AR_right / 3)
			predata.exchange1_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right +
												   predata.exchange1_value_bare_head_mpc[i, j, iR]
		end
	end
	predata.exchange1_value_bare .+= predata.exchange1_value_bare_iq
end
@inbounds function _woeps_not_q0_exchange1(predata, ρR, q)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		for (i, j, iR, R31, ki, R24, jl) in predata.exchange1_idx_chunk[ichunk]
			ρR_left = view(ρR, :, R31, ki)
			ρR_right = view(ρR, :, R24, jl)
			predata.exchange1_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right
		end
	end
	predata.exchange1_value_bare .+= predata.exchange1_value_bare_iq
end
