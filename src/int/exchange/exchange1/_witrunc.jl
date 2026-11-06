function _witrunc_woeps_q0_exchange1(predata, ρR)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		@inbounds for (i, j, iR, R31, ki, R24, jl) in predata.exchange1_idx_chunk[ichunk]
			ρR_left = view(ρR, :, R31, ki)
			ρR_right = view(ρR, :, R24, jl)
			predata.exchange1_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right
		end
	end
	predata.exchange1_value_bare .+= predata.exchange1_value_bare_iq
end
function _witrunc_wieps_q0_exchange1(predata, ρR_eps, ρR_eps_mul, iω, weight)
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		@inbounds for (i, j, iR, R31, ki, R24, jl) in predata.exchange1_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, R31, ki)
			ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
			predata.exchange1_value[i, j, iR, iω] += (predata.exchange1_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right) * weight
		end
	end
end
