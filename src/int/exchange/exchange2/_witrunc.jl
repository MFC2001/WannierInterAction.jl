function _witrunc_woeps_q0_exchange2(predata, ρR)
	# eqR = [cispi(-2 * (q ⋅ R)) for R in predata.exchange_R]
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		@inbounds for (i, j, iR, R31, ki, R24, jl) in predata.exchange2_idx_chunk[ichunk]
			ρR_left = view(ρR, :, R31, ki)
			ρR_right = view(ρR, :, R24, jl)
			predata.exchange2_value_bare_iq[i, j, iR] = ρR_left ⋅ ρR_right
		end
	end
	predata.exchange2_value_bare .+= predata.exchange2_value_bare_iq
end
function _witrunc_wieps_q0_exchange2(predata, ρR_eps, ρR_eps_mul, iω, weight)
	# eqR = [cispi(-2 * (q ⋅ R)) for R in predata.exchange_R]
	Threads.@threads for ichunk in Base.OneTo(predata.JuliaTask_nchunk_int)
		@inbounds for (i, j, iR, R31, ki, R24, jl) in predata.exchange2_idx_chunk[ichunk]
			ρR_eps_left = view(ρR_eps, :, R31, ki)
			ρR_eps_right = view(ρR_eps_mul, :, R24, jl)
			predata.exchange2_value[i, j, iR, iω] += (predata.exchange2_value_bare_iq[i, j, iR] + ρR_eps_left ⋅ ρR_eps_right) * weight
		end
	end
end
