function _int3D_q0_wi_negq0_ijkl0(
	ijkl0_value, ijkl0_idx, epsmat, ρR_data,
	buffer_eps_mul, nG_eps_iq, nG_rho_iq, G0idx_q0, head_correction,
)
	tasks = map(enumerate(ijkl0_idx)) do (i_chunk, idxs)
		Threads.@spawn begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[i_chunk], 1:nG_eps_iq)
			for (i, j, k, l, ρR_left, li, ρR_right, jk) in idxs
				ρR_left_G0 = ρR_data[G0idx_q0, ρR_left, li]
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, ρR_left, li)
				ρR_left_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_left, li)
				ρR_right_G0 = ρR_data[G0idx_q0, ρR_right, jk]
				ρR_right_eps = view(ρR_data, 1:nG_eps_iq, ρR_right, jk)
				ρR_right_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_right, jk)
				head_correction_idx = head_correction * conj(ρR_left_G0) * ρR_right_G0
				for iω in eachindex(epsmat)
					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					ijkl0_value[iω, i, j, k, l] += (ρR_left_eps ⋅ buffer_eps_mul_chunk +
													head_correction_idx * epsmat[iω][G0idx_q0, G0idx_q0]) / 2
				end
				ijkl0_value[end, i, j, k, l] += (ρR_left_eps ⋅ ρR_right_eps + head_correction_idx) / 2
				high_dot = (ρR_left_high ⋅ ρR_right_high) / 2
				ijkl0_value[:, i, j, k, l] .+= high_dot
			end
		end
	end
	wait.(tasks)
end
function _int3D_q0_wo_negq0_ijkl0(
	ijkl0_value, ijkl0_idx, epsmat, ρR_data,
	buffer_eps_mul, nG_eps_iq, nG_rho_iq, G0idx_q0, head_correction,
)
	tasks = map(enumerate(ijkl0_idx)) do (i_chunk, idxs)
		Threads.@spawn begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[i_chunk], 1:nG_eps_iq)
			for (i, j, k, l, ρR_left, li, ρR_right, jk) in idxs
				ρR_left_G0 = ρR_data[G0idx_q0, ρR_left, li]
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, ρR_left, li)
				ρR_left_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_left, li)
				ρR_right_G0 = ρR_data[G0idx_q0, ρR_right, jk]
				ρR_right_eps = view(ρR_data, 1:nG_eps_iq, ρR_right, jk)
				ρR_right_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_right, jk)
				head_correction_idx = head_correction * conj(ρR_left_G0) * ρR_right_G0
				for iω in eachindex(epsmat)
					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					ijkl0_value[iω, i, j, k, l] += (ρR_left_eps ⋅ buffer_eps_mul_chunk +
													head_correction_idx * epsmat[iω][G0idx_q0, G0idx_q0])
				end
				ijkl0_value[end, i, j, k, l] += (ρR_left_eps ⋅ ρR_right_eps + head_correction_idx)
				high_dot = (ρR_left_high ⋅ ρR_right_high)
				ijkl0_value[:, i, j, k, l] .+= high_dot
			end
		end
	end
	wait.(tasks)
end
function _int3D_not_q0_ijkl0(
	ijkl0_value, ijkl0_idx, epsmat, ρR_data,
	buffer_eps_mul, nG_eps_iq, nG_rho_iq,
)
	tasks = map(enumerate(ijkl0_idx)) do (i_chunk, idxs)
		Threads.@spawn begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[i_chunk], 1:nG_eps_iq)
			for (i, j, k, l, ρR_left, li, ρR_right, jk) in idxs
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, ρR_left, li)
				ρR_left_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_left, li)
				ρR_right_eps = view(ρR_data, 1:nG_eps_iq, ρR_right, jk)
				ρR_right_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_right, jk)
				for iω in eachindex(epsmat)
					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					ijkl0_value[iω, i, j, k, l] += (ρR_left_eps ⋅ buffer_eps_mul_chunk)
				end
				ijkl0_value[end, i, j, k, l] += (ρR_left_eps ⋅ ρR_right_eps)
				high_dot = (ρR_left_high ⋅ ρR_right_high)
				ijkl0_value[:, i, j, k, l] .+= high_dot
			end
		end
	end
	wait.(tasks)
end
function _int2D_q0_wi_negq0_ijkl0(
	ijkl0_value, ijkl0_idx, epsmat, ρR_data, factor2D_iq,
	buffer_eps_mul, buffer_rho, nG_eps_iq, nG_rho_iq, G0idx_q0, head_correction,
)
	tasks = map(enumerate(ijkl0_idx)) do (i_chunk, idxs)
		Threads.@spawn begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[i_chunk], 1:nG_eps_iq)
			buffer_rho_chunk = view(buffer_rho[i_chunk], 1:nG_rho_iq)
			for (i, j, k, l, ρR_left, li, ρR_right, jk) in idxs
				ρR_left_G0 = ρR_data[G0idx_q0, ρR_left, li]
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, ρR_left, li)
				ρR_left_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_left, li)
				buffer_rho_chunk .= ρR_data[1:nG_rho_iq, ρR_right, jk] .* factor2D_iq
				ρR_right_G0 = buffer_rho_chunk[G0idx_q0]
				ρR_right_eps = view(buffer_rho_chunk, 1:nG_eps_iq)
				ρR_right_high = view(buffer_rho_chunk, nG_eps_iq+1:nG_rho_iq)
				head_correction_idx = head_correction * conj(ρR_left_G0) * ρR_right_G0
				for iω in eachindex(epsmat)
					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					ijkl0_value[iω, i, j, k, l] += (ρR_left_eps ⋅ buffer_eps_mul_chunk +
													head_correction_idx * epsmat[iω][G0idx_q0, G0idx_q0]) / 2
				end
				ijkl0_value[end, i, j, k, l] += (ρR_left_eps ⋅ ρR_right_eps + head_correction_idx) / 2
				high_dot = (ρR_left_high ⋅ ρR_right_high) / 2
				ijkl0_value[:, i, j, k, l] .+= high_dot
			end
		end
	end
	wait.(tasks)
end
function _int2D_q0_wo_negq0_ijkl0(
	ijkl0_value, ijkl0_idx, epsmat, ρR_data, factor2D_iq,
	buffer_eps_mul, buffer_rho, nG_eps_iq, nG_rho_iq, G0idx_q0, head_correction,
)
	tasks = map(enumerate(ijkl0_idx)) do (i_chunk, idxs)
		Threads.@spawn begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[i_chunk], 1:nG_eps_iq)
			buffer_rho_chunk = view(buffer_rho[i_chunk], 1:nG_rho_iq)
			for (i, j, k, l, ρR_left, li, ρR_right, jk) in idxs
				ρR_left_G0 = ρR_data[G0idx_q0, ρR_left, li]
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, ρR_left, li)
				ρR_left_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_left, li)
				buffer_rho_chunk .= ρR_data[1:nG_rho_iq, ρR_right, jk] .* factor2D_iq
				ρR_right_G0 = buffer_rho_chunk[G0idx_q0]
				ρR_right_eps = view(buffer_rho_chunk, 1:nG_eps_iq)
				ρR_right_high = view(buffer_rho_chunk, nG_eps_iq+1:nG_rho_iq)
				head_correction_idx = head_correction * conj(ρR_left_G0) * ρR_right_G0
				for iω in eachindex(epsmat)
					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					ijkl0_value[iω, i, j, k, l] += (ρR_left_eps ⋅ buffer_eps_mul_chunk +
													head_correction_idx * epsmat[iω][G0idx_q0, G0idx_q0])
				end
				ijkl0_value[end, i, j, k, l] += (ρR_left_eps ⋅ ρR_right_eps + head_correction_idx)
				high_dot = (ρR_left_high ⋅ ρR_right_high)
				ijkl0_value[:, i, j, k, l] .+= high_dot
			end
		end
	end
	wait.(tasks)
end
function _int2D_not_q0_ijkl0(
	ijkl0_value, ijkl0_idx, epsmat, ρR_data, factor2D_iq,
	buffer_eps_mul, buffer_rho, nG_eps_iq, nG_rho_iq,
)
	tasks = map(enumerate(ijkl0_idx)) do (i_chunk, idxs)
		Threads.@spawn begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[i_chunk], 1:nG_eps_iq)
			buffer_rho_chunk = view(buffer_rho[i_chunk], 1:nG_rho_iq)
			for (i, j, k, l, ρR_left, li, ρR_right, jk) in idxs
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, ρR_left, li)
				ρR_left_high = view(ρR_data, nG_eps_iq+1:nG_rho_iq, ρR_left, li)
				buffer_rho_chunk .= ρR_data[1:nG_rho_iq, ρR_right, jk] .* factor2D_iq
				ρR_right_eps = view(buffer_rho_chunk, 1:nG_eps_iq)
				ρR_right_high = view(buffer_rho_chunk, nG_eps_iq+1:nG_rho_iq)
				for iω in eachindex(epsmat)
					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					ijkl0_value[iω, i, j, k, l] += (ρR_left_eps ⋅ buffer_eps_mul_chunk)
				end
				ijkl0_value[end, i, j, k, l] += (ρR_left_eps ⋅ ρR_right_eps)
				high_dot = (ρR_left_high ⋅ ρR_right_high)
				ijkl0_value[:, i, j, k, l] .+= high_dot
			end
		end
	end
	wait.(tasks)
end
