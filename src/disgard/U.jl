@inbounds function _wotrunc_wieps_q0_ijji(::Val{false}, predata, prebuffer, q, eps, nG_eps_iq)
	nρR_eps = length(predata.ijji_idx_ρR)
	ρR_eps = view(prebuffer.ρR_eps_right, 1:nG_eps_iq, 1:nρR_eps)
	ρR_eps_mul = view(prebuffer.ρR_eps_mul, 1:nG_eps_iq, 1:nρR_eps)
	for (idx, (ρR_idx, ij_idx)) in enumerate(predata.ijji_idx_ρR)
		ρR_eps_idx = view(ρR_eps, :, idx)
		ρR_eps_idx .= view(prebuffer.ρR, 1:nG_eps_iq, ρR_idx, ij_idx)
	end
	nq₀ = length(eps.q₀)
	G0idx_q0 = eps.minGidx[eps.q0_idx]
	for iq₀ in 1:nq₀
		epsmat = eps.matrices_q0[iq₀, :]
		eps00 = map(ϵ -> ϵ[G0idx_q0, G0idx_q0], epsmat)
		invϵ = eps00 .+ 1
		for iω in eachindex(epsmat)
			mul!(ρR_eps_mul, epsmat[iω], ρR_eps)
			for (i, j, ρR_left_idx, ρR_right_idx, _, _, _, _) in predata.ijji_idx
				ρR_eps_left = view(ρR_eps, :, ρR_left_idx)
				ρR_eps_right = view(ρR_eps_mul, :, ρR_right_idx)
				predata.ijji_value_head[i, j, iω] += predata.ijji_value_bare_head[i, j] * invϵ[iω]
				predata.ijji_value_iq[i, j, iω] = predata.ijji_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
												  predata.ijji_value_bare_head[i, j] * eps00[iω]
			end
		end
		predata.ijji_value_iq ./= nq₀
		Threads.@threads for iR in eachindex(predata.ijji_R)
			ijji_value_iR = view(predata.ijji_value,:,:,:,iR)
			ijji_value_iR .+= predata.ijji_value_iq
		end
	end
	predata.ijji_value_head ./= nq₀
end
@inbounds function _wotrunc_wieps_q0_ijji(::Val{true}, predata, prebuffer, q, eps, nG_eps_iq)
	nρR_eps = length(predata.ijji_idx_ρR)
	ρR_eps = view(prebuffer.ρR_eps_right, 1:nG_eps_iq, 1:nρR_eps)
	ρR_eps_mul = view(prebuffer.ρR_eps_mul, 1:nG_eps_iq, 1:nρR_eps)
	for (idx, (ρR_idx, ij_idx)) in enumerate(predata.ijji_idx_ρR)
		ρR_eps_idx = view(ρR_eps, :, idx)
		ρR_eps_idx .= view(prebuffer.ρR, 1:nG_eps_iq, ρR_idx, ij_idx)
	end
	nq₀ = length(eps.q₀)
	G0idx_q0 = eps.minGidx[eps.q0_idx]
	for iq₀ in 1:nq₀
		q̂₀ = eps.q̂₀[iq₀]
		epsmat = eps.matrices_q0[iq₀, :]
		eps00 = map(ϵ -> ϵ[G0idx_q0, G0idx_q0], epsmat)
		epsG0 = map(ϵ -> ϵ[:, G0idx_q0], epsmat)
		eps0G = map(ϵ -> ϵ[G0idx_q0, :], epsmat)
		invϵ = eps00 .+ 1
		for iω in eachindex(epsmat)
			mul!(ρR_eps_mul, epsmat[iω], ρR_eps)
			for (i, j, ρR_left_idx, ρR_right_idx, Ridx_left, li, Ridx_right, jk) in predata.ijji_idx
				ρR_eps_left = view(ρR_eps, :, ρR_left_idx)
				ρR_eps_right = view(ρR_eps_mul, :, ρR_right_idx)
				# 注意这里有问题，ρR_eps_right已经乘过节点矩阵了，直接取G0idx_q0即为wing_sumG
				wing_left_sumG, wing_right_sumG =
					_wing_cal(epsG0[iω], eps0G[iω], ρR_eps_left, ρR_eps_right, nG_eps_iq)
				AR_left = view(predata.AR, :, Ridx_left, li)
				AR_right = view(predata.AR, :, Ridx_right, jk)
				wing_mpc = _wing_mpc_cal(Val(predata.wingaverage), AR_left, AR_right, wing_left_sumG, wing_right_sumG, q̂₀)
				predata.ijji_value_head[i, j, iω] += predata.ijji_value_bare_head[i, j] * invϵ[iω]
				predata.ijji_value_iq[i, j, iω] = predata.ijji_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
												  predata.ijji_value_bare_head[i, j] * eps00[iω] + wing_mpc
			end
		end
		predata.ijji_value_iq ./= nq₀
		Threads.@threads for iR in eachindex(predata.ijji_R)
			ijji_value_iR = view(predata.ijji_value,:,:,:,iR)
			ijji_value_iR .+= predata.ijji_value_iq
		end
	end
	predata.ijji_value_head ./= nq₀
end
@inbounds function _wieps_not_q0_ijji(predata, prebuffer, q, epsmat, nG_eps_iq)
	# 该计算要求BLAS多线程，矩阵乘法的计算量远大于向量点乘。
	nρR_eps = length(predata.ijji_idx_ρR)
	ρR_eps = view(prebuffer.ρR_eps_right, 1:nG_eps_iq, 1:nρR_eps)
	ρR_eps_mul = view(prebuffer.ρR_eps_mul, 1:nG_eps_iq, 1:nρR_eps)
	for (idx, (ρR_idx, ij_idx)) in enumerate(predata.ijji_idx_ρR)
		ρR_eps_idx = view(ρR_eps, :, idx)
		ρR_eps_idx .= view(prebuffer.ρR, 1:nG_eps_iq, ρR_idx, ij_idx)
	end
	for iω in eachindex(epsmat)
		mul!(ρR_eps_mul, epsmat[iω], ρR_eps)
		# 大致的测试结果表明，1e6次循环5000维点乘耗时不超过1s，故使用单线程即可，如果BLAS多线程能加速，其会自动启用。
		for (i, j, ρR_left_idx, ρR_right_idx, _, _, _, _) in predata.ijji_idx
			ρR_eps_left = view(ρR_eps, :, ρR_left_idx)
			ρR_eps_right = view(ρR_eps_mul, :, ρR_right_idx)
			predata.ijji_value_iq[i, j, iω] = predata.ijji_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right
		end
	end
	# 这里广播后只余标量乘法，应该不会调用BLAS。
	Threads.@threads for iR in eachindex(predata.ijji_R)
		ijji_value_iR = view(predata.ijji_value,:,:,:,iR)
		phase = cispi(-2 * (q ⋅ predata.ijji_R[iR]))
		ijji_value_iR .+= predata.ijji_value_iq .* phase
	end
end




@inbounds function _witrunc_wieps_q0_ijji(predata, prebuffer, q, eps, nG_eps_iq, head)
	nρR_eps = length(predata.ijji_idx_ρR)
	ρR_eps = view(prebuffer.ρR_eps_right, 1:nG_eps_iq, 1:nρR_eps)
	ρR_eps_mul = view(prebuffer.ρR_eps_mul, 1:nG_eps_iq, 1:nρR_eps)
	for (idx, (ρR_idx, ij_idx)) in enumerate(predata.ijji_idx_ρR)
		ρR_eps_idx = view(ρR_eps, :, idx)
		ρR_eps_idx .= view(prebuffer.ρR, 1:nG_eps_iq, ρR_idx, ij_idx)
	end
	normq = map(q->sqrt(q ⋅ (eps.rldot * q)), eps.qgrid)
	sortq_idx = sortperm(normq)
	minq_idx = sortq_idx[2] # not Γ
	norm_minq = normq[minq_idx]
	G0idx_minq = eps.minGidx[minq_idx]
	ϵ_minq = map(iω -> eps.matrices[eps.qmap_idx[minq_idx], iω][G0idx_minq, G0idx_minq] + 1, 1:eps.nω)
	nq₀ = length(eps.q₀)
	G0idx_q0 = eps.minGidx[eps.q0_idx]
	for iq₀ in 1:nq₀
		q₀ = eps.q₀[iq₀]
		normq₀ = sqrt(q₀ ⋅ q₀)
		epsmat = eps.matrices_q0[iq₀, :]
		eps00 = map(ϵ -> ϵ[G0idx_q0, G0idx_q0], epsmat)
		if normq₀ < 1e-8
			normq_head = norm_minq
			ϵ_head = ϵ_minq
		else
			normq_head = normq₀
			ϵ_head = eps00 .+ 1
		end
		for iω in eachindex(epsmat)
			head_value = head(Val(:wocoff), ϵ_head[iω], normq_head)
			mul!(ρR_eps_mul, epsmat[iω], ρR_eps)
			for (i, j, ρR_left_idx, ρR_right_idx, _, _, _, _) in predata.ijji_idx
				ρR_eps_left = view(ρR_eps, :, ρR_left_idx)
				ρR_eps_right = view(ρR_eps_mul, :, ρR_right_idx)
				predata.ijji_value_head[i, j, iω] += head_value
				predata.ijji_value_iq[i, j, iω] = predata.ijji_value_bare_iq[i, j] + ρR_eps_left ⋅ ρR_eps_right +
												  head_value - predata.ijji_value_bare_head[i, j]
			end
		end
		predata.ijji_value_iq ./= nq₀
		Threads.@threads for iR in eachindex(predata.ijji_R)
			ijji_value_iR = view(predata.ijji_value,:,:,:,iR)
			ijji_value_iR .+= predata.ijji_value_iq
		end
	end
	predata.ijji_value_head ./= nq₀
end
