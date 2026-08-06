function _witrunc_int(ρk::AbstractRhok{RhoT}, eps::DielectricPWBasis, coulomb, predata, ρR, iq) where {RhoT}
	_witrunc_int(ρk, coulomb, predata, ρR, iq)
	nG_eps_iq = eps.nG[iq]
	ρR_eps = ρR[1:nG_eps_iq, :, :]
	ρR_eps_mat = reshape(ρR_eps, nG_eps_iq, :) # a view
	ρR_eps_mul = similar(ρR_eps)
	ρR_eps_mul_mat = reshape(ρR_eps_mul, nG_eps_iq, :) # a view
	if iq == ρk.q0_idx
		norm_minq, ϵ_minq = _witrunc_second_minq_eps(eps)
		G0idx_q0 = ρk.G0idx_q0
		nq₀ = length(eps.q₀)
		for iq₀ in 1:nq₀
			weight = 1 / nq₀
			q₀ = eps.q₀[iq₀]
			normq₀ = sqrt(q₀ ⋅ q₀)
			epsmat = eps.matrices_q0[iq₀, :]
			eps00 = map(ϵ -> ϵ[G0idx_q0, G0idx_q0], epsmat)
			if normq₀ < 1e-8
				head = coulomb(Val(:head), ρk.qgrid.kgrid_size, real.(ϵ_minq), norm_minq)
			else
				head = coulomb(Val(:head), ρk.qgrid.kgrid_size, real.(eps00) .+ 1, normq₀)
			end
			for iω in 1:eps.nω
				epsmat = eps.matrices_q0[iq₀, iω]
				BLAS.set_num_threads(predata.BLAS_nthread_epsilon)
				mul!(ρR_eps_mul_mat, epsmat, ρR_eps_mat)
				BLAS.set_num_threads(1)
				DIRECT ∈ predata.interaction && _witrunc_wieps_q0_direct(predata, ρR_eps, ρR_eps_mul, iω, head[iω], weight)
				EXCHANGE ∈ predata.interaction && _witrunc_wieps_q0_exchange1(predata, ρR_eps, ρR_eps_mul, iω, weight)
				EXCHANGE ∈ predata.interaction && _witrunc_wieps_q0_exchange2(predata, ρR_eps, ρR_eps_mul, iω, weight)
				PAIR ∈ predata.interaction && _witrunc_wieps_q0_pair(predata, ρR_eps, ρR_eps_mul, iω, weight)
			end
		end
	else
		q = ρk.qgrid[iq]
		epsmat = eps[iq, :]
		for iω in 1:eps.nω
			BLAS.set_num_threads(predata.BLAS_nthread_epsilon)
			mul!(ρR_eps_mul_mat, epsmat[iω], ρR_eps_mat)
			BLAS.set_num_threads(1)
			DIRECT ∈ predata.interaction && _wieps_not_q0_direct(predata, ρR_eps, ρR_eps_mul, q, iω)
			EXCHANGE ∈ predata.interaction && _wieps_not_q0_exchange1(predata, ρR_eps, ρR_eps_mul, q, iω)
			EXCHANGE ∈ predata.interaction && _wieps_not_q0_exchange2(predata, ρR_eps, ρR_eps_mul, q, iω)
			PAIR ∈ predata.interaction && _wieps_not_q0_pair(predata, ρR_eps, ρR_eps_mul, q, iω)
		end
	end
end
function _witrunc_second_minq_eps(eps)
	normq = map(q->sqrt(q ⋅ (eps.rldot * q)), eps.qgrid)
	sortq_idx = sortperm(normq)
	minq_idx = sortq_idx[2] # not Γ
	norm_minq = normq[minq_idx]
	G0idx_minq = eps.minGidx[minq_idx]
	ϵ_minq = map(iω -> eps.matrices[eps.qmap_idx[minq_idx], iω][G0idx_minq, G0idx_minq] + 1, 1:eps.nω)
	return norm_minq, ϵ_minq
end
