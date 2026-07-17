function _wotrunc_int(ρk::AbstractRhok{RhoT}, eps::DielectricPWBasis, coulomb, predata, ρR, iq) where {RhoT}
	_wotrunc_int(ρk, coulomb, predata, ρR, iq)
	nG_eps_iq = eps.nG[iq]
	ρR_eps = ρR[1:nG_eps_iq, :, :]
	ρR_eps_mat = reshape(ρR_eps, nG_eps_iq, :) # a view
	ρR_eps_mul = similar(ρR_eps)
	ρR_eps_mul_mat = reshape(ρR_eps_mul, nG_eps_iq, :) # a view
	if iq == ρk.q0_idx
		G0idx_q0 = ρk.G0idx_q0
		nq₀ = length(eps.q₀)
		for iq₀ in 1:nq₀
			weight = 1 / nq₀
			q̂₀ = eps.q̂₀[iq₀]
			for iω in 1:eps.nω
				epsmat = eps.matrices_q0[iq₀, iω]
				BLAS.set_num_threads(predata.BLAS_nthread_epsilon)
				mul!(ρR_eps_mul_mat, epsmat, ρR_eps_mat)
				BLAS.set_num_threads(1)
				eps00 = epsmat[G0idx_q0, G0idx_q0]
				epsG0 = epsmat[:, G0idx_q0]
				predata.ϵ[iω] += (eps00 + 1.0) * weight
				IJJI ∈ predata.interaction && _wotrunc_wieps_q0_ijji(
					Val(predata.havempc), predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
				IJIJ ∈ predata.interaction && _wotrunc_wieps_q0_ijij(
					Val(predata.havempc), predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
				IIJJ ∈ predata.interaction && _wotrunc_wieps_q0_iijj(
					Val(predata.havempc), predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
				PAIR ∈ predata.interaction && _wotrunc_wieps_q0_pair(
					Val(predata.havempc), predata, ρR_eps, ρR_eps_mul, q̂₀, G0idx_q0, iω, eps00, epsG0, weight)
			end
		end
		predata.ϵ .= 1.0 ./ predata.ϵ
	else
		q = ρk.qgrid[iq]
		epsmat = eps[iq, :]
		for iω in 1:eps.nω
			BLAS.set_num_threads(predata.BLAS_nthread_epsilon)
			mul!(ρR_eps_mul_mat, epsmat[iω], ρR_eps_mat)
			BLAS.set_num_threads(1)
			IJJI ∈ predata.interaction && _wieps_not_q0_ijji(predata, ρR_eps, ρR_eps_mul, q, iω)
			IJIJ ∈ predata.interaction && _wieps_not_q0_ijij(predata, ρR_eps, ρR_eps_mul, q, iω)
			IIJJ ∈ predata.interaction && _wieps_not_q0_iijj(predata, ρR_eps, ρR_eps_mul, q, iω)
			PAIR ∈ predata.interaction && _wieps_not_q0_pair(predata, ρR_eps, ρR_eps_mul, q, iω)
		end
	end
end
