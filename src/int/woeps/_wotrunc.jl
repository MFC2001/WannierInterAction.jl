function _wotrunc_int(ρk::AbstractRhok{RhoT}, coulomb, predata, ρR, iq) where {RhoT}
	nG_rho_iq = ρk.nG[iq]
	ρR = view(ρR,(1:nG_rho_iq),:,:)
	if iq == ρk.q0_idx
		head = coulomb(Val(:head), ρk.qgrid.kgrid_size)
		IJJI ∈ predata.interaction && _wotrunc_woeps_q0_ijji(Val(predata.havempc), predata, ρR, head)
		IJIJ ∈ predata.interaction && _wotrunc_woeps_q0_ijij(Val(predata.havempc), predata, ρR)
		IIJJ ∈ predata.interaction && _wotrunc_woeps_q0_iijj(Val(predata.havempc), predata, ρR)
		PAIR ∈ predata.interaction && _wotrunc_woeps_q0_pair(Val(predata.havempc), predata, ρR)
	else
		q = ρk.qgrid[iq]
		IJJI ∈ predata.interaction && _woeps_not_q0_ijji(predata, ρR, q)
		IJIJ ∈ predata.interaction && _woeps_not_q0_ijij(predata, ρR, q)
		IIJJ ∈ predata.interaction && _woeps_not_q0_iijj(predata, ρR, q)
		PAIR ∈ predata.interaction && _woeps_not_q0_pair(predata, ρR, q)
	end
end
