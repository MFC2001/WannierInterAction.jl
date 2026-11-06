function _witrunc_int(ρk::AbstractRhok{RhoT}, coulomb, predata, ρR, iq) where {RhoT}
	nG_rho_iq = ρk.nG[iq]
	ρR = view(ρR,(1:nG_rho_iq),:,:)
	if iq == ρk.q0_idx
		head = coulomb(Val(:head), ρk.qgrid.kgrid_size)
		DIRECT ∈ predata.interaction && _witrunc_woeps_q0_direct(predata, ρR, head)
		EXCHANGE ∈ predata.interaction && _witrunc_woeps_q0_exchange1(predata, ρR)
		EXCHANGE ∈ predata.interaction && _witrunc_woeps_q0_exchange2(predata, ρR)
		PAIR ∈ predata.interaction && _witrunc_woeps_q0_pair(predata, ρR)
	else
		q = ρk.qgrid[iq]
		DIRECT ∈ predata.interaction && _woeps_not_q0_direct(predata, ρR, q)
		EXCHANGE ∈ predata.interaction && _woeps_not_q0_exchange1(predata, ρR, q)
		EXCHANGE ∈ predata.interaction && _woeps_not_q0_exchange2(predata, ρR, q)
		PAIR ∈ predata.interaction && _woeps_not_q0_pair(predata, ρR, q)
	end
end
