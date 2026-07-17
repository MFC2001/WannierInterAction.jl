function _symmetrize_coulomb!(ρk, ρR, coulomb, iq)
	coulomb_iq = coulomb(ρk.qgrid[iq], ρk.Gidx[iq])
	if iq == ρk.q0_idx
		# 完全排除头项和翼项的影响
		coulomb_iq[ρk.G0idx_q0] = 0.0
	end
	# 对称化处理。
	coulomb_iq .= sqrt.(coulomb_iq)
	ρR_iq = view(ρR,(1:ρk.nG[iq]),:,:)
	ρR_iq .*= coulomb_iq
	return ρR_iq
end
