@inline function _wing_mpc_cal(::Val{true}, coff, AR_left, AR_right, wing_left_sumG, wing_right_sumG, q̂₀)
	return coff * (wing_left_sumG * (AR_right[1] + AR_right[2] + AR_right[3]) / 3 +
				   conj(AR_left[1] + AR_left[2] + AR_left[3]) * wing_right_sumG / 3)
end
@inline function _wing_mpc_cal(::Val{false}, coff, AR_left, AR_right, wing_left_sumG, wing_right_sumG, q̂₀)
	return coff * (wing_left_sumG * (AR_right[1] * q̂₀[1] + AR_right[2] * q̂₀[2] + AR_right[3] * q̂₀[3]) +
				   conj(AR_left[1] * q̂₀[1] + AR_left[2] * q̂₀[2] + AR_left[3] * q̂₀[3]) * wing_right_sumG)
end
