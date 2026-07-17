function _preprocessing_buffer(::Val{:wotrunc}, ρk::AbstractRhok{RhoT}, predata) where {RhoT}
	nG_rho_max = ρk.nG_max
	np_ij = length(predata.pair_ij)
	nρR = length(predata.ρRlist)

	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nρR * np_ij)
	println("ρR: $RhoT, dims=($nG_rho_max, $nρR, $np_ij), total=$ρR_size")

	flush(stdout)

	ρR = Array{RhoT}(undef, nG_rho_max, nρR, np_ij)

	return (ρR = ρR,)
end
function _preprocessing_buffer(::Val{:witrunc}, ρk::AbstractRhok{RhoT}, predata) where {RhoT}
	nG_rho_max = ρk.nG_max
	np_ij = length(predata.pair_ij)
	nρR = length(predata.ρRlist)

	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nρR * np_ij)
	println("ρR: $RhoT, dims=($nG_rho_max, $nρR, $np_ij), total=$ρR_size")

	flush(stdout)

	ρR = Array{RhoT}(undef, nG_rho_max, nρR, np_ij)

	return (ρR = ρR,)
end
function _preprocessing_buffer(::Val{:wotrunc}, eps, ρk::AbstractRhok{RhoT}, predata) where {RhoT}
	nG_rho_max = ρk.nG_max
	np_ij = length(predata.pair_ij)
	nρR = length(predata.ρRlist)

	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nρR * np_ij)
	println("ρR: $RhoT, dims=($nG_rho_max, $nρR, $np_ij), total=$ρR_size")

	flush(stdout)

	ρR = Array{RhoT}(undef, nG_rho_max, nρR, np_ij)

	return (ρR = ρR,)
end
function _preprocessing_buffer(::Val{:witrunc}, eps, ρk::AbstractRhok{RhoT}, predata) where {RhoT}
	nG_rho_max = ρk.nG_max
	np_ij = length(predata.pair_ij)
	nρR = length(predata.ρRlist)

	ρR_size = format_memory(sizeof(RhoT) * nG_rho_max * nρR * np_ij)
	println("ρR: $RhoT, dims=($nG_rho_max, $nρR, $np_ij), total=$ρR_size")

	flush(stdout)

	ρR = Array{RhoT}(undef, nG_rho_max, nρR, np_ij)

	return (ρR = ρR,)
end
