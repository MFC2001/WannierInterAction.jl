function _postprocessing_onlyU(fileformat, ρk, outfolder, predata, head, ϵ_∞)
	# e^2/NΩϵ₀
	coff = 1000 * qₑ / (ρk.nq * ρk.Ω * ϵ₀)

	nω = length(ϵ_∞)
	# head term
	head .*= coff
	for iω in 1:nω
		println("iω = $iω, ϵ = $(1/ ϵ_∞[iω]), headϵ = $(head[iω]) eV.")
	end
	println("bare, ϵ = 1, head = $(head[end]) eV.")

	predata.ijji_value .*= coff
	_postprocessing_ijji(Val(fileformat), outfolder, predata.ijji_value, predata.ijji_R, head, ϵ_∞)

	return nothing
end
function _postprocessing_onlyU(fileformat, ρk, outfolder, predata, head)
	# e^2/NΩϵ₀
	coff = 1000 * qₑ / (ρk.nq * ρk.Ω * ϵ₀)

	# head term
	head .*= coff
	println("bare, ϵ = 1, head = $(head) eV.")

	predata.ijji_value .*= coff
	_postprocessing_ijji(Val(fileformat), outfolder, predata.ijji_value, predata.ijji_R, head)

	return nothing
end
