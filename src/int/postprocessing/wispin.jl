function _postprocessing_woeps_int(fileformat::Symbol, ρk::Union{Rhok_collinear, Rhok_noncollinear}, outfolder::AbstractString, predata)

	mkpath(outfolder)

	if DIRECT ∈ predata.interaction
		predata.direct_value_bare ./= ρk.nq
		predata.direct_value_bare_head ./= ρk.nq
		predata.direct_value_bare_head_mpc ./= ρk.nq
		_postprocessing_woeps_direct(Val(fileformat), outfolder, predata)
	end
	if PAIR ∈ predata.interaction
		predata.pair_value_bare ./= ρk.nq
		predata.pair_value_bare_head_mpc ./= ρk.nq
		_postprocessing_woeps_pair(Val(fileformat), outfolder, predata)
	end

	return nothing
end
function _postprocessing_wieps_int(fileformat::Symbol, ρk::Union{Rhok_collinear, Rhok_noncollinear}, outfolder::AbstractString, predata)

	mkpath(outfolder)
	# e^2/NΩϵ₀
	coff = 1000 * qₑ / (ρk.nq * ρk.Ω * ϵ₀)

	if DIRECT ∈ predata.interaction
		predata.direct_value_bare ./= ρk.nq
		predata.direct_value_bare_head ./= ρk.nq
		predata.direct_value_bare_head_mpc ./= ρk.nq
		_postprocessing_woeps_direct(Val(fileformat), outfolder, predata)
		predata.direct_value ./= ρk.nq
		predata.direct_value_head ./= ρk.nq
		predata.direct_value_head_mpc ./= ρk.nq
		predata.direct_value_wing_mpc ./= ρk.nq
		_postprocessing_wieps_direct(Val(fileformat), outfolder, predata)
	end
	if PAIR ∈ predata.interaction
		predata.pair_value_bare ./= ρk.nq
		predata.pair_value_bare_head_mpc ./= ρk.nq
		_postprocessing_woeps_pair(Val(fileformat), outfolder, predata)
		predata.pair_value ./= ρk.nq
		predata.pair_value_head_mpc ./= ρk.nq
		predata.pair_value_wing_mpc ./= ρk.nq
		_postprocessing_wieps_pair(Val(fileformat), outfolder, predata)
	end
	return nothing
end


