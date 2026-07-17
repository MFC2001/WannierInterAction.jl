function _postprocessing_int(::Val{:txt}, ρk::Rhok_collinear, outfolder, interaction, nw, predata, head, ϵ_∞)
	# e^2/NΩϵ₀
	coff = 1000 * qₑ / (ρk.nq * ρk.Ω * ϵ₀)

	nω = length(ϵ_∞)
	# head term
	head .*= coff
	for iω in 1:nω
		println("iω = $iω, ϵ = $(1/ ϵ_∞[iω]), headϵ = $(head[iω]) eV.")
	end
	println("bare, ϵ = 1, head = $(head[end]) eV.")

	IJJI ∈ interaction && begin
		predata.ijji_value .*= coff
		hrs = _int_make_hrs_forU(predata.ijji_value, predata.ijji_R)
		for iω in Base.OneTo(nω)
			write(joinpath(outfolder, "U_$iω.dat"), hrs[iω], wannier90_hr;
				comment = "ϵ = $(1/ ϵ_∞[iω]); head = $(head[iω]) eV. ")
		end
		write(joinpath(outfolder, "U_bare.dat"), hrs[end], wannier90_hr;
			comment = "ϵ = 1, head = $(head[end]) eV. ")
	end
	PAIR ∈ interaction && begin
		_postprocessing_pair!(predata)
		predata.pair_value .*= coff
		for iω in Base.OneTo(nω)
			int_pair_struct = WanIntijklR(view(predata.pair_value,iω,:,:,:),
				predata.pair, predata.pair_pR, predata.pair_R; nw, atol = predata.pair_atol)
			write(joinpath(outfolder, "pair_$iω.dat"), int_pair_struct, WanIntijklR)
		end
		int_pair_struct = WanIntijklR(view(predata.pair_value,(nω+1),:,:,:),
			predata.pair, predata.pair_pR, predata.pair_R; nw, atol = predata.pair_atol)
		write(joinpath(outfolder, "pair_bare.dat"), int_pair_struct, WanIntijklR)
	end

	return nothing
end
function _postprocessing_int(::Val{:jld2}, ρk::Rhok_collinear, outfolder, interaction, nw, predata, head, ϵ_∞)
	# e^2/NΩϵ₀
	coff = 1000 * qₑ / (ρk.nq * ρk.Ω * ϵ₀)

	nω = length(ϵ_∞)
	# head term
	head .*= coff
	for iω in 1:nω
		println("iω = $iω, ϵ = $(1/ ϵ_∞[iω]), headϵ = $(head[iω]) eV.")
	end
	println("bare, ϵ = 1, head = $(head[end]) eV.")

	IJJI ∈ interaction && begin
		predata.ijji_value .*= coff
		hrs = _int_make_hrs_forU(predata.ijji_value, predata.ijji_R)
		for iω in Base.OneTo(nω)
			jldsave(joinpath(outfolder, "U_$iω.jld2"); interaction = hrs[iω], ϵ = 1 / ϵ_∞[iω], head = head[iω])
		end
		jldsave(joinpath(outfolder, "U_bare.jld2"); interaction = hrs[end], ϵ = 1, head = head[end])
	end
	PAIR ∈ interaction && begin
		_postprocessing_pair!(predata)
		predata.pair_value .*= coff
		for iω in Base.OneTo(nω)
			int_pair_struct = WanIntijklR(view(predata.pair_value,iω,:,:,:),
				predata.pair, predata.pair_pR, predata.pair_R; nw, atol = predata.pair_atol)
			jldsave(joinpath(outfolder, "pair_$iω.jld2"); interaction = int_pair_struct)
		end
		int_pair_struct = WanIntijklR(view(predata.pair_value,(nω+1),:,:,:),
			predata.pair, predata.pair_pR, predata.pair_R; nw, atol = predata.pair_atol)
		jldsave(joinpath(outfolder, "pair_bare.jld2"); interaction = int_pair_struct)
	end

	return nothing
end
function _postprocessing_int(::Val{:txt}, ρk::Rhok_collinear, outfolder, interaction, nw, predata, head)
	# e^2/NΩϵ₀
	coff = 1000 * qₑ / (ρk.nq * ρk.Ω * ϵ₀)

	# head term
	head *= coff
	println("bare, ϵ = 1, head = $(head) eV.")

	IJJI ∈ interaction && begin
		predata.ijji_value .*= coff
		hr = _int_make_hr_forU(predata.ijji_value, predata.ijji_R)
		write(joinpath(outfolder, "U_bare.dat"), hr, wannier90_hr; comment = "ϵ = 1, head = $(head) eV. ")
	end
	PAIR ∈ interaction && begin
		_postprocessing_pair!(predata)
		predata.pair_value .*= coff
		int_pair_struct = WanIntijklR(predata.pair_value,
			predata.pair, predata.pair_pR, predata.pair_R; nw, atol = predata.pair_atol)
		write(joinpath(outfolder, "pair_bare.dat"), int_pair_struct, WanIntijklR)
	end

	return nothing
end
function _postprocessing_int(::Val{:jld2}, ρk::Rhok_collinear, outfolder, interaction, nw, predata, head)
	# e^2/NΩϵ₀
	coff = 1000 * qₑ / (ρk.nq * ρk.Ω * ϵ₀)

	# head term
	head *= coff
	println("bare, ϵ = 1, head = $(head) eV.")

	IJJI ∈ interaction && begin
		predata.ijji_value .*= coff
		hr = _int_make_hr_forU(predata.ijji_value, predata.ijji_R)
		jldsave(joinpath(outfolder, "U_bare.jld2"); interaction = hr, ϵ = 1, head = head)
	end
	PAIR ∈ interaction && begin
		_postprocessing_pair!(predata)
		predata.pair_value .*= coff
		int_pair_struct = WanIntijklR(predata.pair_value,
			predata.pair, predata.pair_pR, predata.pair_R; nw, atol = predata.pair_atol)
		jldsave(joinpath(outfolder, "pair_bare.jld2"); interaction = int_pair_struct)
	end

	return nothing
end
