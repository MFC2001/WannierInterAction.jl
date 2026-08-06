function _postprocessing_woeps_direct(::Val{:txt}, outfolder, predata)
	hr = _int_make_hr_nwnwnR(predata.direct_value_bare, predata.direct_R)
	max_imag = maximum(abs∘imag, hr.value)
	max_ratio = maximum(x->abs(imag(x)/real(x)), hr.value)
	comment = @sprintf("ϵ = 1.0. Maximum Im: %10.6f; Im/Re Ratio: %10.6f.", max_imag, max_ratio)
	write(joinpath(outfolder, "U_bare.dat"), hr, wannier90_hr; comment)
end
function _postprocessing_wieps_direct(::Val{:txt}, outfolder, predata)
	hrs = _int_make_hrs_nwnwnRnω(predata.direct_value, predata.direct_R)
	for iω in eachindex(hrs)
		max_imag = maximum(abs∘imag, hrs[iω].value)
		max_ratio = maximum(x->abs(imag(x)/real(x)), hrs[iω].value)
		comment = @sprintf("ϵ = %10.6f. Maximum Im: %10.6f; Im/Re Ratio: %10.6f.", predata.ϵ[iω], max_imag, max_ratio)
		write(joinpath(outfolder, "U_$iω.dat"), hrs[iω], wannier90_hr; comment)
	end
end
function _postprocessing_woeps_direct(::Val{:jld2}, outfolder, predata)
	hr = _int_make_hr_nwnwnR(predata.direct_value_bare, predata.direct_R)
	max_imag = maximum(abs∘imag, hr.value)
	max_ratio = maximum(x->abs(imag(x)/real(x)), hr.value)
	jldsave(joinpath(outfolder, "U_bare.jld2");
		interaction = hr, max_imag, max_ratio,
		head = predata.direct_value_bare_head,
		head_mpc = predata.direct_value_bare_head_mpc,
		ϵ = 1.0)
end
function _postprocessing_wieps_direct(::Val{:jld2}, outfolder, predata)
	hrs = _int_make_hrs_nwnwnRnω(predata.direct_value, predata.direct_R)
	for iω in eachindex(hrs)
		max_imag = maximum(abs∘imag, hrs[iω].value)
		max_ratio = maximum(x->abs(imag(x)/real(x)), hrs[iω].value)
		jldsave(joinpath(outfolder, "U_$iω.jld2");
			interaction = hrs[iω], max_imag, max_ratio,
			head = predata.direct_value_head[:, :, iω],
			head_mpc = predata.direct_value_head_mpc[:, :, iω],
			wing_mpc = predata.direct_value_wing_mpc[:, :, iω],
			ϵ = predata.ϵ[iω])
	end
end
