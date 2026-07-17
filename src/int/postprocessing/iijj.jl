function _postprocessing_woeps_iijj(::Val{:txt}, outfolder, predata)
	hr = _int_make_hr_nwnwnR(predata.iijj_value_bare, predata.iijj_R)
	max_imag = maximum(abs∘imag, hr.value)
	max_ratio = maximum(x->abs(imag(x)/real(x)), hr.value)
	comment = @sprintf("Maximum Im: %10.6f; Im/Re Ratio: %10.6f.", max_imag, max_ratio)
	write(joinpath(outfolder, "iijj_bare.dat"), hr, wannier90_hr; comment)
end
function _postprocessing_wieps_iijj(::Val{:txt}, outfolder, predata)
	hrs = _int_make_hrs_nwnwnRnω(predata.iijj_value, predata.iijj_R)
	for iω in eachindex(hrs)
		max_imag = maximum(abs∘imag, hrs[iω].value)
		max_ratio = maximum(x->abs(imag(x)/real(x)), hrs[iω].value)
		comment = @sprintf("Maximum Im: %10.6f; Im/Re Ratio: %10.6f.", max_imag, max_ratio)
		write(joinpath(outfolder, "iijj_$iω.dat"), hrs[iω], wannier90_hr; comment)
	end
end
function _postprocessing_woeps_iijj(::Val{:jld2}, outfolder, predata)
	hr = _int_make_hr_nwnwnR(predata.iijj_value_bare, predata.iijj_R)
	hr_head_mpc = _int_make_hr_nwnwnR(predata.iijj_value_bare_head_mpc, predata.iijj_R)
	max_imag = maximum(abs∘imag, hr.value)
	max_ratio = maximum(x->abs(imag(x)/real(x)), hr.value)
	jldsave(joinpath(outfolder, "iijj_bare.jld2");
		interaction = hr, max_imag, max_ratio,
		head_mpc = hr_head_mpc)
end
function _postprocessing_wieps_iijj(::Val{:jld2}, outfolder, predata)
	hrs = _int_make_hrs_nwnwnRnω(predata.iijj_value, predata.iijj_R)
	hrs_head_mpc = _int_make_hrs_nwnwnRnω(predata.iijj_value_head_mpc, predata.iijj_R)
	hrs_wing_mpc = _int_make_hrs_nwnwnRnω(predata.iijj_value_wing_mpc, predata.iijj_R)
	for iω in eachindex(hrs)
		max_imag = maximum(abs∘imag, hrs[iω].value)
		max_ratio = maximum(x->abs(imag(x)/real(x)), hrs[iω].value)
		jldsave(joinpath(outfolder, "iijj_$iω.jld2");
			interaction = hrs[iω], max_imag, max_ratio,
			head_mpc = hrs_head_mpc[iω],
			wing_mpc = hrs_wing_mpc[iω])
	end
end
