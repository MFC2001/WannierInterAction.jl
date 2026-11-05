function _postprocessing_woeps_exchange1(::Val{:txt}, outfolder, predata)
    hr = _int_make_hr_nwnwnR(predata.exchange1_value_bare, predata.exchange_R, ones(Int, length(predata.exchange_R)))
    max_imag = maximum(abs∘imag, hr.vals)
    max_ratio = maximum(x->abs(imag(x)/real(x)), hr.vals)
    comment = @sprintf("Maximum Im: %10.6f; Im/Re Ratio: %10.6f.", max_imag, max_ratio)
    write(joinpath(outfolder, "J1_bare.dat"), hr, wannier90_hr; comment)
end
function _postprocessing_wieps_exchange1(::Val{:txt}, outfolder, predata)
    hrs = _int_make_hrs_nwnwnRnω(predata.exchange1_value, predata.exchange_R, ones(Int, length(predata.exchange_R)))
    for iω in eachindex(hrs)
        max_imag = maximum(abs∘imag, hrs[iω].vals)
        max_ratio = maximum(x->abs(imag(x)/real(x)), hrs[iω].vals)
        comment = @sprintf("Maximum Im: %10.6f; Im/Re Ratio: %10.6f.", max_imag, max_ratio)
        write(joinpath(outfolder, "J1_$iω.dat"), hrs[iω], wannier90_hr; comment)
    end
end
function _postprocessing_woeps_exchange1(::Val{:jld2}, outfolder, predata)
    hr = _int_make_hr_nwnwnR(predata.exchange1_value_bare, predata.exchange_R, ones(Int, length(predata.exchange_R)))
    hr_head_mpc = _int_make_hr_nwnwnR(predata.exchange1_value_bare_head_mpc, predata.exchange_R, ones(Int, length(predata.exchange_R)))
    max_imag = maximum(abs∘imag, hr.vals)
    max_ratio = maximum(x->abs(imag(x)/real(x)), hr.vals)
    jldsave(joinpath(outfolder, "J1_bare.jld2");
        interaction=hr, max_imag, max_ratio,
        head_mpc=hr_head_mpc)
end
function _postprocessing_wieps_exchange1(::Val{:jld2}, outfolder, predata)
    hrs = _int_make_hrs_nwnwnRnω(predata.exchange1_value, predata.exchange_R, ones(Int, length(predata.exchange_R)))
    hrs_head_mpc = _int_make_hrs_nwnwnRnω(predata.exchange1_value_head_mpc, predata.exchange_R, ones(Int, length(predata.exchange_R)))
    hrs_wing_mpc = _int_make_hrs_nwnwnRnω(predata.exchange1_value_wing_mpc, predata.exchange_R, ones(Int, length(predata.exchange_R)))
    for iω in eachindex(hrs)
        max_imag = maximum(abs∘imag, hrs[iω].vals)
        max_ratio = maximum(x->abs(imag(x)/real(x)), hrs[iω].vals)
        jldsave(joinpath(outfolder, "J1_$iω.jld2");
            interaction=hrs[iω], max_imag, max_ratio,
            head_mpc=hrs_head_mpc[iω],
            wing_mpc=hrs_wing_mpc[iω])
    end
end
