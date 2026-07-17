function intU(::Val{3}, ρk::AbstractRhok{RhoT},
	head::AbstractHead, truncatefactor::TruncateFactor, iqidx::Vector{Int};
	outfolder::AbstractString = "./", fileformat::Symbol = :jld2,
	interaction_para,
) where {RhoT}

	nq = ρk.nq
	nG_rho = ρk.nG
	qgrid = ρk.qgrid
	kgrid = ρk.kgrid
	q0 = ρk.q0
	q0_idx = ρk.q0_idx

	# head term
	head_value = head(Val(:wocoff))

	predata = _preprocessing_int_onlyU(ρk, 0, Threads.nthreads())

	prebuffer = _preprocessing_buffer_onlyU(Val(3), ρk, Threads.nthreads())
	ρR = RhoR(ρk; nchunk = interaction_para.JuliaTask_nchunk_ρR)

	start_time = time()
	for iq in iqidx
		print("iq = $(iq)/$(nq): ")
		flush(stdout)

		iq_start_time = time()
		ρR(prebuffer.ρR, iq)
		iq_ρR_time = time()

		nG_rho_iq = nG_rho[iq]
		if iq == q0_idx
			q = q0
			G0idx_q0 = ρk.G0idx_q0
			_wotrunc_woeps_q0_onlyU(predata, prebuffer.ρR, nG_rho_iq, G0idx_q0, head_value)
		else
			q = qgrid[iq]
			_wotrunc_woeps_not_q0_onlyU(predata, q, prebuffer.ρR, nG_rho_iq)
		end

		ρR_elapsed = iq_ρR_time - iq_start_time
		elapsed = time() - iq_start_time
		println("$(round(elapsed, digits=1))s, ",
			"ρR costs $(round(ρR_elapsed, digits=1))s.")
	end
	elapsed = time() - start_time
	println("alltime: $(round(elapsed, digits=1))s")
	flush(stdout)

	_postprocessing_onlyU(fileformat, ρk, outfolder, predata, head_value)

	return nothing
end
function _wotrunc_woeps_q0_onlyU(predata, ρR_data, nG_rho_iq, G0idx_q0, head)
	head_correction = head - 1
	Threads.@threads :static for ichunk in Base.OneTo(predata.nchunk)
		@inbounds for (iw1, iw2) in predata.ijji_idx[ichunk]
			ρR_left_G0 = ρR_data[G0idx_q0, iw1]
			ρR_left = view(ρR_data, 1:nG_rho_iq, iw1)
			ρR_right_G0 = ρR_data[G0idx_q0, iw2]
			ρR_right = view(ρR_data, 1:nG_rho_iq, iw2)
			head_correction_idx = head_correction * conj(ρR_left_G0) * ρR_right_G0
			predata.ijji_value_iq[iw1, iw2] = ρR_left ⋅ ρR_right + head_correction_idx
		end
	end
	Threads.@threads :static for iR in eachindex(predata.ijji_R)
		@inbounds begin
			ijji_value_iR = view(predata.ijji_value,:,:,iR)
			ijji_value_iR .+= predata.ijji_value_iq
		end
	end
end
function _wotrunc_woeps_not_q0_onlyU(predata, q, ρR_data, nG_rho_iq)
	Threads.@threads :static for ichunk in Base.OneTo(predata.nchunk)
		@inbounds for (iw1, iw2) in predata.ijji_idx[ichunk]
			ρR_left = view(ρR_data, 1:nG_rho_iq, iw1)
			ρR_right = view(ρR_data, 1:nG_rho_iq, iw2)
			predata.ijji_value_iq[iw1, iw2] = ρR_left ⋅ ρR_right
		end
	end
	wait.(tasks)
	Threads.@threads :static for iR in eachindex(predata.ijji_R)
		@inbounds begin
			ijji_value_iR = view(predata.ijji_value,:,:,iR)
			phase = cispi(-2 * (q ⋅ predata.ijji_R[iR]))
			ijji_value_iR .+= predata.ijji_value_iq .* phase
		end
	end
end
