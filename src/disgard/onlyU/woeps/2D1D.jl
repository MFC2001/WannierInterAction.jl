function intU(::Union{Val{2}, Val{1}}, ρk::AbstractRhok{RhoT},
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

	prebuffer = _preprocessing_buffer_onlyU(Val(2), ρk, Threads.nthreads())
	ρR = RhoR(ρk; nchunk = interaction_para.JuliaTask_nchunk_ρR)

	start_time = time()
	for iq in iqidx
		print("iq = $(iq)/$(nq): ")
		flush(stdout)

		iq_start_time = time()
		ρR(prebuffer.ρR, iq)
		iq_ρR_time = time()
		truncatefactor_iq = truncatefactor(iq)
		iq_trunc_time = time()

		nG_rho_iq = nG_rho[iq]

		if iq == q0_idx
			q = q0
			G0idx_q0 = ρk.G0idx_q0
			truncatefactor_iq[G0idx_q0] = 1.0
			_witrunc_woeps_q0_onlyU(predata, prebuffer.ρR, truncatefactor_iq, prebuffer.rho, nG_rho_iq, G0idx_q0, head_value)
		else
			q = qgrid[iq]
			_witrunc_woeps_not_q0_onlyU(predata, q, prebuffer.ρR, truncatefactor_iq, prebuffer.rho, nG_rho_iq)
		end

		ρR_elapsed = iq_ρR_time - iq_start_time
		trunc_elapsed = iq_trunc_time - iq_ρR_time
		elapsed = time() - iq_start_time
		println("$(round(elapsed, digits=1))s; ",
			"ρR costs $(round(ρR_elapsed, digits=1))s, ",
			"truncate_factor costs $(round(trunc_elapsed, digits=1))s.")
	end
	elapsed = time() - start_time
	println("alltime: $(round(elapsed, digits=1))s")
	flush(stdout)

	_postprocessing_onlyU(fileformat, ρk, outfolder, predata, head_value)

	return nothing
end
function _witrunc_woeps_q0_onlyU(predata, ρR_data, truncatefactor_iq, buffer_rho, nG_rho_iq, G0idx_q0, head)
	head_correction = head - 1
	Threads.@threads :static for ichunk in Base.OneTo(predata.nchunk)
		@inbounds begin
			buffer_rho_chunk = view(buffer_rho[ichunk], 1:nG_rho_iq)
			for (iw1, iw2) in predata.ijji_idx[ichunk]
				ρR_left_G0 = ρR_data[G0idx_q0, iw1]
				ρR_left = view(ρR_data, 1:nG_rho_iq, iw1)
				broadcast!(*, buffer_rho_chunk, view(ρR_data, 1:nG_rho_iq, iw2), truncatefactor_iq)
				ρR_right_G0 = buffer_rho_chunk[G0idx_q0]
				head_correction_idx = head_correction * conj(ρR_left_G0) * ρR_right_G0
				predata.ijji_value_iq[iw1, iw2] = ρR_left ⋅ buffer_rho_chunk + head_correction_idx
			end
		end
	end
	Threads.@threads :static for iR in eachindex(predata.ijji_R)
		@inbounds begin
			ijji_value_iR = view(predata.ijji_value,:,:,iR)
			ijji_value_iR .+= predata.ijji_value_iq
		end
	end
end
function _witrunc_woeps_not_q0_onlyU(predata, q, ρR_data, truncatefactor_iq, buffer_rho, nG_rho_iq)
	Threads.@threads :static for ichunk in Base.OneTo(predata.nchunk)
		@inbounds begin
			buffer_rho_chunk = view(buffer_rho[ichunk], 1:nG_rho_iq)
			for (iw1, iw2) in predata.ijji_idx[ichunk]
				ρR_left = view(ρR_data, 1:nG_rho_iq, iw1)
				broadcast!(*, buffer_rho_chunk, view(ρR_data, 1:nG_rho_iq, iw2), truncatefactor_iq)
				predata.ijji_value_iq[iw1, iw2] = ρR_left ⋅ buffer_rho_chunk
			end
		end
	end
	Threads.@threads :static for iR in eachindex(ijji_R)
		@inbounds begin
			ijji_value_iR = view(predata.ijji_value,:,:,iR)
			phase = cispi(-2 * (q ⋅ ijji_R[iR]))
			ijji_value_iR .+= predata.ijji_value_iq .* phase
		end
	end
end
