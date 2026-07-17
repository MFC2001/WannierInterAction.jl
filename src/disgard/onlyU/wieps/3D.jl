function intU(::Val{3}, eps::DielectricPWBasis, ρk::AbstractRhok{RhoT},
	head::AbstractHead, truncatefactor::TruncateFactor, iqidx::Vector{Int};
	outfolder::AbstractString = "./", fileformat::Symbol = :jld2,
	interaction_para,
) where {RhoT}

	nq = ρk.nq
	nG_rho = ρk.nG
	nG_eps = eps.nG
	qgrid = ρk.qgrid
	kgrid = ρk.kgrid
	nω = eps.nω
	q0 = ρk.q0
	q0_idx = ρk.q0_idx
	have_negq0 = (eps.have_negq0 && ρk.have_negq0)

	# head, wing
	normq0 = sqrt(q0 ⋅ (ρk.rldot * q0))

	predata = _preprocessing_int_onlyU(ρk, nω, Threads.nthreads())

	prebuffer = _preprocessing_buffer_onlyU(Val(3), eps, ρk, Threads.nthreads())
	ρR = RhoR(ρk; nchunk = interaction_para.JuliaTask_nchunk_ρR)

	ϵ_∞ = Vector{ComplexF64}(undef, nω)
	head_value = Vector{ComplexF64}(undef, nω + 1)

	start_time = time()
	for iq in iqidx
		print("iq = $(iq)/$(nq): ")
		flush(stdout)

		iq_start_time = time()
		ρR(prebuffer.ρR, iq)
		iq_ρR_time = time()
		epsmat = eps[iq, :]
		iq_eps_time = time()

		nG_eps_iq = nG_eps[iq]
		nG_rho_iq = nG_rho[iq]
		if iq == q0_idx
			q = q0
			G0idx_q0 = ρk.G0idx_q0
			for iω in Base.OneTo(nω)
				ϵ_∞[iω] = epsmat[iω][G0idx_q0, G0idx_q0]
			end
			head_value .= head(Val(:wocoff), ϵ_∞)
			if have_negq0
				_wotrunc_wieps_q0_wi_negq0_onlyU(
					predata, q, epsmat, prebuffer.ρR,
					prebuffer.eps_mul, prebuffer.eps, nG_eps_iq, nG_rho_iq, G0idx_q0, head_value, normq0,
				)

				ρR(prebuffer.ρR, -1)
				epsmat = eps[-1, :]
				q = -q0
				nG_eps_iq = eps.nG_negq0
				nG_rho_iq = ρk.nG_negq0
				G0idx_q0 = ρk.G0idx_negq0

				_wotrunc_wieps_q0_wi_negq0_onlyU(
					predata, q, epsmat, prebuffer.ρR,
					prebuffer.eps_mul, prebuffer.eps, nG_eps_iq, nG_rho_iq, G0idx_q0, head_value, normq0,
				)
			else
				_wotrunc_wieps_q0_wo_negq0_onlyU(
					predata, q, epsmat, prebuffer.ρR,
					prebuffer.eps_mul, prebuffer.eps, nG_eps_iq, nG_rho_iq, G0idx_q0, head_value,
				)
			end
		else
			q = qgrid[iq]
			_wotrunc_wieps_not_q0_onlyU(
				predata, q, epsmat, prebuffer.ρR,
				prebuffer.eps_mul, nG_eps_iq, nG_rho_iq,
			)
		end

		ρR_elapsed = iq_ρR_time - iq_start_time
		eps_elapsed = iq_eps_time - iq_ρR_time
		elapsed = time() - iq_start_time
		println("$(round(elapsed, digits=1))s; ",
			"ρR costs $(round(ρR_elapsed, digits=1))s, ",
			"epsmat costs $(round(eps_elapsed, digits=1))s.")
	end
	elapsed = time() - start_time
	println("alltime: $(round(elapsed, digits=1))s")
	flush(stdout)

	_postprocessing_onlyU(fileformat, ρk, outfolder, predata, head_value, ϵ_∞)

	return nothing
end
function _wotrunc_wieps_q0_wi_negq0_onlyU(
	predata, q, epsmat, ρR_data,
	buffer_eps_mul, buffer_eps, nG_eps_iq, nG_rho_iq, G0idx_q0, head, normq0,
)
	head_correction = head .- 1
	wing_correction = 1 / normq0 - 1
	epsG0G = map(ϵ -> view(ϵ, G0idx_q0, :), epsmat)
	epsGG0 = map(ϵ -> view(ϵ, :, G0idx_q0), epsmat)
	epsG0G0 = map(ϵ -> ϵ[G0idx_q0, G0idx_q0], epsmat)
	Threads.@threads :static for ichunk in Base.OneTo(predata.nchunk)
		@inbounds begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[ichunk], 1:nG_eps_iq)
			buffer_eps_chunk = view(buffer_eps[ichunk], 1:nG_eps_iq)
			for (iw1, iw2) in predata.ijji_idx[ichunk]
				ρR_left_G0 = ρR_data[G0idx_q0, iw1]
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, iw1)
				ρR_left_high = view(ρR_data, (nG_eps_iq+1):nG_rho_iq, iw1)
				ρR_right_G0 = ρR_data[G0idx_q0, iw2]
				ρR_right_eps = view(ρR_data, 1:nG_eps_iq, iw2)
				ρR_right_high = view(ρR_data, (nG_eps_iq+1):nG_rho_iq, iw2)
				high_dot = ρR_left_high ⋅ ρR_right_high
				ρR_G0_G0 = conj(ρR_left_G0) * ρR_right_G0
				for iω in eachindex(epsmat)
					head_term = epsG0G0[iω] * ρR_G0_G0
					buffer_eps_chunk .= epsG0G[iω] .* ρR_right_eps
					wing1 = conj(ρR_left_G0) * sum(buffer_eps_chunk) - head_term
					wing2 = (ρR_left_eps ⋅ epsGG0[iω]) * ρR_right_G0 - head_term
					wing_correction_idx = (wing1 + wing2) * wing_correction
					head_correction_idx = ρR_G0_G0 * head_correction[iω]

					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					predata.ijji_value_iq[iω, iw1, iw2] =
						ρR_left_eps ⋅ buffer_eps_mul_chunk + head_correction_idx + wing_correction_idx + high_dot
				end
				head_correction_idx = ρR_G0_G0 * head_correction[end]
				predata.ijji_value_iq[end, iw1, iw2] = ρR_left_eps ⋅ ρR_right_eps + head_correction_idx + high_dot
			end
		end
	end
	Threads.@threads :static for iR in eachindex(predata.ijji_R)
		@inbounds begin
			ijji_value_iR = view(predata.ijji_value,:,:,:,iR)
			# phase = cispi(-2 * (q ⋅ predata.ijji_R[iR])) / 2
			# ijji_value_iR .+= predata.ijji_value_iq .* phase
			ijji_value_iR .+= predata.ijji_value_iq ./ 2
		end
	end
end
function _wotrunc_wieps_q0_wo_negq0_onlyU(
	predata, q, epsmat, ρR_data,
	buffer_eps_mul, buffer_eps, nG_eps_iq, nG_rho_iq, G0idx_q0, head,
)
	head_correction = head .- 1
	epsG0G = map(ϵ -> view(ϵ, G0idx_q0, :), epsmat)
	epsGG0 = map(ϵ -> view(ϵ, :, G0idx_q0), epsmat)
	epsG0G0 = map(ϵ -> ϵ[G0idx_q0, G0idx_q0], epsmat)
	Threads.@threads :static for ichunk in Base.OneTo(predata.nchunk)
		@inbounds begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[ichunk], 1:nG_eps_iq)
			buffer_eps_chunk = view(buffer_eps[ichunk], 1:nG_eps_iq)
			for (iw1, iw2) in predata.ijji_idx[ichunk]
				ρR_left_G0 = ρR_data[G0idx_q0, iw1]
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, iw1)
				ρR_left_high = view(ρR_data, (nG_eps_iq+1):nG_rho_iq, iw1)
				ρR_right_G0 = ρR_data[G0idx_q0, iw2]
				ρR_right_eps = view(ρR_data, 1:nG_eps_iq, iw2)
				ρR_right_high = view(ρR_data, (nG_eps_iq+1):nG_rho_iq, iw2)
				high_dot = ρR_left_high ⋅ ρR_right_high
				ρR_G0_G0 = conj(ρR_left_G0) * ρR_right_G0
				for iω in eachindex(epsmat)
					head_term = epsG0G0[iω] * ρR_G0_G0
					buffer_eps_chunk .= epsG0G[iω] .* ρR_right_eps
					wing1 = conj(ρR_left_G0) * sum(buffer_eps_chunk) - head_term
					wing2 = (ρR_left_eps ⋅ epsGG0[iω]) * ρR_right_G0 - head_term
					wing_idx = wing1 + wing2
					head_correction_idx = ρR_G0_G0 * head_correction[iω]

					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					predata.ijji_value_iq[iω, iw1, iw2] =
						ρR_left_eps ⋅ buffer_eps_mul_chunk + head_correction_idx - wing_idx + high_dot
				end
				head_correction_idx = ρR_G0_G0 * head_correction[end]
				predata.ijji_value_iq[end, iw1, iw2] = ρR_left_eps ⋅ ρR_right_eps + head_correction_idx + high_dot
			end
		end
	end
	Threads.@threads :static for iR in eachindex(predata.ijji_R)
		@inbounds begin
			ijji_value_iR = view(predata.ijji_value,:,:,:,iR)
			# phase = cispi(-2 * (q ⋅ predata.ijji_R[iR]))
			# ijji_value_iR .+= predata.ijji_value_iq .* phase
			ijji_value_iR .+= predata.ijji_value_iq
		end
	end
end
function _wotrunc_wieps_not_q0_onlyU(
	predata, q, epsmat, ρR_data,
	buffer_eps_mul, nG_eps_iq, nG_rho_iq,
)
	Threads.@threads :static for ichunk in Base.OneTo(predata.nchunk)
		@inbounds begin
			buffer_eps_mul_chunk = view(buffer_eps_mul[ichunk], 1:nG_eps_iq)
			for (iw1, iw2) in predata.ijji_idx[ichunk]
				ρR_left_eps = view(ρR_data, 1:nG_eps_iq, iw1)
				ρR_left_high = view(ρR_data, (nG_eps_iq+1):nG_rho_iq, iw1)
				ρR_right_eps = view(ρR_data, 1:nG_eps_iq, iw2)
				ρR_right_high = view(ρR_data, (nG_eps_iq+1):nG_rho_iq, iw2)
				high_dot = ρR_left_high ⋅ ρR_right_high
				for iω in eachindex(epsmat)
					mul!(buffer_eps_mul_chunk, epsmat[iω], ρR_right_eps)
					predata.ijji_value_iq[iω, iw1, iw2] = ρR_left_eps ⋅ buffer_eps_mul_chunk + high_dot
				end
				predata.ijji_value_iq[end, iw1, iw2] = ρR_left_eps ⋅ ρR_right_eps + high_dot
			end
		end
	end
	Threads.@threads :static for iR in eachindex(predata.ijji_R)
		@inbounds begin
			ijji_value_iR = view(predata.ijji_value,:,:,:,iR)
			phase = cispi(-2 * (q ⋅ predata.ijji_R[iR]))
			ijji_value_iR .+= predata.ijji_value_iq .* phase
		end
	end
end
