
include("./_witrunc.jl")
include("./_wotrunc.jl")

function int(::Val{3}, ρk::AbstractRhok{RhoT}, iqidx::Vector{Int};
	coulomb::PlaneWaveCoulomb,
	multipolecorrection::Union{AbstractMultipoleCorrection, Nothing} = nothing,
	outfolder::AbstractString = "./", fileformat::Symbol = :jld2,
	interaction_para,
) where {RhoT}

	# If no epsilon, only need cdot of vector.
	BLAS.set_num_threads(1)
	nq = ρk.nq

	predata = _preprocessing_int(ρk, interaction_para, multipolecorrection)
	prebuffer = _preprocessing_buffer(Val(:wotrunc), ρk, predata)
	ρR = RhoR(ρk, predata.ρRlist, predata.pair_ij;
		nchunk = interaction_para.JuliaTask_nchunk_ρR, ρRbuffer = interaction_para.ρR_buffer_level)

	start_time = time()
	for iq in iqidx
		print("iq = $(iq)/$(nq): ")
		flush(stdout)

		iq_start_time = time()
		ρR(prebuffer.ρR, iq)
		iq_ρR_time = time()
		_symmetrize_coulomb!(ρk, prebuffer.ρR, coulomb, iq)
		iq_coulomb_time = time()

		_wotrunc_int(ρk, coulomb, predata, prebuffer.ρR, iq)

		ρR_elapsed = iq_ρR_time - iq_start_time
		coulomb_elapsed = iq_coulomb_time - iq_ρR_time
		elapsed = time() - iq_start_time
		println("$(round(elapsed, digits=1))s; ",
			"ρR costs $(round(ρR_elapsed, digits=1))s, ",
			"coulomb costs $(round(coulomb_elapsed, digits=1))s.")
	end
	elapsed = time() - start_time
	println("alltime: $(round(elapsed, digits=1))s")
	flush(stdout)

	_postprocessing_woeps_int(fileformat, ρk, outfolder, predata)

	return nothing
end
function int(::Union{Val{2}, Val{1}}, ρk::AbstractRhok{RhoT}, iqidx::Vector{Int};
	coulomb::PlaneWaveCoulomb,
	multipolecorrection::Nothing = nothing,
	outfolder::AbstractString = "./", fileformat::Symbol = :jld2,
	interaction_para,
) where {RhoT}

	BLAS.set_num_threads(1)
	nq = ρk.nq

	predata = _preprocessing_int(ρk, interaction_para, multipolecorrection)
	prebuffer = _preprocessing_buffer(Val(:witrunc), ρk, predata)
	ρR = RhoR(ρk, predata.ρRlist, predata.pair_ij;
		nchunk = interaction_para.JuliaTask_nchunk_ρR, ρRbuffer = interaction_para.ρR_buffer_level)

	start_time = time()
	for iq in iqidx
		print("iq = $(iq)/$(nq): ")
		flush(stdout)

		iq_start_time = time()
		ρR(prebuffer.ρR, iq)
		iq_ρR_time = time()
		_symmetrize_coulomb!(ρk, prebuffer.ρR, coulomb, iq)
		iq_coulomb_time = time()

		_witrunc_int(ρk, coulomb, predata, prebuffer.ρR, iq)

		ρR_elapsed = iq_ρR_time - iq_start_time
		coulomb_elapsed = iq_coulomb_time - iq_ρR_time
		elapsed = time() - iq_start_time
		println("$(round(elapsed, digits=1))s; ",
			"ρR costs $(round(ρR_elapsed, digits=1))s, ",
			"coulomb costs $(round(coulomb_elapsed, digits=1))s.")
	end
	elapsed = time() - start_time
	println("alltime: $(round(elapsed, digits=1))s")
	flush(stdout)

	_postprocessing_woeps_int(fileformat, ρk, outfolder, predata)

	return nothing
end
function int(::Val{0}, ρk::Rhok_wospin{RhoT}, iqidx::Vector{Int};
	coulomb::PlaneWaveCoulomb,
	multipolecorrection::Nothing = nothing,
	outfolder::AbstractString = "./", fileformat::Symbol = :jld2,
	interaction_para,
) where {RhoT}

	error("Need to code!")

	# need to test
	FFTW.set_num_threads(Threads.nthreads())


	return nothing
end
