
include("./_witrunc.jl")
include("./_wotrunc.jl")
include("./wing.jl")

function int(::Val{3}, eps::DielectricPWBasis, ρk::AbstractRhok{RhoT}, iqidx::Vector{Int};
	coulomb::PlaneWaveCoulomb,
	multipolecorrection::Union{AbstractMultipoleCorrection, Nothing} = nothing,
	outfolder::AbstractString = "./", fileformat::Symbol = :jld2,
	interaction_para,
) where {RhoT}

	BLAS.set_num_threads(1)
	nq = ρk.nq
	nω = eps.nω

	predata = _preprocessing_int(ρk, nω, interaction_para, multipolecorrection)
	prebuffer = _preprocessing_buffer(Val(:wotrunc), eps, ρk, predata)
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

		_wotrunc_int(ρk, eps, coulomb, predata, prebuffer.ρR, iq)

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

	_postprocessing_wieps_int(fileformat, ρk, outfolder, predata)

	return nothing
end
function int(::Union{Val{2}, Val{1}}, eps::DielectricPWBasis, ρk::AbstractRhok{RhoT}, iqidx::Vector{Int};
	coulomb::PlaneWaveCoulomb,
	multipolecorrection::Nothing = nothing,
	outfolder::AbstractString = "./", fileformat::Symbol = :jld2,
	interaction_para,
) where {RhoT}

	BLAS.set_num_threads(1)
	nq = ρk.nq
	nω = eps.nω

	predata = _preprocessing_int(ρk, nω, interaction_para, multipolecorrection)
	prebuffer = _preprocessing_buffer(Val(:witrunc), eps, ρk, predata)
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

		_witrunc_int(ρk, eps, coulomb, predata, prebuffer.ρR, iq)

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

	_postprocessing_wieps_int(fileformat, ρk, outfolder, predata)

	return nothing
end
