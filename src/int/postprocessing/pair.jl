function _postprocessing_woeps_pair(::Val{:txt}, outfolder, predata)
	pairint = _postprocessing_pair(Val(predata.pair_hermitian_part), predata)
	write(joinpath(outfolder, "pair_bare.dat"), pairint, WanIntijklR)
end
function _postprocessing_wieps_pair(::Val{:txt}, outfolder, predata)
	nω = size(predata.pair_value, 4)
	for iω in Base.OneTo(nω)
		pairint = _postprocessing_pair(Val(predata.pair_hermitian_part), predata, iω)
		write(joinpath(outfolder, "pair_$iω.dat"), pairint, WanIntijklR)
	end
end
function _postprocessing_woeps_pair(::Val{:jld2}, outfolder, predata)
	pairint = _postprocessing_pair(Val(predata.pair_hermitian_part), predata)
	jldsave(joinpath(outfolder, "pair_bare.jld2");
		interaction = pairint,
		head_mpc = predata.pair_value_bare_head_mpc)
end
function _postprocessing_wieps_pair(::Val{:jld2}, outfolder, predata)
	nω = size(predata.pair_value, 4)
	for iω in Base.OneTo(nω)
		pairint = _postprocessing_pair(Val(predata.pair_hermitian_part), predata, iω)
		jldsave(joinpath(outfolder, "pair_$iω.jld2");
			interaction = pairint,
			head_mpc = predata.pair_value_head_mpc[:, :, iω],
			wing_mpc = predata.pair_value_wing_mpc[:, :, iω])
	end
end
_postprocessing_pair(::Val{:all}, predata) =
	WanIntijklR(predata.pair_value_bare,
		predata.pair, predata.pair_pR, predata.pair_R;
		nw = predata.nw, atol = predata.pair_atol,
		dosort = true, dofilter = true)
_postprocessing_pair(::uplo, predata) where {uplo <: Union{Val{:U}, Val{:L}}} =
	WanIntijklR(predata.pair_value_bare,
		predata.pair_uploidx, predata.pair, predata.pair_pR, predata.pair_R;
		nw = predata.nw, atol = predata.pair_atol,
		uplo = predata.pair_hermitian_part, ishermitian = true, dofilter = true,
	)
_postprocessing_pair(::Val{:all}, predata, iω) =
	WanIntijklR(view(predata.pair_value,:,:,:,iω),
		predata.pair, predata.pair_pR, predata.pair_R;
		nw = predata.nw, atol = predata.pair_atol,
		dosort = true, dofilter = true)
_postprocessing_pair(::uplo, predata, iω) where {uplo <: Union{Val{:U}, Val{:L}}} =
	WanIntijklR(view(predata.pair_value,:,:,:,iω),
		predata.pair_uploidx, predata.pair, predata.pair_pR, predata.pair_R;
		nw = predata.nw, atol = predata.pair_atol,
		uplo = predata.pair_hermitian_part, ishermitian = true, dofilter = true,
	)



# _postprocessing_pair!(::Val{:U}, pair_value::Array, pair_U_idx, pair_L_idx, pair, pair_R, pair_pR) =
# 	_postprocessing_pair!(pair_value, pair_U_idx, pair_L_idx, pair, pair_R, pair_pR)
# _postprocessing_pair!(::Val{:L}, pair_value::Array, pair_U_idx, pair_L_idx, pair, pair_R, pair_pR) =
# 	_postprocessing_pair!(pair_value, pair_L_idx, pair_U_idx, pair, pair_R, pair_pR)
# _postprocessing_pair!(::Val{:all}, pair_value::Array, pair_U_idx, pair_L_idx, pair, pair_R, pair_pR) = pair_value
# function _postprocessing_pair!(pair_value::Array{T, 3}, pair_cal_idx, pair_uncal_idx, pair, pair_R, pair_pR) where {T}

# 	pR_minus_pR = [R₁ - R₂ for R₁ in pair_pR, R₂ in pair_pR]
# 	R′ = unique(vec(pR_minus_pR))
# 	R′_to_idx = Dict{ReducedCoordinates{Int}, Int}(RR => i for (i, RR) in enumerate(R′))
# 	pR_minus_pR_to_R′ = map(RR->R′_to_idx[RR], pR_minus_pR)

# 	R_to_idx = Dict{ReducedCoordinates{Int}, Int}(RR => i for (i, RR) in enumerate(pair_R))
# 	nR = length(pair_R)
# 	nR′ = length(R′)
# 	R_add_R′_to_R = Matrix{Int}(undef, nR, nR′)
# 	for iR′ in Base.OneTo(nR′), iR in Base.OneTo(nR)
# 		R_add_R′ = pair_R[iR] + R′[iR′]
# 		R_add_R′_to_R[iR, iR′] = get(R_to_idx, R_add_R′, 0)
# 	end

# 	for (Ical, Iuncal) in zip(pair_cal_idx, pair_uncal_idx)
# 		(p₁_c, p₂_c) = Ical
# 		(p₁, p₂) = Iuncal
# 		iR′ = pR_minus_pR_to_R′[pair[p₂][3], pair[p₁][3]]
# 		for iR in Base.OneTo(nR)
# 			R_add_R′ = R_add_R′_to_R[iR, iR′]
# 			if !iszero(R_add_R′)
# 				pair_value[p₁, p₂, iR] = conj(pair_value[p₁_c, p₂_c, R_add_R′])
# 			end
# 		end
# 	end

# 	return pair_value
# end
# function _postprocessing_pair!(pair_value::Array{T, 4}, pair_cal_idx, pair_uncal_idx, pair, pair_R, pair_pR) where {T}

# 	pR_minus_pR = [R₁ - R₂ for R₁ in pair_pR, R₂ in pair_pR]
# 	R′ = unique(vec(pR_minus_pR))
# 	R′_to_idx = Dict{ReducedCoordinates{Int}, Int}(RR => i for (i, RR) in enumerate(R′))
# 	pR_minus_pR_to_R′ = map(RR->R′_to_idx[RR], pR_minus_pR)

# 	R_to_idx = Dict{ReducedCoordinates{Int}, Int}(RR => i for (i, RR) in enumerate(pair_R))
# 	nR = length(pair_R)
# 	nR′ = length(R′)
# 	R_add_R′_to_R = Matrix{Int}(undef, nR, nR′)
# 	for iR′ in Base.OneTo(nR′), iR in Base.OneTo(nR)
# 		R_add_R′ = pair_R[iR] + R′[iR′]
# 		R_add_R′_to_R[iR, iR′] = get(R_to_idx, R_add_R′, 0)
# 	end

# 	nω = size(pair_value, 1)
# 	for (Ical, Iuncal) in zip(pair_cal_idx, pair_uncal_idx)
# 		(p₁_c, p₂_c) = Ical
# 		(p₁, p₂) = Iuncal
# 		iR′ = pR_minus_pR_to_R′[pair[p₂][3], pair[p₁][3]]
# 		for iR in Base.OneTo(nR)
# 			R_add_R′ = R_add_R′_to_R[iR, iR′]
# 			if !iszero(R_add_R′)
# 				for iω in Base.OneTo(nω)
# 					pair_value[iω, p₁, p₂, iR] = conj(pair_value[iω, p₁_c, p₂_c, R_add_R′])
# 				end
# 			end
# 		end
# 	end

# 	return pair_value
# end
