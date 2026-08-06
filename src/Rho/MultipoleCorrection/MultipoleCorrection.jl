abstract type AbstractMultipoleCorrection end
# 统一存储 Ak 和 Bk
# 在int内获取R变为AR 和 BR

struct MultipoleCorrectionDirectional <: AbstractMultipoleCorrection
	q̂::CartesianCoordinates{Float64}
	berry::Array{ComplexF64} # [3, nk, nw, nw]
	# 取特定方向的结果时，难以计算二阶导内积。
end


export MultipoleCorrection

struct MultipoleCorrectionAveraged <: AbstractMultipoleCorrection
	kpoints::Vector{ReducedCoordinates{Float64}}
	berry::Array{ComplexF64} # [3, nk, nw, nw]
	metric::Array{ComplexF64} # [nk, nw, nw]
end
function MultipoleCorrection(
	M::AbstractArray{<:Number, 4}, # [nw, nw, nntot, nk]
	nnkp::wannier90_nnkp;
	shell_tol::Real = 1e-6, orth_tol::Real = 1e-6,
	kwargs...,
)

	nk = size(M, 4)
	@assert length(nnkp.kpoints) == nk "Mismatched M and nnkp!"
	nntot = size(M, 3)
	@assert size(nnkp.nnkpts, 2) == nk * nntot "Mismatched M and nnkp!"

	# note: 这里遵循wannier90的逻辑，先ib后ik。
	# write (nnkpout, '(a)') 'begin nnkpts'
	# write (nnkpout, '(i4)') nntot
	# do nkp = 1, num_kpts
	#   do nn = 1, nntot
	#     write (nnkpout, '(2i6,3x,3i4)') &
	#       nkp, nnlist(nkp, nn), (nncell(i, nkp, nn), i=1, 3)
	#   end do
	# end do
	# write (nnkpout, '(a/)') 'end nnkpts'
	blists = Matrix{CartesianCoordinates{Float64}}(undef, nntot, nk)
	ikb = 0
	for ik in 1:nk, ib in 1:nntot
		ikb += 1
		kb = nnkp.nnkpts[:, ikb]
		G = ReducedCoordinates{Int}(kb[3], kb[4], kb[5])
		b_frac = nnkp.kpoints[kb[2]] + G - nnkp.kpoints[kb[1]]
		blists[ib, ik] = nnkp.rlattice * b_frac
	end
	weights = Matrix{Float64}(undef, nntot, nk)
	for ik in 1:nk
		# 这个计算很轻量，一般在1s内。
		weights[:, ik] = WannierInterActionBase._blist_to_shell_weights(blists[:, ik]; shell_tol, orth_tol)
	end

	return MultipoleCorrection(M, nnkp.kpoints, blists, weights; kwargs...)
end
function MultipoleCorrection(
	M::AbstractArray{<:Number, 4}, # [nw, nw, nntot, nk]
	kpoints::AbstractVector{<:AbstractVector{<:Real}},
	blists::AbstractMatrix{<:AbstractVector{<:Real}}, # [nntot, nk]
	weights::AbstractMatrix{<:Real}; # [nntot, nk]
	q̂ = nothing,
)

	nntot = size(M, 3)
	nk = size(M, 4)
	@assert length(kpoints) == nk "Mismatched M and kpoints!"
	@assert size(blists) == (nntot, nk) "Mismatched M and blists!"
	@assert size(weights) == (nntot, nk) "Mismatched M and weights!"

	kpoints = ReducedCoordinates{Float64}.(kpoints)
	blists = CartesianCoordinates{Float64}.(blists)

	nw = size(M, 1)
	berry = Array{ComplexF64}(undef, 3, nk, nw, nw)
	metric = Array{ComplexF64}(undef, nk, nw, nw)

	for j in 1:nw, i in 1:nw
		δᵢⱼ = i == j ? 1 : 0
		for ik in 1:nk
			berry[:, ik, i, j] .= zero(ComplexF64)
			metric[ik, i, j] = zero(ComplexF64)
			for ib in 1:nntot
				berry[:, ik, i, j] += weights[ib, ik] * (M[i, j, ib, ik] - δᵢⱼ) * blists[ib, ik]
				metric[ik, i, j] += weights[ib, ik] * (2 * δᵢⱼ - M[i, j, ib, ik] - conj(M[j, i, ib, ik]))
			end
		end
	end

	return MultipoleCorrectionAveraged(kpoints, berry, metric)
end
struct MultipoleCorrectionAveraged_coll <: AbstractMultipoleCorrection
	kpoints::Vector{ReducedCoordinates{Float64}}
	berry_up::Array{ComplexF64}
	metric_up::Array{ComplexF64}
	berry_dn::Array{ComplexF64}
	metric_dn::Array{ComplexF64}
	nw_up::Int
	nw_dn::Int
end
function MultipoleCorrection(M_up::MultipoleCorrectionAveraged, M_dn::MultipoleCorrectionAveraged)
	nw_up = size(M_up.berry, 3)
	nw_dn = size(M_dn.berry, 3)
	return MultipoleCorrectionAveraged_coll(M_up.kpoints, M_up.berry, M_up.metric, M_dn.berry, M_dn.metric, nw_up, nw_dn)
end

function (M::MultipoleCorrectionAveraged)(
	lattice::Union{AbstractMatrix{<:Real}, Lattice},
	ρRlist::AbstractVector{<:AbstractVector{<:Integer}},
	pair_ij::AbstractVector{<:NTuple{2, <:Integer}},
)
	nR = length(ρRlist)
	nij = length(pair_ij)

	ρRlist = ReducedCoordinates{Int}.(ρRlist)
	lattice = Lattice(parent(lattice))
	ρRlist_car = map(x->lattice * x, ρRlist)

	AR = Array{ComplexF64}(undef, 3, nR, nij)
	BR = Array{ComplexF64}(undef, nR, nij)
	CR = Array{ComplexF64}(undef, nR, nij)

	nk = length(M.kpoints)
	ekρR = [cispi(2 * (k ⋅ R)) / nk for k in M.kpoints, R in ρRlist]
	for ij in 1:nij
		(i, j) = pair_ij[ij]
		Ak_ij = view(M.berry,:,:,i,j)
		AR_ij = view(AR,:,:,ij)
		mul!(AR_ij, Ak_ij, ekρR)
	end

	ekρR = permutedims(ekρR, (2, 1))
	for ij in 1:nij
		(i, j) = pair_ij[ij]
		Bk_ij = view(M.metric, :, i, j)
		BR_ij = view(BR, :, ij)
		mul!(BR_ij, ekρR, Bk_ij)
	end

	for ij in 1:nij, iR in 1:nR
		R = ρRlist_car[iR]
		RA = R[1] * AR[1, iR, ij] + R[2] * AR[2, iR, ij] + R[3] * AR[3, iR, ij]
		CR[iR, ij] = complex(imag(RA), -real(RA)) - BR[iR, ij]
	end

	return AR, CR
end
function (M::MultipoleCorrectionAveraged_coll)(
	lattice::Union{AbstractMatrix{<:Real}, Lattice},
	ρRlist::AbstractVector{<:AbstractVector{<:Integer}},
	pair_ij::AbstractVector{<:NTuple{2, <:Integer}},
)
	nR = length(ρRlist)
	nij = length(pair_ij)

	ρRlist = ReducedCoordinates{Int}.(ρRlist)
	lattice = Lattice(parent(lattice))
	ρRlist_car = map(x->lattice * x, ρRlist)

	AR = Array{ComplexF64}(undef, 3, nR, nij)
	BR = Array{ComplexF64}(undef, nR, nij)
	CR = Array{ComplexF64}(undef, nR, nij)

	nk = length(M.kpoints)
	ekρR = [cispi(2 * (k ⋅ R)) / nk for k in M.kpoints, R in ρRlist]
	for ij in 1:nij
		(i, j) = pair_ij[ij]
		AR_ij = view(AR,:,:,ij)
		ispin = i <= M.nw_up ? :up : :down
		jspin = j <= M.nw_up ? :up : :down
		if ispin ≠ jspin
			fill!(AR_ij, zero(ComplexF64))
		elseif ispin == :up
			Ak_ij = view(M.berry_up,:,:,i,j)
			mul!(AR_ij, Ak_ij, ekρR)
		else
			Ak_ij = view(M.berry_dn,:,:,(i-M.nw_up),(j-M.nw_up))
			mul!(AR_ij, Ak_ij, ekρR)
		end
	end

	ekρR = permutedims(ekρR, (2, 1))
	for ij in 1:nij
		(i, j) = pair_ij[ij]
		BR_ij = view(BR, :, ij)
		ispin = i <= M.nw_up ? :up : :down
		jspin = j <= M.nw_up ? :up : :down
		if ispin ≠ jspin
			fill!(BR_ij, zero(ComplexF64))
		elseif ispin == :up
			Bk_ij = view(M.metric_up, :, i, j)
			mul!(BR_ij, ekρR, Bk_ij)
		else
			Bk_ij = view(M.metric_dn,:,:,(i-M.nw_up),(j-M.nw_up))
			mul!(BR_ij, ekρR, Bk_ij)
		end
	end

	for ij in 1:nij, iR in 1:nR
		R = ρRlist_car[iR]
		RA = R[1] * AR[1, iR, ij] + R[2] * AR[2, iR, ij] + R[3] * AR[3, iR, ij]
		CR[iR, ij] = complex(imag(RA), -real(RA)) - BR[iR, ij]
	end

	return AR, CR
end
