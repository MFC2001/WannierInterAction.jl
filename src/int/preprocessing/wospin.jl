function _preprocessing_int(ρk::Rhok_wospin{RhoT}, nω::Integer, interaction_para, multipolecorrection) where {RhoT}

	nw = ρk.nw
	qgrid = ρk.qgrid
	lattice = ρk.lattice
	Ω = ρk.Ω

	pair_ij = Set([(i, i) for i in 1:nw])
	ρRlist = Set([ReducedCoordinates{Int}(0, 0, 0)])

	if IJJI ∈ interaction_para.interaction
		ijji_R = gridindex(Val(:WignerSeitz), lattice, qgrid.kgrid_size)
		ijji_nR = length(ijji_R)
		ijji_idx = vec([(iw1, iw2) for iw1 in 1:nw, iw2 in 1:nw])
		ijji_value_bare = zeros(RhoT, nw, nw, ijji_nR)
		ijji_value_bare_iq = zeros(RhoT, nw, nw)
		ijji_value_bare_head = zeros(RhoT, nw, nw)
		ijji_value_bare_head_mpc = zeros(RhoT, nw, nw)
		ijji_value = iszero(nω) ? nothing : zeros(RhoT, nw, nw, ijji_nR, nω)
		ijji_value_iq = iszero(nω) ? nothing : zeros(RhoT, nw, nw)
		ijji_value_head = iszero(nω) ? nothing : zeros(RhoT, nw, nw, nω)
		ijji_value_head_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, nω)
		ijji_value_wing_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, nω)
		ϵ = iszero(nω) ? nothing : zeros(ComplexF64, nω)
	else
		ijji_R = nothing
		ijji_idx = nothing
		ijji_value_bare = nothing
		ijji_value_bare_iq = nothing
		ijji_value_bare_head = nothing
		ijji_value_bare_head_mpc = nothing
		ijji_value = nothing
		ijji_value_iq = nothing
		ijji_value_head = nothing
		ijji_value_head_mpc = nothing
		ijji_value_wing_mpc = nothing
		ϵ = nothing
	end
	if IJIJ ∈ interaction_para.interaction
		ijij_R = gridindex(interaction_para.ijij_range)
		ijij_nR = length(ijij_R)
		union!(ρRlist, ijij_R)
		pair_ij_add = vec([(i, j) for i in 1:nw, j in 1:nw])
		union!(pair_ij, pair_ij_add)
		ijij_R0idx = findfirst(R -> iszero(R), ijij_R)
		ijij_idx = [(iw1, iw2, iR) for iw1 in 1:nw, iw2 in 1:nw, iR in 1:ijij_nR]
		# without U-type interaction
		if !isnothing(ijij_R0idx)
			for iw in 1:nw
				ijij_idx[iw, iw, ijij_R0idx] = (0, 0, 0)
			end
			ijij_idx = vec(ijij_idx)
			filter!(idx -> !any(iszero, idx), ijij_idx)
		else
			ijij_idx = vec(ijij_idx)
		end
		ijij_value_bare = zeros(RhoT, nw, nw, ijij_nR)
		ijij_value_bare_iq = zeros(RhoT, nw, nw, ijij_nR)
		ijij_value_bare_head_mpc = zeros(RhoT, nw, nw, ijij_nR)
		ijij_value = iszero(nω) ? nothing : zeros(RhoT, nw, nw, ijij_nR, nω)
		ijij_value_head_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, ijij_nR, nω)
		ijij_value_wing_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, ijij_nR, nω)
	else
		ijij_R = nothing
		ijij_idx = nothing
		ijij_value_bare = nothing
		ijij_value_bare_iq = nothing
		ijij_value_bare_head_mpc = nothing
		ijij_value = nothing
		ijij_value_head_mpc = nothing
		ijij_value_wing_mpc = nothing
	end
	if IIJJ ∈ interaction_para.interaction
		iijj_R = gridindex(interaction_para.iijj_range)
		iijj_nR = length(iijj_R)
		union!(ρRlist, iijj_R)
		union!(ρRlist, map(R -> -R, iijj_R))
		pair_ij_add = vec([(i, j) for i in 1:nw, j in 1:nw])
		union!(pair_ij, pair_ij_add)
		iijj_R0idx = findfirst(R -> iszero(R), iijj_R)
		iijj_idx = [(iw1, iw2, iR) for iw1 in 1:nw, iw2 in 1:nw, iR in 1:iijj_nR]
		if !isnothing(iijj_R0idx)
			for iw in 1:nw
				iijj_idx[iw, iw, iijj_R0idx] = (0, 0, 0)
			end
			iijj_idx = vec(iijj_idx)
			filter!(idx -> !any(iszero, idx), iijj_idx)
		else
			iijj_idx = vec(iijj_idx)
		end
		iijj_value_bare = zeros(RhoT, nw, nw, iijj_nR)
		iijj_value_bare_iq = zeros(RhoT, nw, nw, iijj_nR)
		iijj_value_bare_head_mpc = zeros(RhoT, nw, nw, iijj_nR)
		iijj_value = iszero(nω) ? nothing : zeros(RhoT, nw, nw, iijj_nR, nω)
		iijj_value_head_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, iijj_nR, nω)
		iijj_value_wing_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, iijj_nR, nω)
	else
		iijj_R = nothing
		iijj_idx = nothing
		iijj_value_bare = nothing
		iijj_value_bare_iq = nothing
		iijj_value_bare_head_mpc = nothing
		iijj_value = nothing
		iijj_value_head_mpc = nothing
		iijj_value_wing_mpc = nothing
	end
	if PAIR ∈ interaction_para.interaction
		pair, pair_R, pair_pR, pair_idx =
			_preprocessing_pair(interaction_para.pair, interaction_para.pair_range, nw, qgrid.kgrid_size, lattice, interaction_para.pair_hermitian_part)

		pair_ij_add = map(pair) do p
			(p[1], p[2])
		end
		union!(pair_ij, pair_ij_add)
		union!(ρRlist, pair_pR)

		pair_isp0 = map(p -> WannierInterActionBase._iszero_pair(p, pair_pR), pair)

		pair_uploidx = pair_idx
		pair_nR = length(pair_R)
		pair_np = length(pair)

		println("pair_np: ", pair_np)
		println("pair_nR: ", pair_nR)

		pair_value_bare = zeros(RhoT, pair_np, pair_np, pair_nR)
		pair_value_bare_iq = zeros(RhoT, pair_np, pair_np)
		pair_value_bare_head_mpc = zeros(RhoT, pair_np, pair_np)
		pair_value = iszero(nω) ? nothing : zeros(RhoT, pair_np, pair_np, pair_nR, nω)
		pair_value_iq = iszero(nω) ? nothing : zeros(RhoT, pair_np, pair_np)
		pair_value_head_mpc = iszero(nω) ? nothing : zeros(RhoT, pair_np, pair_np, nω)
		pair_value_wing_mpc = iszero(nω) ? nothing : zeros(RhoT, pair_np, pair_np, nω)
	else
		pair = nothing
		pair_R = nothing
		pair_pR = nothing
		pair_idx = nothing
		pair_isp0 = nothing
		pair_uploidx = nothing
		pair_value_bare = nothing
		pair_value_bare_iq = nothing
		pair_value_bare_head_mpc = nothing
		pair_value = nothing
		pair_value_iq = nothing
		pair_value_head_mpc = nothing
		pair_value_wing_mpc = nothing
	end

	pair_ij = collect(pair_ij)
	ρRlist = collect(ρRlist)

	pair_ij_to_idx = Dict{NTuple{2, Int}, Int}(ij => idx for (idx, ij) in enumerate(pair_ij))
	ρRlist_to_idx = Dict{ReducedCoordinates{Int}, Int}(R => i for (i, R) in enumerate(ρRlist))

	JuliaTask_nchunk_int = interaction_para.JuliaTask_nchunk_int
	# _idx分块。后续进行内存预分配需要预分块，如果不涉及大量内存分配，建议使用Julia原生多线程循环，会自动分块。
	# 目前没有任务依赖的内存分配，但是为保持拓展性，暂时保持手动任务分块。
	ρR0idx = ρRlist_to_idx[ReducedCoordinates{Int}(0, 0, 0)]
	if IJJI ∈ interaction_para.interaction
		ijji_idx = map(ijji_idx) do (i, j)
			ii = pair_ij_to_idx[(i, i)]
			jj = pair_ij_to_idx[(j, j)]
			return (i, j, ρR0idx, ii, ρR0idx, jj)
		end
		ijji_idx_chunk = _split_tasks(ijji_idx, JuliaTask_nchunk_int)
	else
		ijji_idx = nothing
		ijji_idx_chunk = nothing
	end
	if IJIJ ∈ interaction_para.interaction
		ijij_iR2ρR = map(R -> ρRlist_to_idx[R], ijij_R)
		ijij_idx = map(ijij_idx) do (i, j, iR)
			ijij_ρR = ijij_iR2ρR[iR]
			ji = pair_ij_to_idx[(j, i)]
			return (i, j, iR, ijij_ρR, ji, ijij_ρR, ji)
		end
		ijij_idx_chunk = _split_tasks(ijij_idx, JuliaTask_nchunk_int)
	else
		ijij_idx = nothing
		ijij_idx_chunk = nothing
	end
	if IIJJ ∈ interaction_para.interaction
		iijj_iR2ρR_left = map(R -> ρRlist_to_idx[R], iijj_R)
		iijj_iR2ρR_right = map(R -> ρRlist_to_idx[-R], iijj_R)
		iijj_idx = map(iijj_idx) do (i, j, iR)
			ρR_left = iijj_iR2ρR_left[iR]
			ρR_right = iijj_iR2ρR_right[iR]
			ji = pair_ij_to_idx[(j, i)]
			ij = pair_ij_to_idx[(i, j)]
			return (i, j, iR, ρR_left, ji, ρR_right, ij)
		end
		iijj_idx_chunk = _split_tasks(iijj_idx, JuliaTask_nchunk_int)
	else
		iijj_idx = nothing
		iijj_idx_chunk = nothing
	end
	if PAIR ∈ interaction_para.interaction
		pair_idx_ρR = map(pair) do (i, j, ipR)
			ji = pair_ij_to_idx[(j, i)] # note!
			iρR = ρRlist_to_idx[pair_pR[ipR]]
			return (iρR, ji)
		end
		pair_idx = map(pair_idx) do (ip1, ip2)
			(ρR_left, li) = pair_idx_ρR[ip1]
			(ρR_right, jk) = pair_idx_ρR[ip2]
			return (ip1, ip2, ρR_left, li, ρR_right, jk)
		end
		pair_idx_chunk = _split_tasks(pair_idx, JuliaTask_nchunk_int)
	else
		pair_idx = nothing
		pair_idx_chunk = nothing
	end

	# multipolecorrection
	if isnothing(multipolecorrection)
		havempc = false
		AR = nothing
		CR = nothing
		head_mpc_coff = 0
		wing_mpc_coff = 0
	else
		havempc = true
		AR, CR = multipolecorrection(lattice, ρRlist, pair_ij)
		head_mpc_coff = CoulombScale * 4π / Ω
		wing_mpc_coff = sqrt(head_mpc_coff)
	end

	return (interaction = interaction_para.interaction, nw = nw, pair_ij = pair_ij, ρRlist = ρRlist,
		JuliaTask_nchunk_int = JuliaTask_nchunk_int, BLAS_nthread_epsilon = interaction_para.BLAS_nthread_epsilon,
		# multipolecorrection for 3D
		havempc = havempc, AR = AR, CR = CR, wingaverage = interaction_para.wingaverage,
		head_mpc_coff = head_mpc_coff, wing_mpc_coff = wing_mpc_coff,
		# ijji
		ijji_R = ijji_R, ijji_idx = ijji_idx, ijji_idx_chunk = ijji_idx_chunk,
		ijji_value_bare = ijji_value_bare, ijji_value_bare_iq = ijji_value_bare_iq,
		ijji_value_bare_head = ijji_value_bare_head, ijji_value_bare_head_mpc = ijji_value_bare_head_mpc,
		ijji_value = ijji_value, ijji_value_iq = ijji_value_iq, ijji_value_head = ijji_value_head,
		ijji_value_head_mpc = ijji_value_head_mpc, ijji_value_wing_mpc = ijji_value_wing_mpc,
		ϵ = ϵ,
		# ijij
		ijij_R = ijij_R, ijij_idx = ijij_idx, ijij_idx_chunk = ijij_idx_chunk,
		ijij_value_bare = ijij_value_bare, ijij_value_bare_iq = ijij_value_bare_iq,
		ijij_value_bare_head_mpc = ijij_value_bare_head_mpc,
		ijij_value = ijij_value, ijij_value_head_mpc = ijij_value_head_mpc, ijij_value_wing_mpc = ijij_value_wing_mpc,
		# iijj
		iijj_R = iijj_R, iijj_idx = iijj_idx, iijj_idx_chunk = iijj_idx_chunk,
		iijj_value_bare = iijj_value_bare, iijj_value_bare_iq = iijj_value_bare_iq,
		iijj_value_bare_head_mpc = iijj_value_bare_head_mpc,
		iijj_value = iijj_value, iijj_value_head_mpc = iijj_value_head_mpc, iijj_value_wing_mpc = iijj_value_wing_mpc,
		# pair
		pair = pair, pair_R = pair_R, pair_pR = pair_pR, pair_isp0 = pair_isp0, pair_uploidx = pair_uploidx,
		pair_idx = pair_idx, pair_idx_chunk = pair_idx_chunk,
		pair_value_bare = pair_value_bare, pair_value_bare_iq = pair_value_bare_iq,
		pair_value_bare_head_mpc = pair_value_bare_head_mpc,
		pair_value = pair_value, pair_value_iq = pair_value_iq,
		pair_value_head_mpc = pair_value_head_mpc, pair_value_wing_mpc = pair_value_wing_mpc,
		pair_atol = interaction_para.pair_atol, pair_hermitian_part = interaction_para.pair_hermitian_part,
	)
end
