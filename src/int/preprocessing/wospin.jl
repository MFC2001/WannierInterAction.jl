function _preprocessing_int(ρk::Rhok_wospin{RhoT}, nω::Integer, interaction_para, multipolecorrection) where {RhoT}

    nw = ρk.nw
    qgrid = ρk.qgrid
    lattice = ρk.lattice
    Ω = ρk.Ω

    pair_ij = Set([(i, i) for i in 1:nw])
    ρRlist = Set([ReducedCoordinates{Int}(0, 0, 0)])

    if DIRECT ∈ interaction_para.interaction
        direct_R = gridindex(Val(:WignerSeitz), lattice, qgrid.kgrid_size)
        direct_nR = length(direct_R)
        direct_idx = vec([(iw1, iw2) for iw1 in 1:nw, iw2 in 1:nw])
        direct_value_bare = zeros(RhoT, nw, nw, direct_nR)
        direct_value_bare_iq = zeros(RhoT, nw, nw)
        direct_value_bare_head = zeros(RhoT, nw, nw)
        direct_value_bare_head_mpc = zeros(RhoT, nw, nw)
        direct_value = iszero(nω) ? nothing : zeros(RhoT, nw, nw, direct_nR, nω)
        direct_value_iq = iszero(nω) ? nothing : zeros(RhoT, nw, nw)
        direct_value_head = iszero(nω) ? nothing : zeros(RhoT, nw, nw, nω)
        direct_value_head_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, nω)
        direct_value_wing_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, nω)
        ϵ = iszero(nω) ? nothing : zeros(ComplexF64, nω)
    else
        direct_R = nothing
        direct_idx = nothing
        direct_value_bare = nothing
        direct_value_bare_iq = nothing
        direct_value_bare_head = nothing
        direct_value_bare_head_mpc = nothing
        direct_value = nothing
        direct_value_iq = nothing
        direct_value_head = nothing
        direct_value_head_mpc = nothing
        direct_value_wing_mpc = nothing
        ϵ = nothing
    end
    if EXCHANGE ∈ interaction_para.interaction
        exchange_R = gridindex(interaction_para.exchange_range)
        exchange_nR = length(exchange_R)
        union!(ρRlist, exchange_R)
        union!(ρRlist, map(R -> -R, exchange_R))
        pair_ij_add = vec([(i, j) for i in 1:nw, j in 1:nw])
        union!(pair_ij, pair_ij_add)
        exchange_R0idx = findfirst(R -> iszero(R), exchange_R)
        exchange_idx = [(iw1, iw2, iR) for iw1 in 1:nw, iw2 in 1:nw, iR in 1:exchange_nR]
        # without U-type interaction
        if !isnothing(exchange_R0idx)
            for iw in 1:nw
                exchange_idx[iw, iw, exchange_R0idx] = (0, 0, 0)
            end
            exchange_idx = vec(exchange_idx)
            filter!(idx -> !any(iszero, idx), exchange_idx)
        else
            exchange_idx = vec(exchange_idx)
        end
        exchange1_value_bare = zeros(RhoT, nw, nw, exchange_nR)
        exchange2_value_bare = zeros(RhoT, nw, nw, exchange_nR)
        exchange1_value_bare_iq = zeros(RhoT, nw, nw, exchange_nR)
        exchange2_value_bare_iq = zeros(RhoT, nw, nw, exchange_nR)
        exchange1_value_bare_head_mpc = zeros(RhoT, nw, nw, exchange_nR)
        exchange2_value_bare_head_mpc = zeros(RhoT, nw, nw, exchange_nR)
        exchange1_value = iszero(nω) ? nothing : zeros(RhoT, nw, nw, exchange_nR, nω)
        exchange2_value = iszero(nω) ? nothing : zeros(RhoT, nw, nw, exchange_nR, nω)
        exchange1_value_head_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, exchange_nR, nω)
        exchange2_value_head_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, exchange_nR, nω)
        exchange1_value_wing_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, exchange_nR, nω)
        exchange2_value_wing_mpc = iszero(nω) ? nothing : zeros(RhoT, nw, nw, exchange_nR, nω)
    else
        exchange_R = nothing
        exchange_idx = nothing
        exchange1_value_bare = nothing
        exchange2_value_bare = nothing
        exchange1_value_bare_iq = nothing
        exchange2_value_bare_iq = nothing
        exchange1_value_bare_head_mpc = nothing
        exchange2_value_bare_head_mpc = nothing
        exchange1_value = nothing
        exchange2_value = nothing
        exchange1_value_head_mpc = nothing
        exchange2_value_head_mpc = nothing
        exchange1_value_wing_mpc = nothing
        exchange2_value_wing_mpc = nothing
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

    pair_ij_to_idx = Dict{NTuple{2,Int},Int}(ij => idx for (idx, ij) in enumerate(pair_ij))
    ρRlist_to_idx = Dict{ReducedCoordinates{Int},Int}(R => i for (i, R) in enumerate(ρRlist))

    JuliaTask_nchunk_int = interaction_para.JuliaTask_nchunk_int
    # _idx分块。后续进行内存预分配需要预分块，如果不涉及大量内存分配，建议使用Julia原生多线程循环，会自动分块。
    # 目前没有任务依赖的内存分配，但是为保持拓展性，暂时保持手动任务分块。
    ρR0idx = ρRlist_to_idx[ReducedCoordinates{Int}(0, 0, 0)]
    if DIRECT ∈ interaction_para.interaction
        direct_idx = map(direct_idx) do (i, j)
            ii = pair_ij_to_idx[(i, i)]
            jj = pair_ij_to_idx[(j, j)]
            return (i, j, ρR0idx, ii, ρR0idx, jj)
        end
        direct_idx_chunk = _split_tasks(direct_idx, JuliaTask_nchunk_int)
    else
        direct_idx = nothing
        direct_idx_chunk = nothing
    end
    if EXCHANGE ∈ interaction_para.interaction
        exchange_iR2ρR = map(R -> ρRlist_to_idx[R], exchange_R)
        exchange_iR2ρR_neg = map(R -> ρRlist_to_idx[-R], exchange_R)
        exchange1_idx = map(exchange_idx) do (i, j, iR)
            ρR_idx = exchange_iR2ρR[iR]
            ji = pair_ij_to_idx[(j, i)]
            return (i, j, iR, ρR_idx, ji, ρR_idx, ji)
        end
        exchange2_idx = map(exchange_idx) do (i, j, iR)
            ρR_left = exchange_iR2ρR[iR]
            ρR_right = exchange_iR2ρR_neg[iR]
            ji = pair_ij_to_idx[(j, i)]
            ij = pair_ij_to_idx[(i, j)]
            return (i, j, iR, ρR_left, ji, ρR_right, ij)
        end
        exchange1_idx_chunk = _split_tasks(exchange1_idx, JuliaTask_nchunk_int)
        exchange2_idx_chunk = _split_tasks(exchange2_idx, JuliaTask_nchunk_int)
    else
        exchange1_idx = nothing
        exchange2_idx = nothing
        exchange1_idx_chunk = nothing
        exchange2_idx_chunk = nothing
    end
    if PAIR ∈ interaction_para.interaction
        pair_idx_ρR = map(pair) do (i, j, ipR)
            ji = pair_ij_to_idx[(j, i)] # note!
            iρR = ρRlist_to_idx[pair_pR[ipR]]
            return (iρR, ji)
        end
        pair_idx = map(pair_idx) do (ip1, ip2)
            (ρR_left, ji_left) = pair_idx_ρR[ip1]
            (ρR_right, ji_right) = pair_idx_ρR[ip2]
            return (ip1, ip2, ρR_left, ji_left, ρR_right, ji_right)
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
        @info "need to check multipolecorrection!!!"
        println("need to check multipolecorrection!!!")
        havempc = true
        AR, CR = multipolecorrection(lattice, ρRlist, pair_ij)
        head_mpc_coff = CoulombScale * 4π / Ω
        wing_mpc_coff = sqrt(head_mpc_coff)
    end

    return (interaction=interaction_para.interaction, nw=nw, pair_ij=pair_ij, ρRlist=ρRlist,
        JuliaTask_nchunk_int=JuliaTask_nchunk_int, BLAS_nthread_epsilon=interaction_para.BLAS_nthread_epsilon,
        # multipolecorrection for 3D
        havempc=havempc, AR=AR, CR=CR, wingaverage=interaction_para.wingaverage,
        head_mpc_coff=head_mpc_coff, wing_mpc_coff=wing_mpc_coff,
        # direct
        direct_R=direct_R, direct_idx=direct_idx, direct_idx_chunk=direct_idx_chunk,
        direct_value_bare=direct_value_bare, direct_value_bare_iq=direct_value_bare_iq,
        direct_value_bare_head=direct_value_bare_head, direct_value_bare_head_mpc=direct_value_bare_head_mpc,
        direct_value=direct_value, direct_value_iq=direct_value_iq, direct_value_head=direct_value_head,
        direct_value_head_mpc=direct_value_head_mpc, direct_value_wing_mpc=direct_value_wing_mpc,
        ϵ=ϵ,
        # exchange
        exchange_R=exchange_R, exchange_idx=exchange_idx,
        exchange1_idx=exchange1_idx, exchange2_idx=exchange2_idx,
        exchange1_idx_chunk=exchange1_idx_chunk, exchange2_idx_chunk=exchange2_idx_chunk,
        exchange1_value_bare=exchange1_value_bare, exchange2_value_bare=exchange2_value_bare,
        exchange1_value_bare_iq=exchange1_value_bare_iq, exchange2_value_bare_iq=exchange2_value_bare_iq,
        exchange1_value_bare_head_mpc=exchange1_value_bare_head_mpc, exchange2_value_bare_head_mpc=exchange2_value_bare_head_mpc,
        exchange1_value=exchange1_value, exchange2_value=exchange2_value,
        exchange1_value_head_mpc=exchange1_value_head_mpc, exchange2_value_head_mpc=exchange2_value_head_mpc,
        exchange1_value_wing_mpc=exchange1_value_wing_mpc, exchange2_value_wing_mpc=exchange2_value_wing_mpc,
        # pair
        pair=pair, pair_R=pair_R, pair_pR=pair_pR, pair_isp0=pair_isp0, pair_uploidx=pair_uploidx,
        pair_idx=pair_idx, pair_idx_chunk=pair_idx_chunk,
        pair_value_bare=pair_value_bare, pair_value_bare_iq=pair_value_bare_iq,
        pair_value_bare_head_mpc=pair_value_bare_head_mpc,
        pair_value=pair_value, pair_value_iq=pair_value_iq,
        pair_value_head_mpc=pair_value_head_mpc, pair_value_wing_mpc=pair_value_wing_mpc,
        pair_atol=interaction_para.pair_atol, pair_hermitian_part=interaction_para.pair_hermitian_part,
    )
end
