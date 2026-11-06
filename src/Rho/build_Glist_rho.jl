"""
	build_Glist_rho(qgrid, fft_grid, rldot; normcutrho::Real = -1) -> (Glist_rho, Gidx_rho)
	build_Glist_rho(eps::DielectricPWBasis, fft_grid, rldot; normcutrho::Real = -1) -> (Glist_rho, Gidx_rho, Gidx_rho_negq0)

Get a proper Glist and Gidx for ρ at each q.
Glist_supp contains G which has at least one norm of q+G less than normcutrho and which isn't in eps.Glist.
Glist_rho = [eps.Glist; Glist_supp].
Gidx_rho is a Vector{Vector{Int}} which gives the Gidx for each q, note that Gidx_rho[iq][1:eps.nGidx[iq]] == eps.Gidx[iq]
"""
#这个函数中不能使用多线程循环，会导致得到存在随机性的Gidx
function build_Glist_rho(eps::DielectricPWBasis, fft_grid, rldot; normcutrho::Real = -1)

	qgrid = eps.qgrid

	# Glist_norm = _Glist_from_norm(rlattice, normcutrho)
	Glist_norm = _Glist_from_fftgrid(fft_grid)

	if normcutrho <= 0
		# 该设置保证所有的Glist_norm都不会被截断。
		normcutrho = maximum(G -> G ⋅ (rldot * G), Glist_norm) * 10
	end

	Glist_eps_set = Set(eps.Glist)
	normcutrho_sq = normcutrho^2
	nG_norm = length(Glist_norm)
	nq = length(qgrid)

	# judge G_supp
	Gbool_supp = falses(nG_norm)
	qG_norm_sq_supp = Vector{Vector{Float64}}(undef, nG_norm)
	for iG in Base.OneTo(nG_norm)
		G = Glist_norm[iG]
		G ∈ Glist_eps_set && continue
		qG_norm_square_iG = map(qgrid) do q
			qG = q + G
			qG ⋅ (rldot * qG)
		end
		if any(x -> x < normcutrho_sq, qG_norm_square_iG)
			Gbool_supp[iG] = true
			qG_norm_sq_supp[iG] = qG_norm_square_iG
		end
	end

	# collect G_supp
	Glist_supp = Glist_norm[Gbool_supp]
	qG_norm_sq_supp = qG_norm_sq_supp[Gbool_supp]
	# sort G_supp
	sort_keys = map(G -> G ⋅ (rldot * G), Glist_supp)
	perm = sortperm(sort_keys)
	permute!(Glist_supp, perm)
	permute!(qG_norm_sq_supp, perm)

	# Glist_rho
	Glist_rho = [eps.Glist; Glist_supp]

	# compute qG_norm_sq
	nG_rho = length(Glist_rho)
	qG_norm_sq = Matrix{Float64}(undef, nG_rho, nq)
	nG_eps = length(eps.Glist)
	for iq in Base.OneTo(nq)
		q = qgrid[iq]
		for iG in Base.OneTo(nG_eps)
			qG = q + Glist_rho[iG]
			qG_norm_sq[iG, iq] = qG ⋅ (rldot * qG)
		end
		for iG_supp in eachindex(qG_norm_sq_supp)
			iG_rho = nG_eps + iG_supp
			qG_norm_sq[iG_rho, iq] = qG_norm_sq_supp[iG_supp][iq]
		end
	end

	# Gidx_rho
	Gidx_rho = Vector{Vector{Int}}(undef, nq)
	for iq in Base.OneTo(nq)
		Gidx_eps_iq_set = Set(eps.Gidx[iq])
		Gidx_rho_iq_else = filter(Base.OneTo(nG_rho)) do iG
			iG ∈ Gidx_eps_iq_set ? false : (qG_norm_sq[iG, iq] < normcutrho_sq)
		end
		sort!(Gidx_rho_iq_else)
		Gidx_rho[iq] = [eps.Gidx[iq]; Gidx_rho_iq_else]
	end

	return Glist_rho, Gidx_rho
end
function build_Glist_rho(qgrid::Union{MonkhorstPack, RedKgrid, AbstractVector{<:AbstractVector{<:Real}}},
	fft_grid, rldot; normcutrho::Real = -1)

	if qgrid isa MonkhorstPack
		qgrid = RedKgrid(qgrid)
	end

	# Glist_norm = _Glist_from_norm(rlattice, normcutrho)
	Glist_norm = _Glist_from_fftgrid(fft_grid)

	nq = length(qgrid)

	if normcutrho <= 0
		Glist_rho = Glist_norm
		nG_rho = length(Glist_rho)
		Gidx_rho = Vector{Vector{Int}}(undef, nq)
		for iq in Base.OneTo(nq)
			Gidx_rho[iq] = collect(1:nG_rho)
		end
		return Glist_rho, Gidx_rho
	end

	normcutrho_sq = normcutrho^2
	nG_norm = length(Glist_norm)
	nq = length(qgrid)

	# judge G_rho
	qG_norm_sq_rho = Vector{Vector{Float64}}(undef, nG_norm)
	Gbool_rho = falses(nG_norm)
	for iG in Base.OneTo(nG_norm)
		G = Glist_norm[iG]
		qG_norm_sq_iG = map(qgrid) do q
			qG = q + G
			qG ⋅ (rldot * qG)
		end
		any(x -> x < normcutrho_sq, qG_norm_sq_iG) && begin
			Gbool_rho[iG] = true
			qG_norm_sq_rho[iG] = qG_norm_sq_iG
		end
	end

	# collect G_rho
	Glist_rho = Glist_norm[Gbool_rho]
	qG_norm_sq_rho = qG_norm_sq_rho[Gbool_rho]
	# sort G_rho
	sort_keys = map(G -> G ⋅ (rldot * G), Glist_rho)
	perm = sortperm(sort_keys)
	permute!(Glist_rho, perm)
	permute!(qG_norm_sq_rho, perm)

	# compute qG_norm_sq
	nG_rho = length(Glist_rho)
	qG_norm_sq = Matrix{Float64}(undef, nG_rho, nq)
	for iG_rho in Base.OneTo(nG_rho)
		qG_norm_sq[iG_rho, :] .= qG_norm_sq_rho[iG_rho]
	end

	# Gidx_rho
	Gidx_rho = Vector{Vector{Int}}(undef, nq)
	for iq in Base.OneTo(nq)
		Gidx_rho_iq = filter(Base.OneTo(nG_rho)) do iG
			qG_norm_sq[iG, iq] < normcutrho_sq
		end
		sort!(Gidx_rho_iq)
		Gidx_rho[iq] = Gidx_rho_iq
	end

	return Glist_rho, Gidx_rho
end
function _Glist_from_fftgrid(fft_grid)::Vector{ReducedCoordinates{Int}}
	total_G = prod(fft_grid)
	Glist = Vector{ReducedCoordinates{Int}}(undef, total_G)
	G1_range = fft_grid[1] / 2
	G1_range = (-floor(Int, G1_range)):(ceil(Int, G1_range)-1)
	G2_range = fft_grid[2] / 2
	G2_range = (-floor(Int, G2_range)):(ceil(Int, G2_range)-1)
	G3_range = fft_grid[3] / 2
	G3_range = (-floor(Int, G3_range)):(ceil(Int, G3_range)-1)
	idx = 0
	for G3 in G3_range, G2 in G2_range, G1 in G1_range
		idx += 1
		Glist[idx] = ReducedCoordinates{Int}(G1, G2, G3)
	end
	return Glist
end
function _Glist_from_norm(rlattice::ReciprocalLattice, maxnorm::Real)::Vector{ReducedCoordinates{Int}}
	a₁ = rlattice[:, 1]
	a₂ = rlattice[:, 2]
	a₃ = rlattice[:, 3]

	Ω = abs((a₁ × a₂) ⋅ a₃)

	h₁ = Ω / norm(a₂ × a₃)
	h₂ = Ω / norm(a₃ × a₁)
	h₃ = Ω / norm(a₁ × a₂)

	max_Grange = Int.(cld.(maxnorm + 1, [h₁, h₂, h₃])) .+ 1
	max_G1 = max_Grange[1]
	max_G2 = max_Grange[2]
	max_G3 = max_Grange[3]

	total_G = (2max_G1 + 1) * (2max_G2 + 1) * (2max_G3 + 1)
	Glist = Vector{ReducedCoordinates{Int}}(undef, total_G)
	idx = 0
	for G3 in (-max_G3):max_G3, G2 in (-max_G2):max_G2, G1 in (-max_G1):max_G1
		idx += 1
		Glist[idx] = ReducedCoordinates(G1, G2, G3)
	end
	return Glist
end
