"""
	int(w₁::wannier90_xsf{<: Number}, w₂::wannier90_xsf{<: Number}, R = [0, 0, 0])
	int(w₁::AbstractArray{<:Number, 3}, w₂::AbstractArray{<:Number, 3}, range, lattice, R)

∫|ψ₁(r₁)|²|ψ₂(r₂)|²1/|r₁ - r₂ - R| dr₁dr₂

!!! unusable!
"""
function int(w₁::wannier90_xsf{<: Number}, w₂::wannier90_xsf{<: Number}, R = [0, 0, 0])
	@assert w₁.grid_size == w₂.grid_size "Two wannier functions must be defined on the same grid"

	minilattice = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		minilattice[:, i] = w₁.grid_vecs[:, i] / (w₁.grid_size[i] - 1)
	end

	grid_size_uc = Vector{Int}(undef, 3)
	for i in 1:3
		l = norm(w₁.cell.lattice[:, i])
		minil = norm(minilattice[:, i])
		grid_size_uc[i] = round(Int, l / minil)
	end
	miniR = R .* grid_size_uc

	return _interaction_realspace(w₁.value, w₂.value, minilattice, miniR)
end
function int(w₁::AbstractArray{<:Number, 3}, w₂::AbstractArray{<:Number, 3}, range, lattice, R)
	@assert size(w₁) == size(w₂) "Two wannier functions must be defined on the same grid"

	grid_size = size(w₁)
	grid_size_uc = Int.(grid_size ./ range)

	minilattice = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		minilattice[:, i] = lattice[:, i] / grid_size_uc[i]
	end
	miniR = R .* grid_size_uc

	return _interaction_realspace(w₁, w₂, minilattice, miniR)
end
"""
ψ₁ and ψ₂ are defined on the same grid, and the grid should be the same.
minilattice is the basis of the grid.
R is the translation vector between ψ₁ and ψ₂, in unit of minilattice.
"""
function _interaction_realspace(
	ψ₁::AbstractArray{<:Number, 3}, ψ₂::AbstractArray{<:Number, 3},
	minilattice, miniR::AbstractVector{<:Integer},
)
	# 网格点自作用能的计算方法为固定一个格点，限制另一个格点在miniWS胞内。
	# 如果两个格点都限制在miniWS胞内，则难以进行数值积分。
	# for small |r₁ - r₂ - R|
	miniWS = WignerSeitz(minilattice)
	rsz2 = minimum(miniWS.offsets) / 2
	min_miniWS = Vector{Float64}(undef, 3)
	max_miniWS = Vector{Float64}(undef, 3)
	for i in 1:3
		(min_miniWS[i], max_miniWS[i]) = extrema(k -> k[i], miniWS.normals)
	end
	min_miniWS = Vec3{Float64}(min_miniWS)
	max_miniWS = Vec3{Float64}(max_miniWS)
	Ω = abs(det(parent(minilattice)))
	# 同一格点
	function integrand_0(r)
		r² = r[1]^2 + r[2]^2 + r[3]^2
		r² < rsz2 && return 0.0
		r ∈ miniWS && return 1 / sqrt(r²)
		return 0.0
	end
	(result, err) = hcubature(integrand_0, min_miniWS, max_miniWS; norm = abs, rtol = 1e-8, atol = 1e-10)
	selfenergy_0 = (result + 2π * rsz2) / Ω
	# 不同格点
	function integrand_Δr(Δr)
		Δr = Vec3{Float64}(Δr)
		function integrand(r)
			if r ∈ miniWS
				r² = (r[1] + Δr[1])^2 + (r[2] + Δr[2])^2 + (r[3] + Δr[3])^2
				return 1 / sqrt(r²)
			else
				return 0.0
			end
		end
		(result, err) = hcubature(integrand, min_miniWS, max_miniWS; norm = abs, rtol = 1e-8, atol = 1e-10)
		return result / Ω
	end

	minildot = parent(minilattice)
	minildot = SMatrix{3, 3}(transpose(minildot) * minildot)
	@assert size(ψ₁) == size(ψ₂) "ψ₁ and ψ₂ must be defined on the same grid"
	(nx, ny, nz) = size(ψ₁)

	ρ₁ = abs2.(ψ₁)
	ρ₂ = abs2.(ψ₂)

	#sum only with selfenergy_0
	value = 0.0
	l = ReentrantLock()
	Threads.@threads for I in CartesianIndices((nx, ny, nz))
		(x₁, y₁, z₁) = Tuple(I)
		xR = x₁ - miniR[1]
		yR = y₁ - miniR[2]
		zR = z₁ - miniR[3]
		is_same_grid = true
		value_I = 0.0
		for z₂ in Base.OneTo(nz)
			is_same_grid = z₂ == zR
			for y₂ in Base.OneTo(ny)
				is_same_grid = is_same_grid && y₂ == yR
				for x₂ in Base.OneTo(nx)
					if is_same_grid && x₂ == xR
						value_I += ρ₁[x₁, y₁, z₁] * ρ₂[x₂, y₂, z₂] * selfenergy_0
					else
						r = ReducedCoordinates(xR - x₂, yR - y₂, zR - z₂)
						r = sqrt(r ⋅ (minildot * r))
						value_I += ρ₁[x₁, y₁, z₁] * ρ₂[x₂, y₂, z₂] / r
					end
				end
			end
		end
		lock(l) do
			value += value_I
		end
	end

	T₁ = sum(ρ₁)
	T₂ = sum(ρ₂)
	value *= CoulombScale / (T₁ * T₂)

	return value
end
