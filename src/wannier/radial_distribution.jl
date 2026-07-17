function radial_distribution(ψ::AbstractArray{<:Number, 3}, minilattice,
	centre::Union{Nothing, AbstractVector{<:Real}} = nothing;
	type::Symbol = :norm2, Δr::Real = 0.1,
)
	# get data of r
	minilattice = Lattice(parent(minilattice))
	(nx, ny, nz) = size(ψ)
	r3D = Array{CartesianCoordinates{Float64}, 3}(undef, nx, ny, nz)
	for iz in Base.OneTo(nz), iy in Base.OneTo(ny), ix in Base.OneTo(nx)
		r3D[ix, iy, iz] = minilattice * ReducedCoordinates{Int}(ix, iy, iz)
	end
	ψ_norm2 = abs2.(ψ)
	if isnothing(centre)
		centre = sum(ψ_norm2 .* r3D) / sum(ψ_norm2)
	end
	minildot = parent(minilattice)
	minildot = SMatrix{3, 3}(transpose(minildot) * minildot)
	n = nx * ny * nz
	r = Vector{Float64}(undef, n)
	for i in Base.OneTo(n)
		Δr = r3D[i] - centre
		r[i] = sqrt(Δr ⋅ (minildot * Δr))
	end
	I = sortperm(r)
	r = r[I]
	ψ_norm2 = ψ_norm2[I]

	if type == :norm2
		return r, ψ_norm2
	elseif type == :norm
		return r, sqrt.(ψ_norm2)
	else
		error("Unsupported $type")
	end
end
