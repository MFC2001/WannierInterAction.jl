function _int_make_hr_nwnwnR(int_value::AbstractArray{VT, 3}, R) where {VT}
	(nw, nw, nR) = size(int_value)

	N = nw * nw * nR
	pathI = Matrix{Int}(undef, 6, N)
	n = 0
	for I in CartesianIndices((nw, nw, nR))
		(i, j, iR) = Tuple(I)
		n += 1
		pathI[1:3, n] .= R[iR]
		pathI[4, n] = i
		pathI[5, n] = j
		pathI[6, n] = n
	end
	pathI = sortslices(pathI; dims = 2)
	perm = pathI[6, :]
	path = pathI[1:5, :]

	return HR(path, int_value[perm]; orbindex = collect(1:nw))
end
function _int_make_hrs_nwnwnRnω(int_value::AbstractArray{VT, 4}, R) where {VT}
	(nw, nw, nR, nω) = size(int_value)

	N = nw * nw * nR
	pathI = Matrix{Int}(undef, 6, N)
	values = Matrix{VT}(undef, N, nω)
	n = 0
	for iω in 1:nω
		for I in CartesianIndices((nw, nw, nR))
			(i, j, iR) = Tuple(I)
			n += 1
			pathI[1:3, n] .= R[iR]
			pathI[4, n] = i
			pathI[5, n] = j
			pathI[6, n] = n
			values[n, iω] = int_value[i, j, iR, iω]
		end
	end
	pathI = sortslices(pathI; dims = 2)
	perm = pathI[6, :]
	path = pathI[1:5, :]

	hrs = Vector{HR{VT}}(undef, nω)
	for iω in Base.OneTo(nω)
		value = values[perm, iω]
		hrs[iω] = HR(path, value; orbindex = collect(1:nw))
	end

	return hrs
end
# function _int_make_hr_nRnwnw(int_value::AbstractArray{VT, 3}, R) where {VT}
# 	(nR, nw, nw) = size(int_value)

# 	N = nR * nw * nw
# 	pathI = Matrix{Int}(undef, 6, N)
# 	n = 0
# 	for I in CartesianIndices((nR, nw, nw))
# 		(iR, i, j) = Tuple(I)
# 		n += 1
# 		pathI[1:3, n] .= R[iR]
# 		pathI[4, n] = i
# 		pathI[5, n] = j
# 		pathI[6, n] = n
# 	end
# 	pathI = sortslices(pathI; dims = 2)
# 	perm = pathI[6, :]
# 	path = pathI[1:5, :]

# 	return HR(path, int_value[perm]; orbindex = collect(1:nw))
# end
# function _int_make_hrs_nRnwnwnω(int_value::AbstractArray{VT, 4}, R) where {VT}
# 	(nR, nw, nw, nω) = size(int_value)

# 	N = nR * nw * nw
# 	pathI = Matrix{Int}(undef, 6, N)
# 	values = Matrix{VT}(undef, N, nω)
# 	n = 0
# 	for I in CartesianIndices((nR, nw, nw))
# 		(iR, i, j) = Tuple(I)
# 		n += 1
# 		pathI[1:3, n] .= R[iR]
# 		pathI[4, n] = i
# 		pathI[5, n] = j
# 		pathI[6, n] = n
# 		for iω in 1:nω
# 			values[n, iω] = int_value[iR, i, j, iω]
# 		end
# 	end
# 	pathI = sortslices(pathI; dims = 2)
# 	perm = pathI[6, :]
# 	path = pathI[1:5, :]

# 	hrs = Vector{HR{VT}}(undef, nω)
# 	for iω in Base.OneTo(nω)
# 		value = values[perm, iω]
# 		hrs[iω] = HR(path, value; orbindex = collect(1:nw))
# 	end

# 	return hrs
# end
