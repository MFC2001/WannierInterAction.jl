function _report_rhok(ρk::Rhok_wospin)
	println("calculation: wospin")
	println("periodicity: [", Int(ρk.period[1]), ", ", Int(ρk.period[2]), ", ", Int(ρk.period[3]), "]")
	println("nk = ", ρk.nk, "; nq = ", ρk.nq)
	println("fft_grid = ", ρk.fft_grid)
	println("nGlist_rho = ", length(ρk.Glist))
	println("nG_max_rho = ", maximum(ρk.nG_max))
	println("lattice:")
	for row in eachcol(parent(ρk.lattice))
		@printf("%13.6f %13.6f %13.6f\n", row...)
	end
	return nothing
end
function _report_rhok(ρk::Rhok_collinear)
	println("calculation: collinear")
	println("periodicity: [", Int(ρk.period[1]), ", ", Int(ρk.period[2]), ", ", Int(ρk.period[3]), "]")
	println("upindex: 1:", ρk.nw_up)
	println("dnindex: ", ρk.nw_up+1, ":", ρk.nw_up+ρk.nw_dn)
	println("nk = ", ρk.nk, "; nq = ", ρk.nq)
	println("fft_grid = ", ρk.fft_grid)
	println("nGlist_rho = ", length(ρk.Glist))
	println("nG_max_rho = ", maximum(ρk.nG_max))
	println("lattice:")
	for row in eachcol(parent(ρk.lattice))
		@printf("%13.6f %13.6f %13.6f\n", row...)
	end
	return nothing
end
function _report_rhok(ρk::Rhok_noncollinear)
	println("calculation: noncollinear")
	println("periodicity: [", Int(ρk.period[1]), ", ", Int(ρk.period[2]), ", ", Int(ρk.period[3]), "]")
	println("upindex: ", ρk.upindex)
	println("dnindex: ", ρk.dnindex)
	println("nk = ", ρk.nk, "; nq = ", ρk.nq)
	println("fft_grid = ", ρk.fft_grid)
	println("nGlist_rho = ", length(ρk.Glist))
	println("nG_max_rho = ", maximum(ρk.nG_max))
	println("lattice:")
	for row in eachcol(parent(ρk.lattice))
		@printf("%13.6f %13.6f %13.6f\n", row...)
	end
	return nothing
end
