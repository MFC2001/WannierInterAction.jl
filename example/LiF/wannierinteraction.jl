
using JLD2
using WannierInterAction

period = [1, 1, 1]
#1. Generate wannier pairs.
cell = read("./POSCAR", POSCAR; period)
# Li_sp3:1:4, F_p:5:7
atom_orb_index = [
	[1, 2, 3, 4],
	[5, 6, 7],
]
# n:2.013; nn:2.84681; nnn:3.4866; nnnn:4.026;
pair = atom2orbpair(cell, 4, atom_orb_index)

#2. We calculate interaction between wannier basis under plane-wave basis.
intfolder = "./interaction"

# If we calculate MLWF by QE+wannier90, only support this method at present.
chkdata = read("./wannier/wannier90.chk", wannier90_chk)
wave_ik = read("./wannier/pwscf.save", chkdata, QE_wannier; period)
# If we calculate dielectric function by BerkeleyGW:
eps = read("./epsilon/", epsmat_bgw; have_negq0 = true, negq0file = "./epsilon_negq0/eps0mat.h5")
# Calculate interaction. The output will be saved in intfolder with jld2 format.
# This program only need the kgrid size of wanner and dielectric function to be the same.
ρk = Rhok(eps, wave_ik)
int(eps, ρk; pair, outfolder = intfolder)

# If you don't need to calculate screened interaction:
# ρk = Rhok(wave_ik)
# int(ρk; pair, outfolder = intfolder)


# wave_ik can be used to plot the wannier functions in real space.
# wannierfolder = "./"
# jldsave(joinpath(wannierfolder, "wannier.jld2"); wannier = wave_ik)
# wanniers = wannier_realspace(wave_ik; range = (3, 3, 3))
# for iw in eachindex(wanniers)
# 	write(joinpath(wannierfolder, "wannier_$iw.xsf"), wannier90_xsf; cell, range = (3, 3, 3))
# end


#3. Do mirror correction.
# Here we will obtain optimized parameters for the mirror correction of bare direct interaction, and then use the same parameters for screened direct interaction.
kgrid = ρk.qgrid
lattice = chkdata.lattice
orblocat = chkdata.wannier_centres

# Do mirror correction for bare direct interaction.
U_bare = load(joinpath(intfolder, "U_bare.jld2"))
U_bare_mirror = MirrorCorrection(:gauss_optim, U_bare["interaction"], kgrid, lattice, orblocat;
	rcut = 13, rcutmin = 8, outfile = joinpath(intfolder, "U_bare_mirror.jld2"))

U_1 = load(joinpath(intfolder, "U_1.jld2"))
ϵ = real(U_1["ϵ"])
# The long-range behavior of screened U may not equal to U_bare/ϵ, so you can use the average of U_1/U_bare at large distance as the effective ϵ.
# But we don't know whether it's correct.
# result = longrange_ϵ(U_1["interaction"], U_bare["interaction"], lattice, orblocat; atol_range = 1e-2, atol_var = 1e-3)
# ϵ = result.ϵ

U_1_mirror = MirrorCorrection(:gauss, U_1["interaction"], kgrid, lattice, orblocat;
	αmat = U_bare_mirror.αmat, ϵ, outfile = joinpath(intfolder, "U_1_mirror.jld2"))


#4. You can truncate or filter interaction between pairs.

# pair_atol = 1e-3
# pairint_bare = load(joinpath(intfolder, "pair_bare.jld2"))
# jldsave(joinpath(intfolder, "pair_truncate", "pair_bare.jld2"); interaction = truncate(pairint_bare, pair_atol))
# pairint_1 = load(joinpath(intfolder, "pair_1.jld2"))
# jldsave(joinpath(intfolder, "pair_truncate", "pair_1.jld2"); interaction = truncate(pairint_1, pair_atol))

# remain_pair = generate_pair(2.5)
# pairint_bare = load(joinpath(intfolder, "pair_bare.jld2"))
# jldsave(joinpath(intfolder, "pair_filter", "pair_bare.jld2"); interaction = filter_pair(pairint_bare, remain_pair))
# pairint_1 = load(joinpath(intfolder, "pair_1.jld2"))
# jldsave(joinpath(intfolder, "pair_filter", "pair_1.jld2"); interaction = filter_pair(pairint_1, remain_pair))
