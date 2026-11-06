
using JLD2
using WannierInterAction

period = [1, 1, 1]
#1. Generate wannier pairs.
cell = read("./POSCAR_orb", POSCAR; period)
atom_orb_index = [
    Int[], [1], [2], [3], [4], [5], [6],
    [7, 8, 9],
]
pair = atom2orbpair(cell, 4, atom_orb_index)

#2. We calculate interaction between wannier basis under plane-wave basis.
intfolder = "./interaction"

# If we calculate MLWF by QE+wannierberri, only support QE at present.
chkdata = read("./wannier/w90data.chk.npz", wannierberri_chk)
wave_ik = read("./wannier/pwscf.save", chkdata, QE_wannier; period)
# If we calculate dielectric function by BerkeleyGW:
# please set `write_vcoul` in epsilon.inp
eps = read("./epsilon/epsmat.h5", epsmat_bgw)
# Calculate interaction. The output will be saved in outfolder with jld2 format.
# This program only need the kgrid size of wannier and dielectric function to be the same.
ρk = Rhok(eps, wave_ik)
int(eps, ρk; pair, outfolder=intfolder)

# If you don't need to calculate screened interaction:
# ρk = Rhok(wave_ik)
# int(ρk; pair, outfolder = intfolder)

# wave_ik can be used to plot the wannier functions in real space.
# wannierfolder = "./"
# jldsave(joinpath(wannierfolder, "wannier.jld2"); wannier = wave_ik)
# wanniers = realspace(wave_ik; cell, range = (3, 3, 3))
# for iw in eachindex(wanniers)
# 	write(joinpath(wannierfolder, "wannier_$iw.xsf"), wanniers[iw], wannier90_xsf)
# end


#3. Do mirror correction.
# Here we will obtain optimized parameters for the mirror correction of bare direct interaction, and then use the same geometric parameters for screened direct interaction.
kgrid = ρk.qgrid
lattice = chkdata.lattice
orblocat = chkdata.wannier_centres

# Do mirror correction for bare direct interaction.
U_bare = load(joinpath(intfolder, "U_bare.jld2"))
U_bare_mirror = MirrorCorrection(:gauss_optim, U_bare["interaction"], kgrid, lattice, orblocat;
    rcut=13, rcutmin=8, outfile=joinpath(intfolder, "U_bare_mirror.jld2"))

U_1 = load(joinpath(intfolder, "U_1.jld2"))
ϵ = real(U_1["ϵ"])
U_1_mirror = MirrorCorrection(:gauss, U_1["interaction"], kgrid, lattice, orblocat;
    αmat=U_bare_mirror.αmat, ϵ, outfile=joinpath(intfolder, "U_1_mirror.jld2"))


#4. You can prune or filter interaction between pairs.

# pair_atol = 1e-3
# pairint_bare = load(joinpath(intfolder, "pair_bare.jld2"))
# jldsave(joinpath(intfolder, "pair_prune", "pair_bare.jld2"); interaction = prune(pairint_bare, pair_atol))
# pairint_1 = load(joinpath(intfolder, "pair_1.jld2"))
# jldsave(joinpath(intfolder, "pair_prune", "pair_1.jld2"); interaction = prune(pairint_1, pair_atol))

# remain_pair = atom2orbpair(cell, 2.4, atom_orb_index)
# pairint_bare = load(joinpath(intfolder, "pair_bare.jld2"))
# jldsave(joinpath(intfolder, "pair_filter", "pair_bare.jld2"); interaction = filter_pair(pairint_bare, remain_pair))
# pairint_1 = load(joinpath(intfolder, "pair_1.jld2"))
# jldsave(joinpath(intfolder, "pair_filter", "pair_1.jld2"); interaction = filter_pair(pairint_1, remain_pair))
