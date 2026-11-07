
using JLD2
using WannierInterAction

# For low-dimensional system, you need to specify `period`.
period = [1, 1, 0]
cell = read("./POSCAR", POSCAR; period)
orbital_cell = Cell(cell.lattice, [[2/3, 1/3, 1/2]], "R"; name=["C"], period)
# 1:up, 2:dn.
atom_orb_index = [[1, 2]]
# shell:2.54, 4.4, 5.07, 6.7, 7.6, 8.8, 9.1, 10.1
pair = atom2orbpair(orbital_cell, 10.5, atom_orb_index)

intfolder = "./interaction"

chkdata_up = read("../wannier/up/wannier90.chk", wannier90_chk)
wave_ik_up = read("../wannier/up/pwscf.save", chkdata_up, QE_wannier; wfcspin=:up, period)
chkdata_dn = read("../wannier/dn/wannier90.chk", wannier90_chk)
wave_ik_dn = read("../wannier/dn/pwscf.save", chkdata_dn, QE_wannier; wfcspin=:dn, period)

# you can obtain the wannier function in real space,
wannier = realspace(wave_ik_up; cell)
write("./wannierup.xsf", wannier[1], wannier90_xsf)
wannier = realspace(wave_ik_dn; cell)
write("./wannierdn.xsf", wannier[1], wannier90_xsf)

# BerkeleyGW don't support magnetic system, so we use yambo to calculate the dielectric function.
eps = read("./from_node31/epsilon", epsmat_yambo)
# We need to provide wave functions for both spin up and spin down.
ρk = Rhok(eps, wave_ik_up, wave_ik_dn)

# Use slab truncated Coulomb interaction.
coulomb = PlaneWaveCoulomb(:slab, ρk)
int(eps, ρk; coulomb, pair, outfolder=intfolder)


kgrid = ρk.qgrid
lattice = chkdata_up.lattice
orblocat_up = chkdata_up.wannier_centres
orblocat_dn = chkdata_dn.wannier_centres
orblocat = [orblocat_up; orblocat_dn]

nk = 15
rcut = 10 * 3 / 2

# Do mirror correction for bare direct interaction.
U_bare = load(joinpath(intfolder, "U_bare.jld2"))
MirrorCorrection(:keldysh_gauss_optim, U_bare["interaction"], kgrid, lattice, orblocat; period=cell.period,
    rcut, outfile=joinpath(intfolder, "U_bare_mirror.jld2"))

U_1 = load(joinpath(intfolder, "U_1.jld2"))
MirrorCorrection(:keldysh_gauss_optim, U_1["interaction"], kgrid, lattice, orblocat; period=cell.period,
    rcut, r₀_optim=true, outfile=joinpath(intfolder, "U_1_mirror.jld2"))
