
using JLD2
using WannierInterAction

intfolder = "./interaction"

period = [1, 1, 0]
# For low-dimensional system, you need to specify `period`.
cell = read("./POSCAR", POSCAR; period)
# C1:pz; C2:pz.
atom_orb_index = [[1], [2]]
# n:1.42; nn:2.46; nnn:2.84; nnnn:4.26
pair = atom2orbpair(cell, 3, atom_orb_index)

chkdata = read("./wannier/wannier90.chk", wannier90_chk)
wave_ik = read("./wannier/pwscf.save/", chkdata, QE_wannier; period)

ρk = Rhok(wave_ik)
int(ρk; pair, outfolder = intfolder)

kgrid = ρk.qgrid
lattice = chkdata.lattice
orblocat = chkdata.wannier_centres

nk = 15
rcut = (nk - 1) * 2.13 / 2

U_bare = load(joinpath(intfolder, "U_bare.jld2"))
U_bare_mirror = MirrorCorrection(:keldysh_gauss_optim, U_bare["interaction"], kgrid, lattice, orblocat; period = cell.period,
	rcut, outfile = joinpath(intfolder, "U_bare_mirror.jld2"))
