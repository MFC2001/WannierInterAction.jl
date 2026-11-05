
# It's extremely inefficient to calculate the interaction in real space.
# So this function is only used to check the validity of mirror correction.

using WannierInterAction

intfolder = "./interaction"
cell = read("./POSCAR", POSCAR)

chkdata = read("./wannier/wannier90.chk", wannier90_chk)
wave_ik = read("./wannier/pwscf.save/", chkdata, QE_wannier; period = Bool[1, 1, 1])
wannier = realspace(wave_ik; range = (3, 3, 3), cell)

write("./wannier/wannier_1.xsf", wannier[1], wannier90_xsf)
# wannier1 = read("./wannier/wannier90_00001.xsf", wannier90_xsf)

U = int(wannier[1], wannier[1], [0, 0, 0])
U = int(wannier[1], wannier[2], [1, 0, 0])
