# WannierInterAction

[![Build Status](https://github.com/MFC2001/WannierInterAction.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MFC2001/WannierInterAction.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/github/MFC2001/WannierInterAction.jl/graph/badge.svg?token=BW8CIEPMQM)](https://codecov.io/github/MFC2001/WannierInterAction.jl)

A Julia package for evaluating the Coulomb interaction between Wannier pairs
(ordered pairs of Wannier functions) in a plane-wave basis,
and for applying mirror corrections to the long-range part of that interaction,
which would otherwise contain the contributions of all periodic images:

$$
F_{p_1\boldsymbol{R}_1,p_2\boldsymbol{R}_2}
= \iint \mathrm{d}\boldsymbol{r}_1 \mathrm{d}\boldsymbol{r}_2
\psi^{*}_{i_{p_1}\boldsymbol{R}_1}(\boldsymbol{r}_1) \psi^{*}_{j_{p_2}\boldsymbol{R}_2+\boldsymbol{R}_{p_2}}(\boldsymbol{r}_2)
F(\boldsymbol{r}_1,\boldsymbol{r}_2) \psi_{j_{p_1}\boldsymbol{R}_1+\boldsymbol{R}_{p_1}}(\boldsymbol{r}_1) \psi_{i_{p_2}\boldsymbol{R}_2}(\boldsymbol{r}_2)
$$

where $F$ is either the bare or the screened Coulomb interaction.
Currently only **3D** and **2D** systems are supported.

The resulting pair interactions can be passed to
[`WannierExcitonModel.jl`](https://github.com/MFC2001/WannierExcitonModel.jl)
to construct the Bethe–Salpeter equation (BSE) Hamiltonian and obtain excitonic band structures.

## Installation

The package requires Julia 1.10 or later and can be installed with the Julia package manager.
From [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), press `]` to enter
the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) and run:

```julia
pkg> add WannierInterAction
```

or, equivalently, via [`Pkg.jl`](https://pkgdocs.julialang.org/v1/):

```julia
julia> import Pkg; Pkg.add("WannierInterAction")
```

## Quick start

Entering `?WannierInterAction` in the Julia REPL prints the built-in documentation;
the same help is available for every exported function (for example `?Rhok`, `?int`
and `?MirrorCorrection`).

The workflow below follows the template scripts in [`examples`](examples). It assumes that
the plane-wave calculation and the construction of the Wannier functions are already done,
and that `POSCAR`, `pwscf.save`, the Wannier checkpoint and — for screened interactions —
a dielectric matrix are available.

### 1. Define the cell and the Wannier pairs

```julia
using WannierInterAction

period = [1, 1, 1]                              # [1, 1, 0] indicates a 2D slab system which is unperiodic in z-direction
cell = read("./POSCAR", POSCAR; period)
atom_orb_index = [[1, 2], [3], [4]]              # Wannier orbitals on each atom (`Int[]` if none)
pair = atom2orbpair(cell, 4.0, atom_orb_index)   # all pairs within 4 Å
```

Each element of `pair` is a tuple `(i, j, R)`: two orbital indices and the lattice vector
that separates them.

### 2. Load the wavefunctions and the dielectric function

```julia
chkdata = read("./wannier/wannier90.chk", wannier90_chk)   # wannierberri_chk for `*.chk.npz`
wave_ik = read("./wannier/pwscf.save", chkdata, QE_wannier; period)

# Only needed for screened interactions; `epsmat_bgw` also reads `vcoul`,
# so set `write_vcoul` in `epsilon.inp`. `epsmat_yambo` is supported as well.
eps = read("./epsilon/epsmat.h5", epsmat_bgw)
```

### 3. Compute the interactions

```julia
ρk = Rhok(eps, wave_ik)                        # Rhok(wave_ik) without screening
int(eps, ρk; pair, outfolder="./interaction")  # int(ρk; ...) without screening
```

`int` evaluates the direct, exchange and pair interactions and writes them to `outfolder`
as JLD2 files: `U_bare.jld2` for the bare Coulomb interaction, `U_<i>.jld2` for each
frequency of the dielectric function (if one is given), `J1_*`/`J2_*` for exchange and
`pair_*` for the pair interaction, together with a `WannierInterAction.log` file. The
k-grid of the wavefunctions and that of the dielectric function must have the same size.
The calculation is threaded, so start Julia with `julia -t auto` (or set
`JULIA_NUM_THREADS`).

### 4. Apply the mirror correction

```julia
using JLD2

kgrid, lattice, orblocat = ρk.qgrid, chkdata.lattice, chkdata.wannier_centres

U_bare = load("./interaction/U_bare.jld2")
U_bare_mirror = MirrorCorrection(:gauss_optim, U_bare["interaction"],
    kgrid, lattice, orblocat; rcut=13, outfile="./interaction/U_bare_mirror.jld2")
```

Use `:gauss_optim`/`:gauss` for 3D systems and `:keldysh_gauss_optim`/`:keldysh_gauss`
for 2D ones, passing `period=cell.period` to the latter. The `_optim` variants fit the
model parameters to the interaction beyond `rcut`, and the fitted parameters can be reused
for the screened interactions (3D case):

```julia
U_1 = load("./interaction/U_1.jld2")
U_1_mirror = MirrorCorrection(:gauss, U_1["interaction"], kgrid, lattice, orblocat;
    αmat=U_bare_mirror.αmat, ϵ=real(U_1["ϵ"]), outfile="./interaction/U_1_mirror.jld2")
```

`MirrorCorrection` returns a `NamedTuple` with `Uwithlr` (the corrected interaction, ready
to be used in later calculations), `Umirror` (the real-space correction) and `Ushort` (the
interaction with the periodic tail removed only); the `_optim` variants additionally return
the fitted parameters. Afterwards, `prune` and `filter_pair` can be used to drop negligible
entries or to restrict the interaction to a smaller set of pairs.

Note: If the k-grid used in the calculation is too small, the direct interaction may not converge to its theoretical behavior.
This can cause the mirror correction to take an excessively long time and the resulting `Ushort` not to decay to zero.
Users can find more details in `MirrorCorrection.log`, which is stored by default in the same directory as the specified output file.

In addition, the only effective Coulomb potential currently used for 2D systems is `keldysh_gauss`,
which is derived from a zero-thickness continuous dielectric model. Therefore, the 2D mirror
correction with keldysh_gauss is only suitable for sufficiently thin 2D materials. If the material
is not thin enough, a larger k-grid may be required to compensate for this effect.

### Examples

| Template | System | Highlights |
| --- | --- | --- |
| [`examples/LiF`](examples/LiF) | 3D bulk | full workflow with BerkeleyGW screening and `:gauss_optim`; `wannierberri` preprocessing (`pre.py`); reuse of fitted parameters; commented `prune`/`filter_pair` examples |
| [`examples/MoS2`](examples/MoS2) | 2D slab, spinor Wannier functions | noncollinear spin treatment (`upindex`/`dnindex`), `:keldysh_gauss_optim`, `r₀_optim=true`, `pair_uplo=:U` |
| [`examples/Graphene`](examples/Graphene) | 2D, bare interaction | interaction without screening (`Rhok(wave_ik)`, `int(ρk; ...)`) |
| [`examples/Graphone`](examples/Graphone) | 2D, spin-polarized | separate spin channels, Yambo dielectric function, explicit `PlaneWaveCoulomb(:slab, ρk)`, real-space Wannier plots via `realspace`/`wannier90_xsf` |
