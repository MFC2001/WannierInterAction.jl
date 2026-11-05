# WannierInterAction

[![Build Status](https://github.com/MFC2001/WannierInterAction.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MFC2001/WannierInterAction.jl/actions/workflows/CI.yml?query=branch%3Amaster)

A Julia package for evaluating Coulomb interactions between Wannier pairs
(ordered pairs of Wannier functions) and performing mirror corrections for long-range terms.
The Coulomb interaction between pairs can be used 
with [`WannierExcitonModel.jl`](https://github.com/MFC2001/WannierExcitonModel.jl) 
to construct BSE (Bethe–Salpeter equation) Hamiltonian and calculate excitonic band structure.

For a quick start, enter `?WannierInterAction` in the Julia REPL.
If you are using this package for the first time, we recommend referring to the template scripts in examples, 
which demonstrate how to call the functions in typical cases.
