"""
	WannierInterAction

A Julia package used to evaluate the interaction between Wannier pairs and do the mirror correction.
See also [`Rhok`](@ref), [`int`](@ref) and [`MirrorCorrection`](@ref).
"""
module WannierInterAction

using WannierInterActionBase

using WannierInterActionBase.StaticArrays
using WannierInterActionBase.StructEquality

using WannierInterActionBase.JLD2

Base_exports = names(WannierInterActionBase; all = false)
for sym in Base_exports
	@eval export $sym
end

using LinearAlgebra
using Printf
using Dates: Dates
using Logging
using LoggingExtras

# optimization
using Optim

# fft
using AbstractFFTs
using FFTW

# interpolation
using Interpolations

include("./wannier/wannier.jl")
include("./Rho/Rho.jl")
include("./Coulomb/Coulomb.jl")
include("./int/int.jl")
include("./MirrorCorrection/MirrorCorrection.jl")

end
