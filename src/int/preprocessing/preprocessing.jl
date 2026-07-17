
include("./_symmetrize_coulomb!.jl")
include("./buffer.jl")
include("./pair.jl")
include("./wispin.jl")
include("./wospin.jl")

function _preprocessing_int(ρk::AbstractRhok, interaction_para, multipolecorrection)
	return _preprocessing_int(ρk, 0, interaction_para, multipolecorrection)
end
