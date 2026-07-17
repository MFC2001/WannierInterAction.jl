struct RhoR_nobuffer{RhoT <: Number, RhokT <: AbstractRhok{RhoT}} <: AbstractRhoR{RhoT}
	ρk::RhokT
end
function RhoR(ρk::AbstractRhok{RhoT}) where {RhoT}
	return RhoR_nobuffer(ρk)
end
