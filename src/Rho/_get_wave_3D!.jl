@inline function _get_wave_3D!(wavek3D::AbstractArray{complexWaveT, 3}, wavek::AbstractVector{WaveT}, wavek_3DGmap) where {complexWaveT, WaveT}
	@inbounds begin
		fill!(wavek3D, zero(complexWaveT))
		for (idx, Gmap) in enumerate(wavek_3DGmap)
			wavek3D[Gmap] = wavek[idx]
		end
	end
	return wavek3D
end
