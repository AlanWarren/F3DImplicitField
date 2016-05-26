class simple_vol_density(varying float density = 1;
			uniform float densityMultiplier = 1;
			uniform float occlusionSamples = 16;
			uniform float maxVariation = 0.3;)
{
	public void surface(output color Ci, Oi)
	{
        color Cdiff = 0;

        illuminance(P) {
            Cdiff += Cl;
        }

		Oi = density*VolumeField*densityMultiplier;
		Ci = Cdiff * Oi;//*(transmission(P, transform("world", "current", point(20, 50, 0)))[0]);
	}
}
