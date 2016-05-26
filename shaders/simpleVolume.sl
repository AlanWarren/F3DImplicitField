class simpleVolume(	uniform float densityMultiplier = 1;
			varying float density = 1;)
{
	public void opacity(	output varying color Oi)
	{
		Oi = Os*VolumeField*densityMultiplier*density;
	}
	public void surface(	output varying color Ci, Oi)
	{
		uniform string raytype = "";
		rayinfo("type", raytype);

		this->opacity(Oi);
		if(raytype == "transmission") //early out
			return;

		varying color lightCol = 1;
		/*illuminance(P)
		{
			lightCol += Cl;
		}*/

		Ci = Oi*lightCol;
	}
}
