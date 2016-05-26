	imager
	rps_bg( color background = 1; )
	{
        color sky1 = color(0.67, 0.78, 0.88);
        color sky2 = color(0.57, 0.68, 0.78);
        color sky3 = color(0.23, 0.53, 0.79);
        color sky4 = color(0.23, 0.53, 0.89);

        color sp = spline("linear", v, sky1, sky2, sky3, sky4);
		Ci += sp;
                Oi = 1;
	}
