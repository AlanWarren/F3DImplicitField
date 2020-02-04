# Rendering Volumes from Houdini in PRMan – Part 2
[prman – RSL pyro shader test](http://vimeo.com/47593515)

I'll provide a more thorough breakdown with pretty pics soon. For now, a quick rundown.

My ImplicitField DSO has access to the same scalar fields present in Houdini's smoke simulation. Of these fields, I'm using burn, heat, density, temperature and fuel as separate AOV's.

The most important scalar fields for the visual appearance of an explosion are temperature, heat, and density. In order to get that "fire" look, I map a ramp ( spline in RSL speak ) to the heat field, blending from black, orange, yellow, and white. I then scale the result by the temperature field. This gets me 90% of the above look. The rest is done by small modulations of the field data, and adding noise.

I'm using raytraced transmission shadows from my key light in the above video. In the future, I'm going to add support for Pixar's new areashadow() shadeop, and physically plausible shading.

I apologize for the low resolution. This isn't a final result though. The above video was very quick to render. So far I'm very impressed with prman's RiVolume() primitive compared to Blobbies. It's naturally a bit quicker, and it has access to the radiosity cache, which RiBlobby does not.
