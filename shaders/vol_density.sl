#include "pyro/blackbody.h"
// Float Ramp
//---------------------------------------
#define RAMP_PARMSDFLT(PRE,name,POST) \
   string  PRE##name##w##POST##b[] = {"linear", "linear"} ; \
   float   PRE##name##w##POST##k[] = {0,1};   \
   float   PRE##name##w##POST##v[] = {0,1}

#define RAMP_PARMS(PRE,name,POST)    \
   string  PRE##name##w##POST##b[] ; \
   float   PRE##name##w##POST##k[] ; \
   float   PRE##name##w##POST##v[]

#define RAMP_ARGS(PRE,name,POST)   \
           PRE##name##w##POST##b , \
           PRE##name##w##POST##k , \
           PRE##name##w##POST##v

// Color Ramp
//---------------------------------------
#define RAMPC_PARMSDFLT(PRE,name,POST) \
   string  PRE##name##w##POST##b[] = {"linear", "linear"} ; \
   float   PRE##name##w##POST##k[] = {0,1};   \
   vector  PRE##name##w##POST##v[] = {0,1}

#define RAMPC_PARMS(PRE,name,POST)    \
   string  PRE##name##w##POST##b[] ; \
   float   PRE##name##w##POST##k[] ; \
   color   PRE##name##w##POST##v[]

#define RAMPC_ARGS(PRE,name,POST)   \
           PRE##name##w##POST##b , \
           PRE##name##w##POST##k , \
           PRE##name##w##POST##v



// Field shaping parameters
//---------------------------------------
#define FSHAPE_PARMSDFLT(PRE,POST)       \
   float   PRE##tlabel##POST    = 0    ; \
   string  PRE##label##POST     = ""; \
   string  PRE##comment##POST   = ""; \
   float   PRE##enable##POST    = 1    ; \
   string  PRE##field##POST     = "density"; \
   float   PRE##tsharp##POST    = 0    ; \
   float   PRE##sharpe##POST    = 0.01 ; \
   float   PRE##sharpk##POST    = 4    ; \
   float   PRE##tsoft##POST     = 0    ; \
   float   PRE##softe##POST     = 0.8  ; \
   float   PRE##softk##POST     = 1    ; \
   float   PRE##trng##POST      = 0    ; \
   float   PRE##rngsl##POST     = 0    ; \
   float   PRE##rngsh##POST     = 1    ; \
   float   PRE##rngtl##POST     = 0    ; \
   float   PRE##rngth##POST     = 1    ; \
   RAMP_PARMSDFLT(PRE,rng,POST)

#define FSHAPE_PARMS(PRE,POST)    \
   float   PRE##tlabel##POST    ; \
   string  PRE##label##POST     ; \
   string  PRE##comment##POST   ; \
   float   PRE##enable##POST    ; \
   string  PRE##field##POST     ; \
   float   PRE##tsharp##POST    ; \
   float   PRE##sharpe##POST    ; \
   float   PRE##sharpk##POST    ; \
   float   PRE##tsoft##POST     ; \
   float   PRE##softe##POST     ; \
   float   PRE##softk##POST     ; \
   float   PRE##trng##POST      ; \
   float   PRE##rngsl##POST     ; \
   float   PRE##rngsh##POST     ; \
   float   PRE##rngtl##POST     ; \
   float   PRE##rngth##POST     ; \
   RAMP_PARMS(PRE,rng,POST)

#define FSHAPE_ARGS(PRE,POST)  \
           PRE##tlabel##POST , \
           PRE##label##POST  , \
           PRE##comment##POST, \
           PRE##enable##POST , \
           PRE##field##POST  , \
           PRE##tsharp##POST , \
           PRE##sharpe##POST , \
           PRE##sharpk##POST , \
           PRE##tsoft##POST  , \
           PRE##softe##POST  , \
           PRE##softk##POST  , \
           PRE##trng##POST   , \
           PRE##rngsl##POST  , \
           PRE##rngsh##POST  , \
           PRE##rngtl##POST  , \
           PRE##rngth##POST  , \
           RAMP_ARGS(PRE,rng,POST)

// Noise parameters
// XXX: no vector4
//---------------------------------------
#define NOISE_PARMSDFLT(PRE,POST)          \
   float   PRE##enable##POST    = 1;       \
   string  PRE##field##POST     = "none";  \
   float   PRE##rngsl##POST     = 0;       \
   float   PRE##rngsh##POST     = 0.3;     \
   RAMP_PARMSDFLT(PRE,rng,POST);           \
   float   PRE##amp##POST       = 1;       \
   string  PRE##ampm##POST      = "none";  \
   float   PRE##is4d##POST      = 0;       \
   vector  PRE##freq##POST      = 1;       \
   vector  PRE##off##POST       = 0;       \
   float   PRE##fw##POST        = 1;       \
   float   PRE##oct##POST       = 8;       \
   string  PRE##octm##POST      = "none";  \
   float   PRE##lac##POST       = 2.01234; \
   string  PRE##lacm##POST      = "none";  \
   float   PRE##gain##POST      = 0.5;     \
   string  PRE##gainm##POST     = "none";  \
   float   PRE##expon##POST     = 1;       \
   string  PRE##exponm##POST    = "none";  \
   float   PRE##rngol##POST     = 0;       \
   float   PRE##rngoh##POST     = 1

#define NOISE_PARMS(PRE,POST)  \
   float   PRE##enable##POST ; \
   string  PRE##field##POST  ; \
   float   PRE##rngsl##POST  ; \
   float   PRE##rngsh##POST  ; \
   RAMP_PARMS(PRE,rng,POST)  ; \
   float   PRE##amp##POST    ; \
   string  PRE##ampm##POST   ; \
   float   PRE##is4d##POST   ; \
   vector  PRE##freq##POST   ; \
   vector  PRE##off##POST    ; \
   float   PRE##fw##POST     ; \
   float   PRE##oct##POST    ; \
   string  PRE##octm##POST   ; \
   float   PRE##lac##POST    ; \
   string  PRE##lacm##POST   ; \
   float   PRE##gain##POST   ; \
   string  PRE##gainm##POST  ; \
   float   PRE##expon##POST  ; \
   string  PRE##exponm##POST ; \
   float   PRE##rngol##POST  ; \
   float   PRE##rngoh##POST

#define NOISE_ARGS(PRE,POST) \
   PRE##enable##POST  ,      \
   PRE##field##POST   ,      \
   PRE##rngsl##POST   ,      \
   PRE##rngsh##POST   ,      \
   RAMP_ARGS(PRE,rng,POST),  \
   PRE##amp##POST     ,      \
   PRE##ampm##POST    ,      \
   PRE##is4d##POST    ,      \
   PRE##freq##POST    ,      \
   PRE##off##POST     ,      \
   PRE##fw##POST      ,      \
   PRE##oct##POST     ,      \
   PRE##octm##POST    ,      \
   PRE##lac##POST     ,      \
   PRE##lacm##POST    ,      \
   PRE##gain##POST    ,      \
   PRE##gainm##POST   ,      \
   PRE##expon##POST   ,      \
   PRE##exponm##POST  ,      \
   PRE##rngol##POST   ,      \
   PRE##rngoh##POST

#define NOISE_MODS(PRE,POST) \
   PRE##ampm##POST    ,      \
   PRE##octm##POST    ,      \
   PRE##lacm##POST    ,      \
   PRE##gainm##POST   ,      \
   PRE##exponm##POST

#define C_PI               3.14159265358979323846  // pi 
#define C_PI2              6.28318530717958647692       // pi*2
#define VONE               (1,1,1)
#define CONE               color(1,1,1)
#define VZERO              (0,0,0)
#define C_RLUM             (0.2125)
#define C_GLUM             (0.7154)
#define C_BLUM             (0.0721)
#define C_1_SQRT2          0.70710678118654752440  // 1/sqrt(2)
#define C_SQRT2            1.41421356237309504880  // sqrt(2) 
#define C_1_SQRT2          0.70710678118654752440  // 1/sqrt(2) 
#define C_SQRT3            1.73205080756887729353       // sqrt(3)
#define C_1_SQRT3          0.57735026918962576451       // 1/sqrt(3)

color modc(color a, b)
{
    return color(mod(a[0], b[0]),
                 mod(a[1], b[1]),
                 mod(a[2], b[2]));
}

color smoothstepc(color a, b, bias)
{
    return color(smoothstep(a[0], b[0], bias[0]),
                 smoothstep(a[1], b[1], bias[1]),
                 smoothstep(a[2], b[2], bias[2]));
}

color powfc(float b; color e) 
{
    return color(pow(b, e[0]), pow(b, e[1]), pow(b, e[2]));
}
color powcc(color b, e) {
    return color(pow(b[0], e[0]), pow(b[1], e[1]), pow(b[2], e[2]));
}

float inputbias(float x, bias) {
    return pow(x, -log(0.5 + bias) / log(2));
}

color color2opac(color c)
{
    color h = ctransform("hsv", c);
    return ctransform("rgb", color(mod(h[0]+0.5, 1), h[1], h[2]));
}

color rgbhue(color c)
{
    color h = ctransform("hsv", c);
    return ctransform("rgb", color(h[0],h[1],1));
}

color maxfv(float a; color b) {
   return color(max(a,b[0]), max(a,b[1]), max(a,b[2]));
}

color powvv(color b,e) {
   return color(pow(b[0],e[0]), pow(b[1],e[1]), pow(b[2],e[2]));
}

color powfv(float b; color e) {
   return color(pow(b,e[0]), pow(b,e[1]), pow(b,e[2]));
}

color powvf(color b; float e) {
   return color(pow(b[0],e), pow(b[1],e), pow(b[2],e));
}

color ccSat(color col; float sat)
{
    color C = col;
    if(sat != 1.0) 
    {
        float R, G, B; 
        col = color(R, G, B);
        float csat = 1.0-sat;
        float x = csat*C_RLUM, Rx = R*x;
        float y = csat*C_GLUM, Gy = G*y;
        float z = csat*C_BLUM, Bz = B*z;
        C = color(R*(x+sat) + Gy + Bz,
                  Rx + G*(y+sat) + Bz,
                  Rx + Gy + B*(z+sat));
    }
    return C;
}

// Hue Rotation Preserving Luminance
// Rotation is normalized, so [0,1] -> [0,360 degs]
color ccHue(color col; float rot) {
        vector c=(col[0], col[1], col[2]);

   if(rot!=0.) {

        vector l = (C_RLUM,C_GLUM,C_BLUM);
   
        //set x/y rotation matrices for grey vector
        float xrs=C_1_SQRT2;
        float xrc=xrs;
        matrix mat=matrix(1,0,0,0, 0,xrc,xrs,0, 0,-xrs,xrc,0, 0,0,0,1);
        float yrs=-C_1_SQRT3;
        float yrc=C_SQRT2/C_SQRT3;
        mat*=matrix(yrc,0,-yrs,0, 0,1,0,0, yrs,0,yrc,0, 0,0,0,1);
   
        //shear space to make the luminance plane horizontal
        vector ptmp=vtransform(mat, l);  
        float dx=ptmp[0]/ptmp[2];
        float dy=ptmp[1]/ptmp[2];
        mat*=matrix(1,0,dx,0, 0,1,dy,0, 0,0,1,0, 0,0,0,1);
   
        //rotate the hue
        float angle=rot*C_PI2;
        float zrs=sin(angle);
        float zrc=cos(angle);
        mat*=matrix(zrc,zrs,0,0, -zrs,zrc,0,0, 0,0,1,0, 0,0,0,1);
   
        //unshear
        mat*=matrix(1,0,-dx,0, 0,1,-dy,0, 0,0,1,0, 0,0,0,1);
   
        //un-rotate
        mat*=matrix(yrc,0,yrs,0, 0,1,0,0, -yrs,0,yrc,0, 0,0,0,1);
        mat*=matrix(1,0,0,0, 0,xrc,-xrs,0, 0,xrs,xrc,0, 0,0,0,1);
   
        c = vtransform(mat, c);
   }

   return color(c[0], c[1], c[2]);

}

// Brightness
float ccBright(float col; float b) {
   return b!=1.0 ? col*b : col;
}
color ccBrightvf(color col; float b) {
   return b!=1.0 ? col*b : col;
}
color ccBrightvv(color col; color b) {
   return b!=CONE ? col*b : col;
}


// Contrast (pivot = 1)
float ccContrast(float val, cont) {
   return max(0,((val-1.0)*cont)+1.0);
}

color ccContrastvf(color val; float cont) {
   return maxfv(float(0.0),(((val-VONE)*cont)+VONE));
   //return max(0, (((val-VONE)*cont)+VONE));
}

color ccContrastvv(color val, cont) {
   return maxfv(0,((val-VONE)*cont)+VONE);
}


// Gamma
float ccGamma(float val, gam) {
   return pow(val,1.0/gam);
}
color ccGammavf(color val; float gam) {
   return powvf(val,1.0/gam);
}
color ccGammavv(color val, gam) {
   return powvv(val,VONE/gam);
}

// HSV Correction
color ccHSV(color col, hsv) {
   return ccHue(ccBrightvf(ccSat(col,hsv[1]),hsv[2]),hsv[0]);
}

// Full CC
color cc(color col, c_hsv, c_contrast, c_gamma, c_tint) {
   color C = col;
   color redgrn = color(0,1,1);
   if(c_hsv != redgrn)   C = ccHSV(C,c_hsv);
   if(c_contrast!=CONE) C = ccContrastvv(C,c_contrast);
   if(c_gamma!=CONE)    C = ccGammavv(C,c_gamma);
   if(c_tint!=CONE)     C *= c_tint;
   return C;
}

/* look at evaluators */
class vol_density(varying float density = 1;
                  uniform color density_color = color(1, 1, 1);
                  varying float fuel = 1;
                  uniform color fuel_color = color(1, 0.5, 0.01);
                  varying float temperature = 1;
                  uniform color temperature_color = color(1, 0.776, 0.702);
                  varying float heat = 1;
                  uniform color heat_color = color(1, 0.776, 0.702);
                  varying float burn = 1;
                  uniform color burn_color = color(1, 1, 1);
			      uniform float Km = 1;
                  uniform float samples = 4;
                  uniform float bias = 0.01;
                  uniform string envMap = "/home/alan/data/data_root/envmap/snow_machine/snow_cube.ptx";
                  uniform string envSpace = "current";
			      uniform float maxVariation = 0.1;
                  varying color fx1 = color(1, 1, 1);
                  varying color fx2 = color(1, 0.894, 0.517);
                  varying color fx3 = color(0.78, 0.27, 0);
                  varying color fx4 = color(0.67, 0.078, 0);
                  varying color fx5 = color(0.207, 0, 0.003);
                  varying color fx6 = color(0.1, 0.1, 0.1);
                  varying color fx7 = color(0.05, 0.05, 0.05);
                  varying color fx8 = color(0.01, 0.01, 0.01);
                  output varying color o_temperature = 0;
                  output varying color o_fuel = 0;
                  output varying color o_burn = 0;
                  output varying color o_heat = 0;
                  output varying color o_density = 0;
                  )
{
	public void surface(output color Ci, Oi)
	{
        float D = 0;
        float T = 0;
        float F = 0;
        float H = 0;
        float B = 0;
        color Cdiff = 0;
        point lightPos = transform("world", "current", point(11.8670253754, 21.1028060913, 10.6780443192));
        //D = density * VolumeField * Km;
        //D = Km * density;
        T = Km * temperature;
        F = Km * fuel;
        H = Km * heat;
        B = Km * burn;
        D = density * VolumeField * Km;

        // output
        o_temperature = color(T, T, T) * temperature_color;
        o_fuel = color(F, F, F) * fuel_color;
        o_heat = color(H, H, H) * heat_color;
        o_burn = color(B, B, B) * burn_color;

        color fire = spline(heat, fx8, fx7, fx6, fx5, fx4, fx3, fx2, fx5, fx6, fx7);
        
        color black = color(0,0,0);
        //fire = mix(fire, color(density,density,density), heat);
        //fire *= color(T, T, T);
        color env = 0;
        /*color indirect = indirectdiffuse(P, normalize(N),*/
                                         /*samples, "pointbased", 1,*/
                                        /*"environmentmap", envMap,*/
                                        /*"environmentspace", envSpace,*/
                                        /*"environmentcolor", env,*/
                                         /*"maxvariation", maxVariation,*/
                                         /*"clamp", 1, "volume", 1);*/


        illuminance(P) {
            Cdiff += Cl;
        }

        //Cdiff += indirect;
        Cdiff += fire;
        //Cdiff += o_fuel;
        //Cdiff += o_temperature;
        //Cdiff += o_heat;
        //Cdiff += o_burn;
        /*if (temperature != 1) {*/
            /*Cdiff *= color(1, 0, 0);*/
        /*}*/

        float Opac = transmission(P, lightPos)[0];

		Oi = D;
        Ci = Cs * Cdiff;
		Ci = Ci * Oi *Opac;
	}
}
