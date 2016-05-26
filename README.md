Field3d Implicit Field plugin for Pixar's Renderman
==========================================================

F3DImplicitField is the DSO

f3d2prman is a runprogram that facilitates the use of F3DImplicitField.so You will use it to generate the RIB that calls our DSO. 

Notes
-----

It's been a good little while since I wrote this (2012), so please excuse any missing information. I did blog about the project when it was originally developed.
1. http://www.alan-warren.com/2012/07/26/539/ 
2. http://www.alan-warren.com/2012/08/15/580/  
3. http://www.alan-warren.com/2012/08/30/586/ 

This has only been built and tested on LINUX. If you're trying to build it on another platform, good luck! 

Requirements
------------

Field3D 
Boost
Pixar's Renderman Pro Server

Few compilation notes: 
I compiled HDF5 with +static-libs and -fortran -cxx.  
For openmpi, I used the default config, except for +romio and +cxx. 
Boost was compiled with +static-libs +mpi 
ilmbase was also compiled with +static-libs
