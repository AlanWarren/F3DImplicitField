Field3d Implicit Field plugin for Pixar's Renderman
==========================================================

F3DImplicitField.cpp is the Implicit Field DSO.

f3d2prman.cpp is a runprogram that facilitates the use of F3DImplicitField. 
You will use it to generate the RIB that calls our DSO.

Notes
-----

It's been a good little while since I wrote this (2012), so please excuse any missing information. I did blog about the project when it was originally developed.

(archived to markdown since my site is no longer around)
1. [Rendering Volumes 1](Rendering_Volumes1.md)

2. [Rendering Volumes 2](Rendering_Volumes2.md)

3. [Rendering Volumes 3](Rendering_Volumes3.md)

This has only been built and tested on LINUX using prman 17+18 with Houdini. If you're trying to build it on another platform, or with another version of prman good luck!

I haven't touched this code since 2012, and have no plans to maintain it with current versions of prman or field3d. This was a learning project for me, and the
main reason I'm putting on github, is to hopefully help somebody else who's looking to learn from it. Feel free to contact me if you have any questions. 

Requirements
------------

* Field3D
* OpenEXR (ilmbase)
* Boost
* Pixar's Renderman Pro Server

Few compilation notes: 
* I compiled HDF5 with +static-libs and -fortran -cxx.  
* For openmpi, I used the default config, except for +romio and +cxx. 
* Boost was compiled with +static-libs +mpi 
* ilmbase was also compiled with +static-libs

Build
----

I'm not providing a makefile, but this should build and link.

g++ -fPIC -I$RMANTREE/include -c

g++ F3DImplicitField.o -lField3D -L$RMANTREE/lib -lprman -shared -o F3DImplicitField.so
