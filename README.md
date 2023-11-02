========
open-pic
========

open-pic is a 3D [Particle-In-Cell (PIC) code](https://en.wikipedia.org/wiki/Particle-in-cell)
for collisionless plasma simulation.
Depends on luabind and boost libraries.


Build steps for Ubuntu 18:

Install lua
   sudo apt-get install lua5.1

1. Extract /3rdparty/3rdparty.tar.xz ti the same directory
2. In Makefile and setup.sh set the OPIC_DIR to your local repo directory path
3. Run
      make
   in your local repo directory path
4. Make sure that the "opic" binary eecutable is created.
5. Copy "opic" file to /sim/m15d1 
      cp ./opic ./sim/m15d1/
6. Run comands listed setup.sh commands in your terminal
7. Run
    cd ./sim/m15d1
8. Run open-pic sumulation
    opic main.lua


See "main.lua" and other Lua scripts for simulation flow control.

General information about open-pic can be found on comphys.narod.ru.
