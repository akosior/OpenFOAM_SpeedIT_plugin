                                  INTRODUCTION

This code allows to integrate CUDA accelerated linear systems of equations 
solvers from SpeedIT Classic toolkit with OpenFOAM 1.6. The code is ditributed
under terms of the GNU General Public License as published by the Free Software
 Foundation, either version 3 of the License, or (at your option) any later
 version.


                                PLUGIN INSTALATION

1. Go to your path for FOAM applications

  cd $WM_PROJECT_USER_DIR


2. Download plugin source

  Plugin source can be downloaded from svn
  
  svn checkout https://62.87.249.40/repos/speedIT/branches/1.0/OpenFOAM_SpeedIT_plugin  


3. Go into source directory

  cd OpenFOAM_SpeedIT_plugin/


4. Edit solver source files 

  PBiCG_accel.C
  PCG_accel.C

to allow use of your solver function. By default functions from SpeedIT Classic
library are used. If you wish to use other solver functions, remember to copy 
or link apropriate header files to 

 $WM_PROJECT_USER_DIR/OpenFOAM_SpeedIT_plugin/

directory OR to add path to your header files to file

 $WM_PROJECT_USER_DIR/OpenFOAM_SpeedIT_plugin/Make/options


5. If your CUDA installation is NOT in default path, update path to CUDA header
files in file

  $WM_PROJECT_USER_DIR/OpenFOAM_SpeedIT_plugin/Make/options


6. Compile

  wmake libso

After this step you should have library libexternalsolv.so in directory 
$FOAM_USER_LIBBIN



                                PLUGIN USE

1. Copy libraries 

  libSpeedIT_Classic.so
  libcublas.so
  libcudart.so

to $FOAM_USER_LIBBIN directory. After proper plugin installation you should
already have library libexternalsolv.so in $FOAM_USER_LIBBIN directory. If you
use solver functions from other library than libSpeedIT_Classic.so, copy your
library instead of libSpeedIT_Classic.so to $FOAM_USER_LIBBIN directory.

NOTE: libcublas.so and libcudart.so are parts of your CUDA installation. 
Remember to copy proper library version: 32-bit library in 32-bit operating
systems and 64-bit library for 64-bit operating systems.


2. Go in to the directiory with your FOAM case, i.e.

  cd $FOAM_RUN/tutorials/incompressible/icoFoam/cavity


3. Append


libs
    (
        "libcudart.so"
        "libcublas.so"
        "libSpeedIT_Classic.so"
        "libexternalsolv.so"
    );


to the end of your system/controDict file for every FOAM case, for which you
want to use external, accelerated solvers.


4. In file system/fvSolution change solver names for solvers, for which you
are going to enable acceleration. Remember to use proper names for accelerated
solvers. You may replace:

  PBiCG with PBiCG_accel
  PCG   with PCG_accel


5. For accelerated solvers choose an apriopriate preconditioner in file 
system/fvSolution. If you use unsupported preconditioner, warning will be
printed and no preconditioner will be used.


6. Run FOAM solver, i.e. icoFoam. Accelerated solvers should be used by now.


