#!/usr/bin/env python
# coding: utf-8

r"""Tutorial 1 for Cassiopee

This tutorial shows how to run Cassiopee scripts using python.

This tutorial supposed that Cassiopee is installed in the directory CASSIOPEE.
Environement
Under Windows: In the windows menu, open the Cassiopee folder.
Click on the "Command Shell" of Cassiopee. This will open a dos interpreter
with the ready environment to run Cassiopee scripts.
If you don't want to use this short cut, you can load by yourself the file
CASSIOPEE/Dist/env_Cassiopee.bat depending on the version you have downloaded.

Cassiopee can take advantage of multi cores using (here 4 cores):
set OMP_NUM_THREADS=4

Running python script is done with:
python tuto1.py

Under Unix shell (sh/bash/ksh): Loading environment is done with:
. $CASSIOPEE/Dist/sh_Cassiopee

Cassiopee can take advantage of multi cores using (here 4 cores):
export OMP_NUM_THREADS=4

Running python script is done with:
python tuto1.py

If you use mpi in your script (C.mpi, etc...), you must run:
mpirun -np 2 python tuto.py

Under Unix C-shell (csh/tcsh): Loading environment is done with:
source $CASSIOPEE/Dist/env_Cassiopee

Cassiopee can take advantage of multi cores using: setenv OMP_NUM_THREADS=4

Running python script is done with:
python tuto1.py

If you use mpi in your script (C.mpi, etc...), you must run:
mpirun -np 2 python tuto.py

"""

from __future__ import print_function

# import Converter.PyTree as C
import Generator.PyTree as G


a = G.cart(Xo=(0, 0, 0), H=(1, 1, 1), N=(10, 10, 10))
print("Creating a cartesian grid.")
