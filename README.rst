cassiopee
=========

CFD Advanced Set of Services in an OpenPython Environment

Install
-------

Ubuntu 14.04

Requirements: you will need dev version of Python, numpy, scons, gcc, gfortran and X11. If you don't have them:

sudo apt-get install python-dev
sudo apt-get install python-numpy
sudo apt-get install scons
sudo apt-get install gfortran
sudo apt-get install xorg-dev

Then, first install KCore module needed by all modules:

gunzip KCore-X.X.tar.gz
tar xvf KCore-X.X.tar
cd KCore

You have to set a few things in the environment:
export CASSIOPEE=/home/guillaume/_Repositories/github/guillaume-florent/cassiopee/build/
If you are using a 32 version of ubuntu:
export LD_LIBRARY_PATH=/usr/lib/i386-linux-gnu
If you are using a 64 version of ubuntu:
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu

Then type:
./install

Installation will end with Failed (dont worry). Then:
source $CASSIOPEE/Dist/env_Cassiopee.sh

./install

It should now be correctly installed. Proceed identically with other modules.
This procedure may also work for other debian linux flavours. 
