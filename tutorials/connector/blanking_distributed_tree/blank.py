#!/usr/bin/env python
# coding: utf-8

r"""Blanking en distribue

This tutorial demonstrates how to blank a distributed tree, that is a pyTree
with zones distributed on different processors.

Blanking a distributed pyTree, by C. Benoit.

The first script (case.py) is only used to setup the case. Second script
performs blanking.
It is executed on each processor and must be launched
with mpirun -np 10 python blank.py. It loads walls.cgns on every processors.
This file contains the definition of bodies. Then it loads only the skeleton
of the pyTree to be blanked.
A distribution is performed to determine which
zone must be loaded on each processor.
Since blanking doesnt require any communications once the walls are defined,
useCom=0 and fast algorithm are chosen.
Then, only certain zones are loaded on processor according to the
distribution using C.readZones. Skeleton tree is then converted to a partial
tree, that is a tree where only loaded zones are defined.
Blanking is performed on local zones.
Finally, each processor write its zones to a common file
using Cmpi.convertPyTree2File.

"""

import Converter.PyTree as C
import Distributor2.PyTree as Distributor2
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Connector.PyTree as X
import Converter.Internal as Internal
import numpy

rank = Cmpi.rank
size = Cmpi.size

# lecture des corps servant a masquer
bodies = C.convertFile2PyTree('walls.cgns')
bodies = Internal.getNodesFromType(bodies, 'Zone_t')

# lecture du squelette
a = Cmpi.convertFile2SkeletonTree('in.cgns')

# equilibrage
(a, dic) = Distributor2.distribute(a, NProc=size, algorithm='fast', useCom=0)

# load des zones locales dans le squelette
a = Cmpi.readZones(a, 'in.cgns', proc=rank)

# Passage en arbre partiel
a = Cmpi.convert2PartialTree(a)

# Blanking local
BM = numpy.array([[1]])
a = X.blankCells(a, [bodies], BM, blankingType='node_in', dim=3)

# Reconstruit l'arbre complet a l'ecriture
Cmpi.convertPyTree2File(a, 'out.cgns')
