#!/usr/bin/python
# coding: utf-8

r"""Generator : cut an octree along X/Y/Z planes"""

# cut an octree along X/Y/Z planes
#
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import Geom.PyTree as D
import Post.PyTree as P
import Converter.Internal as Internal

# dictionary of the X/Y/Z planes
dico = {'CoordinateX': 2., 'CoordinateZ': 1.}


a = D.sphere((0, 0, 0), 1., 20)
a = C.convertArray2Tetra(a)
a = G.close(a)
surf = C.newPyTree(['SURF', 2])
surf[2][1][2].append(a)

snear = [0.5]

# symetrize the bodies wrt. X/Y/Z planes
for coord in dico.keys():
    for zone in Internal.getByType(surf, 'Zone_t')[2]:
    if coord=='CoordinateX':
        surf[2][1][2].append(T.symetrize(zone,
                                         (dico[coord], 0., 0.),
                                         (0, 1, 0),
                                         (0, 0, 1)))
    if coord=='CoordinateY':
        surf[2][1][2].append(T.symetrize(zone,
                                         (0., dico[coord], 0.),
                                         (1, 0, 0),
                                         (0, 0, 1)))
    if coord=='CoordinateZ':
        surf[2][1][2].append(T.symetrize(zone,
                                         (0., 0., dico[coord]),
                                         (1, 0, 0),
                                         (0, 1, 0)))
    snear+=snear

octr = G.octree(surf, snear, dfar=5.)

# select cells below X/Y/Z planes (can be modified for cells beyond planes)
for coord in dico.keys():
    octr = P.selectCells(octr, '{%s}<=%s'%(coord,dico[coord]), strict=1)

octr=G.octree2Struct(octr,
                     vmin=11,
                     ext=0,
                     optimized=0,
                     merged=1,
                     sizeMax=250000)
carttree=C.newPyTree(['Cart', 3])
carttree[2][1][2]+=octr

C.convertPyTree2File(carttree,'cart.cgns')
