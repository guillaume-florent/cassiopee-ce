#!/usr/bin/env python
# coding: utf-8

r"""Computes the matching connectivity for a portion of cylinder

This tutorial demonstrates how to compute the periodic match connectivity
in the case of an azimuthal periodicity.

Computes the periodic match connectivity of a quarter of cylinder by S. Péron.

Suppose you have a periodic mesh in azymuth and want to compute the
periodic match connectivity.
The idea is to duplicate the meshes by rotation and then perform
a X.connectMatch on the full mesh. Since zone names in the duplicated mesh
are the same as in original mesh, this will result in a correct connectivity.
Then, we perform the same thing for the other direction.

"""

import Converter.PyTree as C
import Connector.PyTree as X
import Transform.PyTree as T
import Generator.PyTree as G

a = G.cylinder((0., 0., 0.), 0.1, 1., 0., 90., 5., (11, 11, 11))
t = C.newPyTree(['Cylindre'])
t[2][1][2] += [a]

# Duplicate mesh by a rotation of +90 degrees
ap = T.rotate(a, (0., 0., 0.), (0., 0., 1.), 90.)
t = C.addBase2PyTree(t, 'CylindreDup')
t[2][2][2] += [ap]
t = X.connectMatch(t)

# Duplicate mesh by a rotation of +90 degrees
# Replace duplicated zones
am = T.rotate(a, (0., 0., 0.), (0., 0., 1.), -90.)
t[2][2][2] = [am]
t = X.connectMatch(t)

# Remove duplicated basis
t[2][1:] = t[2][1:2]
C.convertPyTree2File(t, 'out.cgns')
