# - deform (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

# Zone structuree
a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,2))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a = C.addVars(a, 'F'); a = C.addVars(a, 'centers:G')
vect = ['hx','hy','hz'];a = C.addVars(a, vect)
a = C.initVars(a, 'hx', 10.)
a = T.deform(a, vect)
test.testT(a, 1)

# Zone non structuree
a = G.cartTetra((0.,0.,0.),(1.,1.,1.),(10,10,2))
a = C.addVars(a, 'F'); a = C.addVars(a, 'centers:G')
vect = ['hx','hy','hz']; a = C.addVars(a, vect)
a = C.initVars(a, 'hx', 10.)
a = T.deform(a, vect)
test.testT(a, 2)

# Arbre
a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,1))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
b = G.cartTetra((1.5,1.5,0.),(1.,1.,1.),(10,10,1)); b[0] = 'cart2'
t = C.newPyTree(['Base',2]); t[2][1][2] += [a, b]
t = C.initVars(t, 'F',2); t = C.addVars(t, 'centers:G')
vect = ['hx','hy','hz']; t = C.addVars(t, vect)
t = C.initVars(t, 'hx', 10.)
t = T.deform(t, vect)
test.testT(t, 3)
