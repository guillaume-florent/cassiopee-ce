# - patch (pyTree) -
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Converter.Internal as Internal
import KCore.test as test

def dens(x,y): return 3*x*y

# cas 2D
c1 = G.cart((0,0,0), (0.01,0.01,1), (201,101,2))
c1 = C.addBC2Zone(c1,'wall1','BCWall','imin')
c1 = C.addBC2Zone(c1,'match1','BCMatch','imax',c1,'imin',[1,2,3])
c1 = C.addBC2Zone(c1, 'overlap1', 'BCOverlap', 'jmax')
c1 = C.initVars(c1,'centers:celln',1.)
c1 = C.initVars(c1,'Density',dens,['CoordinateX','CoordinateY'])

c2 = G.cart((0,0,0), (0.01,0.01,1), (51,81,2))
c2 = T.rotate(c2, (0,0,0),(0,0,1),0.2)
c2 = C.addBC2Zone(c2,'wall1','BCWall','imin')
c2 = C.addBC2Zone(c2,'overlap1','BCOverlap','imax')
c2 = C.initVars(c2, 'centers:celln',1.)
c2 = C.initVars(c2,'Density',dens,['CoordinateX','CoordinateY'])
a = T.patch(c1, c2,(1,1,1))
t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
test.testT(t, 1)

# multizone
a = D.sphere6((0,0,0), 1, N=20)
t = C.newPyTree(['Base',3]); t[2][1][2] += a
k=0
t2 = T.subzone(t,(1,1,k+1),(-1,-1,k+1))
t2 = T.smooth(t2, eps=0.5, niter=20, projConstraints=Internal.getNodesFromType(t2,'Zone_t'))
t = T.patch(t, t2, (1,1,k+1))
test.testT(t, 3)

# cas 3D
c1 = G.cart((0,0,0), (0.01,0.01,1), (201,101,20))
c1 = C.addBC2Zone(c1,'wall1','BCWall','imin')
c1 = C.addBC2Zone(c1,'match1','BCMatch','imax',c1,'imin',[1,2,3])
c1 = C.addBC2Zone(c1, 'overlap1', 'BCOverlap', 'jmax')
c1 = C.initVars(c1,'centers:celln',1.)
c1 = C.initVars(c1,'Density',dens,['CoordinateX','CoordinateY'])

c2 = G.cart((0,0,0), (0.01,0.01,1), (51,81,20))
c2 = T.rotate(c2, (0,0,0),(0,0,1),0.2)
c2 = C.addBC2Zone(c2,'wall1','BCWall','imin')
c2 = C.addBC2Zone(c2,'overlap1','BCOverlap','imax')
c2 = C.initVars(c2, 'centers:celln',1.)
c2 = C.initVars(c2,'Density',dens,['CoordinateX','CoordinateY'])
a = T.patch(c1, c2,(1,1,1))
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t, 2)

# multi-zone 3D
a = D.sphere6((0,0,0), 1, N=20,)
d = G.cart((0.1,0.,0.), (0.1,1,1),(5,1,1)) 
a = G.addNormalLayers(a, d)
t = C.newPyTree(['Base',3]); t[2][1][2] += a
t2 = T.subzone(t,(1,1,2),(-1,-1,4))
t2 = T.smooth(t2, eps=0.5, niter=20)
t = T.patch(t,t2,(1,1,2))
test.testT(t, 4)
