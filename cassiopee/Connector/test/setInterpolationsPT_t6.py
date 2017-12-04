import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import Converter.Internal as Internal
import Transform.PyTree as T
import KCore.test as test

LOCAL = test.getLocal()

NK = 2
DIM = 3
if NK == 2: DIM = 2
a = G.cylinder( (0,0,0), 1, 2, 10, 130, 1, (60, 20, NK)); a[0] = 'cyl1'
b = G.cylinder( (0,0,0), 1, 1.5, 0, 30, 1, (30, 20, NK)); b[0] = 'cyl2'
a = C.addBC2Zone(a, 'wall', 'BCWall', 'jmin')
b = C.addBC2Zone(b, 'wall', 'BCWall', 'jmin')
b = C.fillEmptyBCWith(b,'overlap','BCOverlap',dim=DIM)
t = C.newPyTree(['Base','Base2']); t[2][1][2] = [a]; t[2][2][2] = [b]
t = C.initVars(t,'Density',1.); t = C.initVars(t,'centers:cellN',1)
t = X.connectMatchPeriodic(t,rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,120.], dim=DIM)
t = C.fillEmptyBCWith(t,'nref','BCFarfield',dim=DIM)
t[2][1] = C.addState(t[2][1], 'EquationDimension',DIM)
t = X.applyBCOverlaps(t,depth=1)
t = X.setInterpolations(t,double_wall=1,storage='direct',prefixFile=LOCAL+'/chm')
test.testT(t,1)
# test avec Chimere periodique
NK = 51
a = G.cylinder( (0,0,0), 1, 2, 10, 130, 4., (60, 20, NK)); a[0] = 'cyl1'
b = G.cart((0.4,1.2,-0.3),(0.04,0.04,0.1),(11,11,21))
a = X.connectMatchPeriodic(a,rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,120.])
a = C.fillEmptyBCWith(a,"wall","BCWall")
a = C.addBC2Zone(a,'overlapdd','BCOverlap','kmin',zoneDonor=[b],rangeDonor='doubly_defined')#dd
# 
b = C.addBC2Zone(b,'overlap','BCOverlap','kmax')
b = C.fillEmptyBCWith(b,"wall","BCWall")
for rangel in ['imin','imax','jmin','jmax']:
    b = C.addBC2Zone(b,'overlapdd','BCOverlap',rangel,zoneDonor=[a],rangeDonor='doubly_defined')
Internal.createChild(b,'.Solver#Param','UserDefinedData_t',value=None,children=[])
solverParam = b[2][len(b[2])-1]
Internal.createChild(solverParam,'axis_ang_1'  ,'DataArray_t',value=3,children=[])
Internal.createChild(solverParam,'axis_ang_2'  ,'DataArray_t',value=1,children=[])
Internal.createChild(solverParam,'axis_pnt_x'  ,'DataArray_t',value=0.,children=[])
Internal.createChild(solverParam,'axis_pnt_y'  ,'DataArray_t',value=0.,children=[])
Internal.createChild(solverParam,'axis_pnt_z'  ,'DataArray_t',value=0.,children=[])
Internal.createChild(solverParam,'axis_vct_x'  ,'DataArray_t',value=0.,children=[])
Internal.createChild(solverParam,'axis_vct_y'  ,'DataArray_t',value=0.,children=[])
Internal.createChild(solverParam,'axis_vct_z'  ,'DataArray_t',value=1.,children=[])
Internal.createChild(solverParam,'periodic_dir' ,'DataArray_t',value=3,children=[])
b = T.rotate(b,(0,0,0),(0,0,1),60.)
#
t = C.newPyTree(['Base','Base2']); t[2][1][2] = [a]; t[2][2][2] = [b]
t = C.initVars(t,'Density',1.); t = C.initVars(t,'centers:cellN',1)
t = X.applyBCOverlaps(t,depth=1)
t = X.setDoublyDefinedBC(t,depth=1)
t = X.setInterpolations(t,double_wall=1,storage='direct',prefixFile=LOCAL+'/chm')
test.testT(t,2)
