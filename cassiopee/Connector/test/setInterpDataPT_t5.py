# - setInterpData (pyTree)-
# cas structure double wall 
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder( (0,0,0), 1, 2, 0, 360, 1, (60, 20, 3) )
b = G.cylinder( (0,0,0), 1, 2, 3, 160, 1, (30, 20, 3) )
a = C.addBC2Zone(a, 'wall', 'FamilySpecified:SKN1', 'jmin') 
a = C.addBC2Zone(a, 'match', 'BCMatch', 'imin', a, 'imax', trirac=[1,2,3])
a = C.addBC2Zone(a, 'match', 'BCMatch', 'imax', a, 'imin', trirac=[1,2,3])
b = C.addBC2Zone(b, 'wall', 'BCWall', 'jmin')
b = C.addBC2Zone(b, 'wall', 'FamilySpecified:SKN2', 'jmin') 
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imin')
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imax')
t = C.newPyTree(['Base','Base2']); t[2][1][2] = [a]; t[2][2][2] = [b]
t = C.fillEmptyBCWith(t,'nref','BCFarfield')
t = C.initVars(t,'Density',1.); t = C.initVars(t,'centers:G',10.)
t[2][1] = C.addState(t[2][1], 'EquationDimension',2)
t[2][1] = C.addFamily2Base(t[2][1], 'SKN1', bndType="BCWall") 
t[2][2] = C.addFamily2Base(t[2][2], 'SKN2', bndType="BCWall") 

t1 = X.applyBCOverlaps(t, depth=2) 
t1[2][2] = X.setInterpData(t1[2][2],t1[2][1], double_wall=1,loc='centers')
test.testT(t1,1)
