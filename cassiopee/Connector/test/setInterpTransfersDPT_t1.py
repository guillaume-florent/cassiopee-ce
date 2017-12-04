# - setInterpTransfers (pyTree)-
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder( (0,0,0), 1, 2, 0, 360, 1, (60, 20, 3) )
b = G.cylinder( (0,0,0), 1, 2, 3, 160, 1, (30, 10, 3) )
a = C.addBC2Zone(a, 'wall', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match', 'BCMatch', 'imin', a, 'imax', trirac=[1,2,3])
a = C.addBC2Zone(a, 'match', 'BCMatch', 'imax', a, 'imin', trirac=[1,2,3])
b = C.addBC2Zone(b, 'wall', 'BCWall', 'jmin')
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imin')
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imax')
tD = C.newPyTree(['Base']); tD[2][1][2] = [a]; 
tR = C.newPyTree(['Base']); tR[2][1][2] = [b]
tD = C.fillEmptyBCWith(tD, 'nref', 'BCFarfield')
tR = C.fillEmptyBCWith(tR, 'nref', 'BCFarfield')

tD = C.initVars(tD, '{Density}=1.')
tD = C.initVars(tD, '{cellN}=1.')
tD = C.initVars(tD, '{MomentumX}= -0.1')
tD = C.initVars(tD, '{MomentumY}= -0.2')
tR = X.applyBCOverlaps(tR, depth=1) 
tD = X.setInterpData(tR, tD, double_wall=1, loc='centers',
                     storage='inverse', order=3)
info = X.setInterpTransfersD(tD, variables=['MomentumX'])
test.testO(info,1)
test.testA([info[0][1]],11)

# Noeuds
tD = C.newPyTree(['Base']); tD[2][1][2] = [a]; 
tR = C.newPyTree(['Base']); tR[2][1][2] = [b]
tD = C.fillEmptyBCWith(tD, 'nref', 'BCFarfield')
tR = C.fillEmptyBCWith(tR, 'nref', 'BCFarfield')

tD = C.initVars(tD, '{Density}=1.')
tD = C.initVars(tD, '{cellN}=1.')
tD = C.initVars(tD, '{MomentumX}= -0.1')
tD = C.initVars(tD, '{MomentumY}= -0.2')
tR = X.applyBCOverlaps(tR, depth=1,loc='nodes') 
tD = X.setInterpData(tR, tD, double_wall=1, loc='nodes',
                     storage='inverse', order=2)
info = X.setInterpTransfersD(tD,variables=['MomentumX','Density'])
test.testO(info,2)
test.testA([info[0][1]],21)
