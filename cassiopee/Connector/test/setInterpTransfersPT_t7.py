# - setInterpTransfers IBC + extraction (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import Post.PyTree as P
import numpy as N
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Initiator.PyTree as I
import Converter.Internal as Internal
import KCore.test as test

a = G.cart((-1,-1,-1),(0.04,0.04,1),(51,51,3))
s = G.cylinder((0,0,-1), 0, 0.4, 360, 0, 4, (30,30,5)) 
s = C.convertArray2Tetra(s); s = T.join(s); s = P.exteriorFaces(s)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
# Blanking
bodies = [[s]]
BM = N.array([[1]],N.int32)
t = X.blankCells(t,bodies,BM,blankingType='center_in')
t = X.setHoleInterpolatedPoints(t,depth=-2)
# Dist2Walls
t = DTW.distance2Walls(t,[s],type='ortho',loc='centers',signed=1)
t = C.center2Node(t,'centers:TurbulentDistance')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
t = I.initConst(t,MInf=0.2,loc='centers')
tc = C.node2Center(t)
t = X.setIBCData(t, tc, loc='centers', storage='direct')
no = 1
for bcType in xrange(4):
    for varType in xrange(1,4):
        t2 = X.setInterpTransfers(t,tc,bcType=bcType,varType=varType)
        test.testT(t2,no)
        no+=1

# Stockage inverse
Internal._rmNodesByName(t,"IBCD_*")
tc = X.setIBCData(t, tc, loc='centers', storage='inverse')
for bcType in xrange(4):
    for varType in xrange(1,4):
        X._setInterpTransfers(t,tc,bcType=bcType,varType=varType)
        test.testT(tc,no)
        no += 1
