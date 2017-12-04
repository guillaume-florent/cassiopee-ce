# - identifySolutions (array) -
import Converter.PyTree  as C
import Post.PyTree  as P
import Generator.PyTree  as G
import Geom.PyTree  as D
import KCore.test as test
ni = 21; nj = 21; nk = 21
m = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))
m = C.initVars(m, 'G={CoordinateY}')
m = C.initVars(m,'centers:F={centers:CoordinateX}')
m = C.fillEmptyBCWith(m, 'nref','BCFarfield')

# Create receptor mesh
a = D.sphere((0,0,0),0.1)
# Structure 
hookN = C.createGlobalHook([m],'nodes')
hookC = None
a2 = C.identifySolutions(a,m,hookN,hookC,tol=1000.)
test.testT([a2],1)

# TRI
a1 = C.convertArray2Tetra(a)
a2 = C.identifySolutions(a1,m,hookN,hookC,tol=1000.)
test.testT([a2],2)

# NGON
a1 = C.convertArray2NGon(a)
a2 = C.identifySolutions(a1,m,hookN,hookC,tol=1000.)
test.testT([a2],3)
