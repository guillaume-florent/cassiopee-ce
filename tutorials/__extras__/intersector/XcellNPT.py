# - XcellN (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Intersector.PyTree as XOR
import Geom.PyTree as D

# Test 1
# Mask
masking = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
masking = C.convertArray2NGon(masking)
# Mesh to blank
bgm = G.cart((-3.,-3.,-3.), (0.5,0.5,0.5), (20,20,20))
t = C.newPyTree(['Cart', bgm])
t = C.convertArray2NGon(t)

# celln init
t = C.initVars(t, 'centers:cellN', 1.)
# Blanking with floating cellN computation
t = XOR.XcellN(t, [[masking]], [])
C.convertPyTree2File(t, 'out1.cgns')

# Test 2
# Tet mask
masking = D.sphere((0,0,0), 15., 30)
masking = C.convertArray2Tetra(masking)
masking = G.close(masking)
masking = G.tetraMesher(masking, algo=1)
#C.convertPyTree2File(masking, 'sph.cgns')
# Mesh to blank
bgm = G.cart((-5.,-5.,-5.), (0.8,0.8,0.8), (40,40,40))
t = C.newPyTree(['Cart', bgm])
t = C.convertArray2NGon(t)
# celln init
t = C.initVars(t, 'centers:cellN', 1.)
# Blanking
t = XOR.XcellN(t, [[masking]], [])
C.convertPyTree2File(t, 'out2.cgns')
