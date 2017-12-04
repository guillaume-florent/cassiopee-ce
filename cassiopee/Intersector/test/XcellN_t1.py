# - XcellN (array) -
import Converter as C
import Connector as X
import Generator as G
import Geom as D
import Intersector as XOR
import KCore.test as test
import sys

# Test 1
# Mask
masking = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
masking = C.convertArray2NGon(masking)
# Mesh to blank
a = G.cart((-3.,-3.,-3.), (0.5,0.5,0.5), (20,20,20))
a = C.convertArray2NGon(a)

# celln init
ca = C.node2Center(a)
ca = C.initVars(ca, 'cellN', 1.)
ca = C.extractVars(ca, ['cellN'])

# Blanking
celln = XOR.XcellN([a], [ca], masking)
test.testA(celln,1)



