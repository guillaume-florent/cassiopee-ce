# - nearestFaces (array) -
# doit etre exact
import Converter as C
import Generator as G
import Transform as T
import time
eps = 1.e-45
N = 50
# structure
a = G.cart((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cart((100000,0,0), (1,1.e-8,1), (N,N,N))
hook = C.createHook(a, function='faceCenters')
faces = C.identifyFaces(hook,b,tol=eps)
ret = faces[faces<0].shape[0]
if ret>0: print 'identifyFaces (Struct) FAILED: FPU is not correct [Struct/FACES]. Check compilation options.'
# HEXA
a = G.cartHexa((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cartNGon((100000,0,0), (1,1.e-8,1), (N/2,N/2,N/2))
hook = C.createHook(a, function='faceCenters')
faces = C.identifyFaces(hook,b,tol=eps)
ret = faces[faces<0].shape[0]
if ret>0: print 'identifyFaces (BE) FAILED: FPU is not correct [BE/FACES]. Check compilation options.'
# NGON
a = G.cartNGon((100000,0,0), (1,1.e-8,1), (N,N,N))
b = G.cartNGon((100000,0,0), (1,1.e-8,1), (N/2,N/2,2))
hook = C.createHook(a, function='faceCenters')
faces = C.identifyFaces(hook,b,tol=eps)
ret = faces[faces<0].shape[0]
if ret >0: print 'identifyFaces (NGON) FAILED: FPU is not correct [NGON/FACES]. Check compilation options.'
print 'done.'
