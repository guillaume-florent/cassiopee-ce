# - close (pyTree)-
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# test structure
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0.01, 10., (20,20,10)) 
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a, 'wall','BCWall')
a = C.addVars(a, 'Density'); a = C.initVars(a, 'centers:cellN', 1)
a1 = G.close(a, 1.e-3)
test.testT(a1, 1)

# test non structure hexa
a2 = C.convertArray2Hexa(a)
a2 = G.close(a2, 1.e-3)
test.testT(a2, 2)

# test close non structure tetra
a3 = C.convertArray2Tetra(a)
a3 = G.close(a3, 1.e-3)
test.testT(a3, 3)
