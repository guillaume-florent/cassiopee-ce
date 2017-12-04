import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

t1 = C.newPyTree(['Base']); t2 = C.newPyTree(['Base'])
z1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
t1[2][1][2].append(z1); t2[2][1][2].append(z1)
t1 = C.initVars(t1,'centers:cellN',1.)
t2 = C.initVars(t2,'centers:cellN',0.)
t1 = C.initVars(t1,'centers:Density',1.)
t1 = C.initVars(t1,'Pressure',10.)

t2 = P.importVariables(t1, t2)
C.convertPyTree2File(t2, 'out.cgns')
