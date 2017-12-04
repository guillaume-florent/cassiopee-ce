# - integMoment (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

m = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
m = C.initVars(m,'vx',1.); m = C.initVars(m,'vy',0.); m = C.initVars(m,'vz',0.)
res = P.integMoment(m, center=(5.,5., 0.),vector=['vx','vy','vz']); print res
