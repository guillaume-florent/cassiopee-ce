# - addReferenceState (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base']); t[2][1][2].append(a)
tp = elsAProfile.addReferenceState(t, conservative=[1.,0,0,0,1.7])
test.testT(tp, 1)
