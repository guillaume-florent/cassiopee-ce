# - connectMatch (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11)) 
t = C.newPyTree(['Base',a])
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],
                           translation=[0.,0.,5.])
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],
                           rotationAngle=[0.,0.,90.])
tp = elsAProfile.adaptPeriodicMatch(t)
test.testT(tp, 1)
