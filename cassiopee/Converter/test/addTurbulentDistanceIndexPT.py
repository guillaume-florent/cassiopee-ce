# - distance2Walls (pyTree) -
import Dist2Walls.PyTree as Dist2Walls
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import Converter.elsAProfile as elsAProfile

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
sphere = D.sphere((1.2,0.,0.),0.2,100)
t = C.newPyTree(['Base',a])
t = Dist2Walls.distance2Walls(t, sphere)
Internal._renameNode(t,'FlowSolution#Centers','FlowSolution#Init')
tp = elsAProfile.addTurbulentDistanceIndex(t)
C.convertPyTree2File(tp, 'out.cgns')
