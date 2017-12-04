# - scale (pyTree) -
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,10))

# scale in all directions
T._scale(a, factor=0.1)

# scale with different factors following directions
T._scale(a, factor=(0.1,0.2,0.3))

C.convertPyTree2File(a, 'out.cgns')
