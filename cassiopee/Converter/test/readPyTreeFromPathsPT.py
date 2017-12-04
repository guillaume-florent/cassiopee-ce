# - readPyTreeFromPaths (pyTree) -
import Converter.PyTree as C
import Converter.converter as converter
import Converter.Internal as Internal
import Generator.PyTree as G

# Cree le fichier test
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, 'test.hdf')


# Relit les noeuds par leur paths
nodes = converter.readPyTreeFromPaths('test.hdf', ['/Base/cart/GridCoordinates'], 'bin_hdf')
Internal.printTree(nodes)
