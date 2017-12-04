# - copyNode (pyTree) - 
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = Internal.newDataArray(value=[1,2,3])
# Copy only the numpy of this node
b = Internal.copyNode(a)
# Modify numpy of b
b[1][0]=5
# a is not modified
print a
#>> ['Data', array([1, 2, 3], dtype=int32), [], 'DataArray_t']
print b
#>> ['Data', array([5, 2, 3], dtype=int32), [], 'DataArray_t']
