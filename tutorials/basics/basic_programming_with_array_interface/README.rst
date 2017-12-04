This tutorial shows some basics using Cassiopee array API. The Cassiopee array interface is really close to the numpy interface.

We need first to provide access to the module functions (here Converter and Generator):
import Converter as C
import Generator as G

Create a structured array corresponding to a regular Cartesian grid:
a = G.cart( (0,0,0), (1,1,1), (10,11,12) )

a is then a structured Cassiopee array: ['x,y,z', n, ni, nj, nk].
This is basically a python list containing a string describing the variables, a numpy storing the data, and 3 integers denoting the number of points in each direction.
You can easily access the numpy with:
n = a[1]
This numpy array can be then classically manipulated with numpy functions.
Convert the previous array as an HEXA unstructured array:
b = C.convertArray2Hexa(a)
It should be noted that Cassiopee array functions always return a copy of the input array.
b is then a unstructured Cassiopee array: ['x,y,z', n, c, 'HEXA'].
This is also a python list containing a string describing the variables, a numpy storing the data, a numpy storing the connectivity and a string desiging the type of elements.

[Download python script].