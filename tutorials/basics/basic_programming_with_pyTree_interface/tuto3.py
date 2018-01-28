#!/usr/bin/env python
# coding: utf-8

r"""Tutorial 3 for Cassiopee

This tutorial shows some basics using Cassiopee pyTree API.
The pyTree is a mapping of CGNS standard.

We need first to provide access to the module functions
(here Converter and Generator):
import Converter.PyTree as C
import Generator.PyTree as G

Create a structured array corresponding to a regular Cartesian grid:
a = G.cart( (0,0,0), (1,1,1), (10,11,12) )

a is then a zone node. Each node of a pyTree is a python list of
of type: ['Name', n, [], 'NodeType_t'] where
  'Name' is the name if the node,
  n is a numpy containing the data,
  [] is a list of sons of this node,
  'NodeType_t' describes the type of node.
You can of course manipulate the nodes with direct access:
print a[2][0]
But, it is easier to use the Internal module.

Print to screen the zone node:
Internal.printTree(a)

You can always access the numpy storing data, for example 'CoordinateX' by:
n = Internal.getNodeFromName(a, 'CoordinateX')[1]

Convert the previous array as an HEXA unstructured array:
b = C.convertArray2Hexa(a)
Cassiopee functions return a copy of zone a. It is possible to modify directly
a without copy (the so called in-place treatment).
Simply call the function preceded by a _:
C._convertArray2Hexa(a)
a is then modified.

"""

import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

# Create a cartesian grid
a = G.cart((0, 0, 0), (1, 1, 1), (10, 11, 12))

# A node is a list of type: ['Name', n, [], 'NodeType_t']
# 'Name' is the name if the node, n is a numpy containing the data,
# [] is a list of sons of this node, 'NodeType_t' describes the type of node.

# This is a zone node
Internal.printTree(a)

# To manipulate nodes, we use the Internal module
node = Internal.getNodeFromName(a, 'CoordinateX')

# n is a numpy you can manipulate classically
n = node[1]

# Convert a zone in unstructured HEXA
# This function returns a copy of the zone
b = C.convertArray2Hexa(a)

# You can directly modify a without a copy (in place functions), using
# the function prefixed with a _
C._convertArray2Hexa(a)

# You can save as a pyTree
C.convertPyTree2File([a], 'out.cgns')
