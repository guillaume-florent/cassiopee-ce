/*    
    Copyright 2013-2017 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#define K_ARRAY_UNIQUE_SYMBOL
#include "intersector.h"

int __activation__;

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyintersector [] =
{
  
  {"conformUnstr", K_INTERSECTOR::conformUnstr, METH_VARARGS},
  
  {"booleanIntersection", K_INTERSECTOR::booleanIntersection, METH_VARARGS},
  {"booleanUnion", K_INTERSECTOR::booleanUnion, METH_VARARGS},
  {"booleanMinus", K_INTERSECTOR::booleanMinus, METH_VARARGS},
  {"booleanIntersectionBorder", K_INTERSECTOR::booleanIntersectionBorder, METH_VARARGS},
  {"booleanModifiedSolid", K_INTERSECTOR::booleanModifiedSolid, METH_VARARGS},
  {"XcellN", K_INTERSECTOR::XcellN, METH_VARARGS},
  {"P1ConservativeChimeraCoeffs", K_INTERSECTOR::P1ConservativeChimeraCoeffs, METH_VARARGS},
  {"selfX", K_INTERSECTOR::selfX, METH_VARARGS},
  {"triangulateExteriorFaces", K_INTERSECTOR::triangulateExteriorFaces, METH_VARARGS},
  {"convexifyFaces", K_INTERSECTOR::convexifyFaces, METH_VARARGS},
  {"prepareCellsSplit", K_INTERSECTOR::prepareCellsSplit, METH_VARARGS},
  {"simplifyCells", K_INTERSECTOR::simplifyCells, METH_VARARGS},
  {"splitNonStarCells", K_INTERSECTOR::splitNonStarCells, METH_VARARGS},
  {"collapseUncomputableFaces", K_INTERSECTOR::collapseUncomputableFaces, METH_VARARGS},
  {"agglomerateSmallCells", K_INTERSECTOR::agglomerateSmallCells, METH_VARARGS},
  //{"agglomerateUncomputableCells", K_INTERSECTOR::agglomerateUncomputableCells, METH_VARARGS},
  {"extractUncomputables", K_INTERSECTOR::extractUncomputables, METH_VARARGS},
  {"extractPathologicalCells", K_INTERSECTOR::extractPathologicalCells, METH_VARARGS},
  {"extractOuterLayers", K_INTERSECTOR::extractOuterLayers, METH_VARARGS},
  {"extractNthCell", K_INTERSECTOR::extractNthCell, METH_VARARGS},
  {"extractNthFace", K_INTERSECTOR::extractNthFace, METH_VARARGS},
  {"removeNthCell", K_INTERSECTOR::removeNthCell, METH_VARARGS},

  {"diffMesh", K_INTERSECTOR::diffMesh, METH_VARARGS},

  { "checkCellsClosure", K_INTERSECTOR::checkCellsClosure, METH_VARARGS },

  { "extrudeUserDefinedBC", K_INTERSECTOR::extrudeUserDefinedBC, METH_VARARGS },

  { "reorientExternalFaces", K_INTERSECTOR::reorientExternalFaces, METH_VARARGS },

  {NULL, NULL}
};

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
  void initintersector();
  void initintersector()
  {
    __activation__ = K_KCORE::activation();
    Py_InitModule("intersector", Pyintersector);
    import_array();
  }
}
