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
#include "post.h"

using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Dictionary of all functions of the python module */
// ============================================================================
static PyMethodDef Pypost [] =
{
  {"extractPoint", K_POST::extractPoint, METH_VARARGS},
  {"extractPlane", K_POST::extractPlane, METH_VARARGS},
  {"projectCloudSolution2Triangle", K_POST::projectCloudSolution2Triangle, METH_VARARGS},
  {"extractMesh", K_POST::extractMesh, METH_VARARGS},
  {"coarsen", K_POST::coarsen, METH_VARARGS},
  {"refine", K_POST::refine, METH_VARARGS},
  {"refineButterfly", K_POST::refineButterfly, METH_VARARGS},
  {"selectCells", K_POST::selectCells, METH_VARARGS},
  {"selectCells3", K_POST::selectCells3, METH_VARARGS},
  {"selectCellCenters", K_POST::selectCellCenters, METH_VARARGS},
  {"interiorFaces", K_POST::selectInteriorFaces, METH_VARARGS},
  {"exteriorFaces", K_POST::selectExteriorFaces, METH_VARARGS},
  {"exteriorFacesStructured", K_POST::selectExteriorFacesStructured, METH_VARARGS},
  {"exteriorElts", K_POST::selectExteriorElts, METH_VARARGS},
  {"frontFaces", K_POST::frontFaces, METH_VARARGS},
  {"integ", K_POST::integ, METH_VARARGS},
  {"integ2", K_POST::integ2, METH_VARARGS},
  {"integNorm", K_POST::integNorm, METH_VARARGS},
  {"integNormProduct", K_POST::integNormProduct, METH_VARARGS},
  {"integMoment", K_POST::integMoment, METH_VARARGS},
  {"integMomentNorm", K_POST::integMomentNorm, METH_VARARGS},
  {"zipper", K_POST::zipperF, METH_VARARGS},
  {"usurp", K_POST::usurpF, METH_VARARGS},
  {"computeVariables", K_POST::computeVariables, METH_VARARGS},
  {"computeGrad", K_POST::computeGrad, METH_VARARGS},
  {"computeGrad2NGon", K_POST::computeGrad2NGon, METH_VARARGS},
  {"computeGrad2Struct", K_POST::computeGrad2Struct, METH_VARARGS},
  {"computeNormGrad", K_POST::computeNormGrad, METH_VARARGS},
  {"computeCurl", K_POST::computeCurl, METH_VARARGS},
  {"computeNormCurl", K_POST::computeNormCurl, METH_VARARGS},
  {"computeDiff", K_POST::computeDiff, METH_VARARGS},
  {"perlinNoise", K_POST::perlinNoise, METH_VARARGS},
  {"compStreamLine", K_POST::compStreamLine, METH_VARARGS},
  {"compStreamRibbon", K_POST::compStreamRibbon, METH_VARARGS},
  {"compStreamSurf", K_POST::compStreamSurf, METH_VARARGS},
  {"isoLine", K_POST::isoLine, METH_VARARGS},
  {"isoSurf", K_POST::isoSurf, METH_VARARGS},
  {"isoSurfMC", K_POST::isoSurfMC, METH_VARARGS},
  {"isoSurfMC_opt", K_POST::isoSurfMC_opt, METH_VARARGS},
  {"enforceIndicatorNearBodies", K_POST::enforceIndicatorNearBodies, METH_VARARGS},
  {"enforceIndicatorForFinestLevel", K_POST::enforceIndicatorForFinestLevel, METH_VARARGS},
  {"enforceIndicatorForCoarsestLevel", K_POST::enforceIndicatorForCoarsestLevel, METH_VARARGS},
  {"sharpEdges", K_POST::sharpEdges, METH_VARARGS},
  {"silhouette", K_POST::silhouette, METH_VARARGS},
  {"computeIndicatorValue",K_POST::computeIndicatorValue,METH_VARARGS},
  {NULL, NULL}
};

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
  void initpost();
  void initpost()
  {
    Py_InitModule("post", Pypost);
    import_array();
  }
}
