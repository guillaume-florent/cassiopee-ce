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

# include "transform.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Translate (zone, in place) */
// ============================================================================
PyObject* K_TRANSFORM::_translate(PyObject* self, PyObject* args)
{
  E_Float vx, vy, vz;
  PyObject* zone;
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  if (!PYPARSETUPLEF(args,
                    "O(ddd)sss", "O(fff)sss",
                    &zone, &vx, &vy, &vz, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
  {
      return NULL;
  }

  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  E_Int res = K_PYTREE::getFromZone(
    zone, 1, 0, varString, fields, locs, im, jm, km, 
    cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
    FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "translate: cannot find coordinates in zone.");
    return NULL;
  }
  // translation du vecteur (vx,vy,vz) des coordonnees
  E_Int npts;
  if (res == 1) npts = im*jm*km;
  else npts = im;
  E_Float* xt = fields[posx]; 
  E_Float* yt = fields[posy];
  E_Float* zt = fields[posz]; 
#pragma omp parallel default(shared)
  {
#pragma omp for 
    for (E_Int ind = 0; ind < npts; ind++)
    {
      xt[ind] += vx; yt[ind] += vy; zt[ind] += vz;
    }
  }
  delete [] eltType;
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  Py_INCREF(Py_None); 
  return Py_None;
}

// ============================================================================
// Translate array2 - in place
//=============================================================================
PyObject* K_TRANSFORM::_translate2(PyObject* self, PyObject* args)
{
  E_Float vx, vy, vz;
  PyObject* array;
  if (!PYPARSETUPLEF(args,
                    "O(ddd)", "O(fff)",
                    &array, &vx, &vy, &vz))
      return NULL;
  
  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                                  cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "translate: invalid array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
   
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "translate: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
    
  E_Int npts = f->getSize();

  E_Float* xt = f->begin(posx); 
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

 #pragma omp parallel default(shared)
  {
#pragma omp for 
    for (E_Int ind = 0; ind < npts; ind++)
    {
      xt[ind] += vx; yt[ind] += vy; zt[ind] += vz;
    }
  }
  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}
