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

// normalize a list of vars

# include "converter.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
// normalize a list of vars
//=============================================================================
PyObject* K_CONVERTER::normalize(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* varsO;
  if (!PyArg_ParseTuple(args, "OO", &array, &varsO )) return NULL;

  // Check vars
  if (PyList_Check(varsO) == 0)
  {
    PyErr_SetString(PyExc_TypeError, "normalize: vars must be a list.");
    return NULL;
  }
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, eltType, true);
  
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "normalize: array is invalid.");
    return NULL;
  }

  // Extraction des variables 
  vector<E_Int> pos;
  E_Int n = 0;
  char* var;
  E_Int m;
  for (E_Int v  = 0 ; v < PyList_Size(varsO); v++)
  {
    if (PyString_Check(PyList_GetItem(varsO, v)) == 0)
    {
      printf("Warning: normalize: invalid string for variable %d. Skipped...\n", v);
    }
    else
    {
      var = PyString_AsString(PyList_GetItem(varsO, v));
      m = K_ARRAY::isNamePresent(var, varString);
      if (m == -1)
        printf("Warning: normalize: variable %d not present in array. Skipped...\n", v);
      else {m++; pos.push_back(m);}
    }
  }

  n = pos.size();
  if (n == 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    printf("Warning: normalize: no variable in result.\n");;
    Py_INCREF(Py_None); return Py_None;
  }
  E_Int npts = f->getSize(); E_Int nfld = f->getNfld();
  PyObject* tpl;

  if (res == 1)
  {
    tpl = K_ARRAY::buildArray(nfld, varString, 
                              ni, nj, nk);
  }
  else //unstr
  {
    tpl = K_ARRAY::buildArray(nfld, varString,
                              npts, cn->getSize(),
                              -1, eltType, false, cn->getSize()*cn->getNfld());
  }
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fp(npts, nfld, fnp, true); fp = *f;
  if (res == 2)
  {
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
  }

  E_Float* fpt0 = fp.begin(pos[0]);

#pragma omp parallel default(shared)
  {
#pragma omp for nowait
    for (E_Int i = 0; i < npts; i++) fpt0[i] = 0.;

    for (E_Int v = 0; v < n; v++)
    {
      E_Float* ft = f->begin(pos[v]);
#pragma omp for nowait
      for (E_Int i = 0; i < npts; i++) fpt0[i] += ft[i]*ft[i];
    }
#pragma omp for nowait
    for (E_Int i = 0; i < npts; i++)
      fpt0[i] = 1./K_FUNC::E_max(sqrt(fpt0[i]), 1.e-12);

    for (E_Int v  = n-1; v >= 0; v--)
    {
      E_Float* ft = f->begin(pos[v]);
      E_Float* fpt = fp.begin(pos[v]);
#pragma omp for
      for (E_Int i = 0; i < npts; i++)
      {
        fpt[i] = ft[i] * fpt0[i];
      }
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  return tpl;
} 
    
