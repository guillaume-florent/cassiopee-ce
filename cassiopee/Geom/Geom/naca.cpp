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
# include "geom.h"

using namespace std;
using namespace K_FLD;

extern "C"
{
  void k6naca2_(const E_Float& e, E_Int& N,
                E_Float* x, E_Float* y, E_Float* z);
}

// ============================================================================
/* Create a naca profile line of N points */
// ============================================================================
PyObject* K_GEOM::nacaMesh(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float e;
  if (!PYPARSETUPLE(args,
                    "dl", "di",
                    "fl", "fi",
                    &e, &N))
  {
      return NULL;
  }

  // Data check
  if (e < 0)
  {
    PyErr_SetString(PyExc_ValueError, "naca: thickness must be positive.");
    return NULL;
  }

  if ((N/2)*2-N == 0)
  {
    printf("Warning: naca: number of points must be odd.\n");
    printf("Warning: naca: number of points set to %d.\n", N+1);
    N = N+1;
  }
  
  // Create a naca profile
  E_Int n = E_Int(N);
  FldArrayF coord(N, 3);
  coord.setAllValuesAtNull();
  k6naca2_(e, n, coord.begin(1), coord.begin(2), coord.begin(3));
  coord.reAllocMat(n, 3);
  
  // Build array
  PyObject* tpl = K_ARRAY::buildArray(coord, "x,y,z", n, 1, 1);
  return tpl;
}
