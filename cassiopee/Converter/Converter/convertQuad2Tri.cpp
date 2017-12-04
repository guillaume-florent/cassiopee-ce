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

// convertit un maillage quad en maillage tri en coupant suivant l'angle max

# include "converter.h" 
# include "kcore.h"
# include <string.h>
# include <stdio.h>

using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Conversion du maillage quad en maillage tri.
   On coupe suivant l'angle alpha maximum. */
//=============================================================================
PyObject* K_CONVERTER::convertQuad2Tri(PyObject* self, PyObject* args)
{
  PyObject* pArray;
  if (!PYPARSETUPLEF(args, "O", "O", &pArray)) return NULL;
 
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(pArray, varString,
                                    f, ni, nj, nk, cn, eltType, true);

  // Test non structure ?
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(pArray,f);
    PyErr_SetString(PyExc_TypeError,
                    "convertQuad2tri: input array must be unstructured.");
    return NULL;
  }

  // Test QUAD
  if (strcmp(eltType, "QUAD") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertQuad2Tri: unstructured array must be QUAD.");
    RELEASESHAREDU(pArray, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertQuad2Tri: coord must be present in array.");
    RELEASESHAREDU(pArray, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  // Connectivite
  E_Int nv = f->getSize();
  E_Int ne = cn->getSize();
  vector< vector<E_Int> > cEEN(ne);
  K_CONNECT::connectEV2EENbrs(eltType, nv, *cn, cEEN);

  // TRI de sortie
  E_Int nt = 2*ne;
  FldArrayI ct(nt, 3);
  E_Int* ct1 = ct.begin(1);
  E_Int* ct2 = ct.begin(2);
  E_Int* ct3 = ct.begin(3);
  
  E_Float alpha1, alpha2;
  E_Int ind1, ind2, ind3, ind4;
  E_Float ptA1[3], ptB1[3], ptC1[3];
  E_Float ptA2[3], ptB2[3], ptC2[3];
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  FldArrayI& cnp = *cn;
  E_Int* cnp1 = cnp.begin(1);
  E_Int* cnp2 = cnp.begin(2);
  E_Int* cnp3 = cnp.begin(3);
  E_Int* cnp4 = cnp.begin(4);

  for (E_Int i = 0; i < ne; i++)
  {
    ind1 = cnp1[i]-1; ind2 = cnp2[i]-1; ind3 = cnp3[i]-1; ind4 = cnp4[i]-1;
    // test la premiere coupe (1-2-3 / 1-3-4)
    ptA1[0] = x[ind1]; ptA1[1] = y[ind1]; ptA1[2] = z[ind1];
    ptB1[0] = x[ind2]; ptB1[1] = y[ind2]; ptB1[2] = z[ind2];
    ptC1[0] = x[ind3]; ptC1[1] = y[ind3]; ptC1[2] = z[ind3];

    ptA2[0] = x[ind1]; ptA2[1] = y[ind1]; ptA2[2] = z[ind1];
    ptB2[0] = x[ind3]; ptB2[1] = y[ind3]; ptB2[2] = z[ind3];
    ptC2[0] = x[ind4]; ptC2[1] = y[ind4]; ptC2[2] = z[ind4];

    alpha1 = K_COMPGEOM::getAlphaAngleBetweenTriangles(ptA1, ptB1, ptC1, ptA2, ptB2, ptC2);

    // test la deuxieme coupe (1-2-4 / 2-3-4)
    ptA1[0] = x[ind1]; ptA1[1] = y[ind1]; ptA1[2] = z[ind1];
    ptB1[0] = x[ind2]; ptB1[1] = y[ind2]; ptB1[2] = z[ind2];
    ptC1[0] = x[ind4]; ptC1[1] = y[ind4]; ptC1[2] = z[ind4];

    ptA2[0] = x[ind1]; ptA2[1] = y[ind1]; ptA2[2] = z[ind1];
    ptB2[0] = x[ind3]; ptB2[1] = y[ind3]; ptB2[2] = z[ind3];
    ptC2[0] = x[ind4]; ptC2[1] = y[ind4]; ptC2[2] = z[ind4];

    alpha2 = K_COMPGEOM::getAlphaAngleBetweenTriangles(ptA1, ptB1, ptC1, ptA2, ptB2, ptC2);

    if (alpha1 > alpha2) // decoupage suivant 1
    {
      ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind3+1;
      ct1[2*i+1] = ind1+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;      
    }
    else // decoupage suivant 2
    {
      ct1[2*i] = ind1+1; ct2[2*i] = ind2+1; ct3[2*i] = ind4+1;
      ct1[2*i+1] = ind2+1; ct2[2*i+1] = ind3+1; ct3[2*i+1] = ind4+1;      
    }
  }

  PyObject* tpl = K_ARRAY::buildArray(*f, varString, ct, 2, NULL);
  return tpl;
}
