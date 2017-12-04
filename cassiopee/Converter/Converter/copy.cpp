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

# include <stdio.h>
# include <stdlib.h>

# include "converter.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Copy the contain of an array in another array */
//=============================================================================
PyObject* K_CONVERTER::copy(PyObject* self, PyObject* args)
{
  PyObject* array;  
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;

  res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, 
                              eltType, true);

  if (res == 1)
  { 
    PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString, 
                                        ni, nj, nk);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    K_KCORE::memcpy__(fnp, f->begin(), f->getSize()*f->getNfld());
    //memcpy(fnp, f->begin(), f->getSize()*f->getNfld()*sizeof(E_Float));
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  {
    E_Int csize = cn->getSize()*cn->getNfld();
    PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString,
                                        f->getSize(), cn->getSize(), 
                                        -1, eltType, false, csize);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    //memcpy(fnp, f->begin(), f->getSize()*f->getNfld()*sizeof(E_Float));
    K_KCORE::memcpy__(fnp, f->begin(), f->getSize()*f->getNfld());
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    return NULL;
  }
}
