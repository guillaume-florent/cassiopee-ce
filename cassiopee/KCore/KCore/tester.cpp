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

#include "kcore.h"
#include <stdlib.h>

//==============================================================================
PyObject* K_KCORE::tester(PyObject* self, PyObject* args)
{
#define TESTARRAYN

#ifdef TEST1
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  E_Float* t = new E_Float[1]; t[0] = 12.;
  PyObject* p = K_PYTREE::createChild(o, "coucou", "DataArray_t", 
                                      t, 1, 1);
  delete [] t;
  return p;
#endif

#ifdef TEST2
  PyObject* o; char* path;
  if (!PyArg_ParseTuple(args, "Os", &o, &path)) return NULL;
  PyObject* p = K_PYTREE::getNodeFromPath(o, path);
  return p;
#endif

#ifdef TEST3
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  f->print();
  E_Float* x = f->begin();
  x[0] = -0.05;
  if (ret == 2) c->print();
  RELEASESHAREDB(ret, o, f, c);
  Py_INCREF(Py_None);
  return Py_None;
#endif

#ifdef TESTARRAY2
  // Structured array1 - build - ni=5,nj=5,nk=5
  PyObject* o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 5,5,5, 1);
  // Structured array1 - get
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType; E_Int ni,nj,nk;
  E_Int ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // Structured array1 - usage
  f->setAllValuesAt(1);
  // Getting information from Fld
  E_Int nfld = f->getNfld(); // nbre de champs
  E_Int size = f->getSize(); // nbre de pts (ni x nj x nk)
  E_Float* x = f->begin(1); // ptrs sur les champs
  E_Float* y = f->begin(2);

  // Structured array1 - free
  RELEASESHAREDB(ret, o, f, c);

  // Structured array2 - build - ni=2,nj=2,nk=2
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 2,2,2, 2);
  // Classic array2 - get
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // Structured array2 - usage
  f->setAllValuesAt(2.2);
  // Getting information from Fld f
  nfld = f->getNfld(); // nbre de champs
  size = f->getSize(); // nbre de pts
  x = f->begin(1);
  y = f->begin(2);

  // Structured array2 - free
  RELEASESHAREDB(ret, o, f, c);

  // Hexa array1 - build - nvertex=10, nelts=5
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 10, 5,-1, "HEXA", false, 0, 0, 0, 1);
  // Hexa array1 - get
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // Hexa array1 - usage
  f->setAllValuesAt(1);
  c->setAllValuesAt(1);
  // Getting information from Fld c
  E_Int ne = c->getSize(); // nbre d'elements
  E_Int nv = c->getNfld(); // nbre de vertex par element
  // Hexa array1 - free
  RELEASESHAREDB(ret, o, f, c);

  // Hexa array2 - build - nvertex=10, nelts=5
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 10, 5,-1, "HEXA", false, 0, 0, 0, 2);
  // Hexa array2 - get
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // Hexa array1 - usage
  f->setAllValuesAt(1);
  c->setAllValuesAt(1);
  // Getting information from Fld c
  ne = c->getSize(); // nbre d'elements
  nv = c->getNfld(); // nbre de vertex par element
  // Hexa array2 - free
  RELEASESHAREDB(ret, o, f, c);

  // NGon array1 - build - nvertex=10, nelts=5, nface=6, sizeNGon=20, sizeNFace=30
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 10, 5,-1, "NGON", false, 20, 30, 6, 1);
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // NGon array1 - usage
  f->setAllValuesAt(1);
  c->setAllValuesAt(1);
  E_Int nelts = c->getNElts();
  E_Int nfaces = c->getNFaces();
  E_Int* ngon = c->getNGon();
  E_Int* nface = c->getNFace();
  E_Int* indPG = c->getIndPG();
  E_Int* indPH = c->getIndPH();
  RELEASESHAREDB(ret, o, f, c);
  //printf("nelts=%d, nfaces=%d\n", nelts, nfaces);

  // NGon array2 - build
  o = K_ARRAY::buildArray2(5, "x,y,z,F,G", 10, 5,-1, "NGON", false, 20, 30, 6, 2);
  ret = K_ARRAY::getFromArray2(o, varString, f, ni, nj, nk, c, eltType);
  // NGon array1 - usage
  f->setAllValuesAt(1);
  nelts = c->getNElts();
  nfaces = c->getNFaces();
  ngon = c->getNGon();
  nface = c->getNFace();
  indPG = c->getIndPG();
  indPH = c->getIndPH();
  //printf("nelts=%d, nfaces=%d\n", nelts, nfaces);
  RELEASESHAREDB(ret, o, f, c);
  return o;
#endif

#ifdef TESTFLD
  // test nouveau FldArray
  K_FLD::FldArrayF t1; 

  //=========================
  // compact, fortranOrdered
  //========================
  K_FLD::FldArrayF t2(120, 3, true);
  t2 = 0.;
  t2.setAllValuesAt(1);
  t2.setAllValuesAtNull();
  for (E_Int i = 0; i < 120; i++)
  {
    t2(i,1) = i;
    t2(i,2) = 0.1*i;
    t2(i,3) = 0.01*i;
  }
  K_FLD::FldArrayF t2p(t2);
  printf("indmin %d %d %d\n", t2p.indMin(1,1), t2p.indMin(1,2), t2p.indMin(1,3));
  printf("indmax %d %d %d\n", t2p.indMax(1,1), t2p.indMax(1,2), t2p.indMax(1,3));
  printf("1: %f %f %f\n", t2(10,1), t2(10,2), t2(10,3));
  printf("1: %f %f %f\n", t2[10], t2[10+120], t2[10+240]);
  E_Float* pt1 = t2.begin(1);
  E_Float* pt2 = t2.begin(2);
  E_Float* pt3 = t2.begin(3);
  printf("1: %f %f %f\n", pt1[10], pt2[10], pt3[10]);
  
  K_FLD::FldArrayF t2t(120, 3, t2.begin(), true); // shared with t2
  t2t(10,2) = 12.;

  // Build shared
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z", 10, 10, 10);
  E_Float* sp = K_ARRAY::getFieldPtr(tpl);
  K_FLD::FldArrayF s(10*10*10, 3, sp, true);
  
  for (E_Int n = 1; n <= 3; n++)
  {
    E_Float* spi = s.begin(n);
    for (E_Int i = 0; i < 10*10*10; i++) spi[i] = 12.;
  }

  //============================
  // non compact, fortranOrdered
  //===========================
  K_FLD::FldArrayF t3(120, 3, false);
  t3 = 0.;
  t3.setAllValuesAt(1);
  t3.setAllValuesAtNull();
  for (E_Int i = 0; i < 120; i++)
  {
    t3(i,1) = i;
    t3(i,2) = 0.1*i;
    t3(i,3) = 0.01*i;
  }
  K_FLD::FldArrayF t3p(t3);
  printf("indmin %d %d %d\n", t3p.indMin(1,1), t3p.indMin(1,2), t3p.indMin(1,3));
  printf("indmax %d %d %d\n", t3p.indMax(1,1), t3p.indMax(1,2), t3p.indMax(1,3));
  printf("2: %f %f %f\n", t3(10,1), t3(10,2), t3(10,3));
  printf("2: %f %f %f\n", t3[10], t3[10+120], t3[10+240]);
  pt1 = t3.begin(1);
  pt2 = t3.begin(2);
  pt3 = t3.begin(3);
  printf("2: %f %f %f\n", pt1[10], pt2[10], pt3[10]);

  //===================
  // compact, C ordered
  //===================
  K_FLD::FldArrayF t4(120, 3, true, false);
  t4 = 0.;
  t4.setAllValuesAt(1);
  t4.setAllValuesAtNull();
  for (E_Int i = 0; i < 120; i++)
  {
    t4(i,1) = i;
    t4(i,2) = 0.1*i;
    t4(i,3) = 0.01*i;
  }
  K_FLD::FldArrayF t4p(t4);
  printf("indmin %d %d %d\n", t4p.indMin(1,1), t4p.indMin(1,2), t4p.indMin(1,3));
  printf("indmax %d %d %d\n", t4p.indMax(1,1), t4p.indMax(1,2), t4p.indMax(1,3));
  printf("3: %f %f %f\n", t4(10,1), t4(10,2), t4(10,3));
  printf("3: %f %f %f\n", t4[10], t4[10+120], t4[10+240]);
  pt1 = t4.begin(1);
  pt2 = t4.begin(2);
  pt3 = t4.begin(3);
  printf("3: %f %f %f\n", pt1[10*3], pt2[10*3], pt3[10*3]);
#endif

  Py_INCREF(Py_None);
  return Py_None;
}
