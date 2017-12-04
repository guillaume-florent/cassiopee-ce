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
#include <string.h>
#include "Array/Array.h"
#include "String/kstring.h"

using namespace K_FLD;

//=============================================================================
// Extrait les donnees (shared) d'un objet python struct array
// defini par: Array1: [ 'vars', a, ni, nj, nk ]
//             Array2: [ 'vars', [a], ni, nj, nk ]
// ou d'un objet python unstruct array
// defini par: Array1: [ 'vars', a, c, "ELTTYPE"]
//             Array2: [ 'vars', [a], [c], "ELTTYPE"]
// ou ELTTYPE vaut: NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON
// avec ou sans star.
// return 1: valid struct array
//           f: field en stockage Fld (champs par variable)
//           ni, nj, nk: number of points
//           varString
// return 2: valid unstruct array
//           f: field en stockage Fld (champs par variable)
//           c: connectivity (champs des indices de noeuds commencants a 1)
//           eltType: type of elements
//           varString
// Dans tous les cas sauf NGON: c est la connectivite elt-noeuds 
// dimensionnee (elt, nbre de noeuds par elements)
// Dans le cas NGON: c est la connectivite Face-noeuds, suivi de la 
// connectivite elt-faces + 4 entiers indiquant les tailles
//
// return -1: given object is not a list.
// return -2: not a valid number of elts in list.
// return -3: first element is not a var string.
// return -4: a is not a valid numpy array.
// return -5: array is structured but ni, nj, nk unvalid.
// return -6: array is unstructured but connectivity is unvalid.
// return -7: array is unstructured but elt type is unknown.
// C'est la responsabilite de l'appelant de liberer la memoire de f et 
// eventuellement de c
//=============================================================================
E_Int K_ARRAY::getFromArray2(PyObject* o,
                             char*& varString,
                             FldArrayF*& f,
                             E_Int& ni, E_Int& nj, E_Int& nk,
                             FldArrayI*& c,
                             char*& eltType)
{
  PyObject* tpl; PyObject* p;
  PyArrayObject* a; PyArrayObject* ac;
  PyObject* ref; //PyObject* ref2;
  IMPORTNUMPY;

  // -- list --
  if (PyList_Check(o) == false)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
    return -1;
  }
  E_Int size = PyList_Size(o);
  if (size != 4 && size != 5)
  {
    PyErr_Warn(PyExc_Warning, 
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check number of elements in list.");
    return -2;
  }
 
  // -- varString --
  if (PyString_Check(PyList_GetItem(o,0)) == false)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. First element must be a string.");
    return -3;
  }
  // pointeur sur la chaine python
  varString = PyString_AsString(PyList_GetItem(o,0));
  E_Int nvar = getNumberOfVariables(varString);

  // -- field --
  tpl = PyList_GetItem(o, 1);
  //ref = tpl; Py_INCREF(ref); // trick
  Py_INCREF(tpl); ref = tpl;

  if (PyArray_Check(tpl) == true) // -- Array1 --
  {
    a = (PyArrayObject*)tpl;

    if (PyArray_NDIM(a) != 2)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: field must have two dimensions.");
      Py_DECREF(tpl);
      return -4;
    }
  
    E_Int s = PyArray_DIMS(a)[1];
    E_Int nfld = PyArray_DIMS(a)[0];

    if (nfld != nvar)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: number of variables different in varString and field.");
      Py_DECREF(tpl);
      return -4;
    }    
    f = new FldArrayF(s, nfld, (E_Float*)PyArray_DATA(a), true, true);
  }
  else if (PyList_Check(tpl) == true) // -- Array2 --
  {
    E_Int nfld = PyList_Size(tpl);
    E_Float** acu = new E_Float* [nfld];
    E_Int s = 0;
    for (E_Int i = 0; i < nfld; i++)
    {
      p = PyList_GetItem(tpl, i);
      if (PyArray_Check(p) == true)
      {
        a = (PyArrayObject*)p; //Py_INCREF(a);
        s = PyArray_SIZE(a); // must be constant
        acu[i] = (E_Float*)PyArray_DATA(a);
      }
      else
      {
        PyErr_Warn(PyExc_Warning,
                   "getFromArray: field must a list of numpys.");
        Py_DECREF(tpl);
        return -4;
      }
    }
    f = new FldArrayF(s, nfld, acu, true, true);
    delete [] acu;
  }
  else
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: second arg in array must be a numpy array.");
    return -4;
  }
  
  if (size == 4) // unstruct array
  {
    // -- connectivity --
    tpl = PyList_GetItem(o, 2);
    //ref2 = tpl; Py_INCREF(ref2);
    Py_INCREF(tpl);
    if (PyArray_Check(tpl) == true) // -- Array1 --
    {
      ac = (PyArrayObject*)tpl;
  
      if (PyArray_NDIM(ac) > 2)
      {
        PyErr_Warn(PyExc_Warning,
                   "getFromArray: connectivity must have two dimensions max.");
        Py_DECREF(ref); Py_DECREF(tpl);
        return -4;
      }
      E_Int s, nfld;
      if (PyArray_NDIM(ac) == 2)
      { s = PyArray_DIMS(ac)[1]; nfld = PyArray_DIMS(ac)[0]; }
      else { s = PyArray_DIMS(ac)[0]; nfld = 1; }
      c = new FldArrayI(s, nfld, (E_Int*)PyArray_DATA(ac), true);
    }
    else if (PyList_Check(tpl) == true) // -- Array2 --
    {
      E_Int nc = PyList_Size(tpl);
      if (nc == 1) // BE => compact + stride
      {
        ac = (PyArrayObject*)PyList_GetItem(tpl,0);
        E_Int s, nfld;
        if (PyArray_NDIM(ac) == 2)
        { s = PyArray_DIMS(ac)[0]; nfld = PyArray_DIMS(ac)[1]; }
        else { s = PyArray_DIMS(ac)[0]; nfld = 1; }
        c = new FldArrayI(s, nfld, (E_Int*)PyArray_DATA(ac), true, false);
      }
      else // suppose NGON
      {
        if (nc == 4) // NGON/NFACE
        {
          PyArrayObject* ac1 = (PyArrayObject*)PyList_GetItem(tpl,0); // NGON
          PyArrayObject* ac2 = (PyArrayObject*)PyList_GetItem(tpl,1); // NFACE
          PyArrayObject* ac3 = (PyArrayObject*)PyList_GetItem(tpl,2); // indPG
          PyArrayObject* ac4 = (PyArrayObject*)PyList_GetItem(tpl,3); // indPH
          E_Int nfaces = PyArray_SIZE(ac3);
          E_Int nelts = PyArray_SIZE(ac4);
          E_Int sizeNGon = PyArray_SIZE(ac1);
          E_Int sizeNFace = PyArray_SIZE(ac2);
          c = new FldArrayI(nfaces, nelts, (E_Int*)PyArray_DATA(ac1), (E_Int*)PyArray_DATA(ac2),
                            (E_Int*)PyArray_DATA(ac3), (E_Int*)PyArray_DATA(ac4), sizeNGon, sizeNFace);
        }
      }
    }
    else
    {
      PyErr_Warn(PyExc_Warning, 
                 "getFromArray: third arg in array must be a numpy array.");
      Py_DECREF(ref); Py_DECREF(tpl);
      return -6;
    }

    // -- element type --
    if (PyString_Check(PyList_GetItem(o,3)) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: an unstruct array must be of list of type ['vars', a, c, 'ELTTYPE']. Last element must be a string.");
      Py_DECREF(ref); Py_DECREF(tpl);
      return -7;
    }
    eltType = PyString_AsString(PyList_GetItem(o,3));

    if (K_STRING::cmp(eltType, "NODE") != 0 &&
        K_STRING::cmp(eltType, "BAR") != 0 &&
        K_STRING::cmp(eltType, "TRI") != 0 &&
        K_STRING::cmp(eltType, "QUAD") != 0 &&
        K_STRING::cmp(eltType, "TETRA") != 0 &&
        K_STRING::cmp(eltType, "PYRA") != 0 &&
        K_STRING::cmp(eltType, "PENTA") != 0 &&
        K_STRING::cmp(eltType, "HEXA") != 0 &&
        K_STRING::cmp(eltType, "NGON") != 0 &&
        K_STRING::cmp(eltType, "NODE*") != 0 &&
        K_STRING::cmp(eltType, "BAR*") != 0 &&
        K_STRING::cmp(eltType, "TRI*") !=0 &&
        K_STRING::cmp(eltType, "QUAD*") != 0 &&
        K_STRING::cmp(eltType, "TETRA*") !=0 &&
        K_STRING::cmp(eltType, "PYRA*") != 0 &&
        K_STRING::cmp(eltType, "PENTA*") != 0 &&
        K_STRING::cmp(eltType, "HEXA*") != 0 &&
        K_STRING::cmp(eltType, "NGON*") != 0)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: element type unknown: %s. Must be in NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON or NODE*, BAR*, TRI*, QUAD*, TETRA*, PYRA*, PENTA*, HEXA*, NGON*.");
      Py_DECREF(ref); Py_DECREF(tpl);
      return -7;
    }
    return 2;
  }
  else // struct array
  {
    tpl = PyList_GetItem(o,2);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: third arg must be an integer.");
      Py_DECREF(ref);
      return -5;
    }
    ni = PyLong_AsLong(tpl);
    tpl = PyList_GetItem(o,3);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: fourth arg must be an integer.");
      Py_DECREF(ref);
      return -5;
    }
    nj = PyLong_AsLong(tpl);
    tpl = PyList_GetItem(o,4);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: fifth arg must be an integer.");
      Py_DECREF(ref);
      return -5;
    }
    nk = PyLong_AsLong(tpl);
    return 1;
  }
}
