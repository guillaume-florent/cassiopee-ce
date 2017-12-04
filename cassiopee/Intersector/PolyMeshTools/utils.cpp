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


# include <string>
# include <sstream> 
# include "intersector.h"
# include "Fld/ngon_t.hxx"
# include "Nuga/Delaunay/Triangulator.h"

//#include <iostream>
#include <memory>

using namespace std;
using namespace K_FLD;


//=============================================================================
/* Creates 4 zones : 1) uncomputable polygons 2) uncomputable polyhedra 3) uncomputable polyhedra & neighbors 4) complementary of 3)*/
//=============================================================================
PyObject* K_INTERSECTOR::extractUncomputables(PyObject* self, PyObject* args)
{
  E_Int neigh_level(1);
  PyObject *arr;

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &neigh_level)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), uphs/*uncomputable phs*/, uphs_wv1/*uncomputable phs with neighbors*/, remaining;
  ngon_unit upgs; //uncomputable polygons

  err = ngon_type::extract_uncomputables<DELAUNAY::Triangulator>(crd, ngi, neigh_level, upgs, uphs, uphs_wv1, remaining);
  
  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "extract_uncomputables failed.");
    delete f; delete cn;
    return NULL;
  } 
  
  PyObject *l(PyList_New(0)), *tpl;

  if (upgs.size() == 0)
  {
    std::cout << "OK : there are no uncomputable polygons." << std::endl;
    PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnt, -1, eltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  else
  {
    // zone 1 : uncomputable pgs
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(upgs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      ngon_type ngo(upgs, false/*one ph per pg*/);
      ngo.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    // zone 2 : uncomputable phs
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(uphs.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      uphs.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    // zone 3 : uncomputable phs and first neighborhood
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(uphs_wv1.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      uphs_wv1.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    // zone 4 : complementary of above selection
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(remaining.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      remaining.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
  }
  
  
  delete f; delete cn;
  return l;
}

//=============================================================================
/* XXX*/
//=============================================================================
PyObject* K_INTERSECTOR::extractPathologicalCells(PyObject* self, PyObject* args)
{
  E_Int neigh_level(2);
  PyObject *arr;
  E_Float concave_threshold(E_EPSILON);

  if (!PYPARSETUPLE(args, 
                    "Odl", "Odi", "Ofl", "Ofi", 
                    &arr, &concave_threshold, &neigh_level)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), neigh_phs/*patho neighbors*/, remaining_phs;
  std::vector<ngon_type> phsv;

  //std::cout << "neigh level : " <<  neigh_level << std::endl;
  std::cout << "concave_threshold : " <<  concave_threshold << std::endl;
 
  err = ngon_type::extract_pathological_PHs<DELAUNAY::Triangulator>(crd, ngi, concave_threshold, neigh_level, 
                                                                    phsv, neigh_phs, remaining_phs);
  
  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "extract_pathological_phs failed.");
    delete f; delete cn;
    return NULL;
  } 
  
  PyObject *l(PyList_New(0)), *tpl;

  if (phsv.empty())
  {
    std::cout << "OK : all the cells are star-shaped regarding there centroids." << std::endl;
    tpl = K_ARRAY::buildArray(crd, varString, cnt, -1, eltType, false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  else
  {
    E_Int nphs = phsv[0].PHs.size();
    //if (nphs)
    // zone 1 : 
    {
      if (nphs) std::cout << "there are " << nphs << " open cells (bug somewhere)." << std::endl;

      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(phsv[0].PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      phsv[0].export_to_array(cnto);
      //patho_name = "open_cells";
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    nphs = phsv[1].PHs.size();
    //if (nphs)
    // zone 2 : 
    {
      if (nphs) std::cout << "there are " << nphs << " non handled cells (showstopper?)." << std::endl;

      K_FLD::FloatArray crdtmp(crd);

      ngon_type::compact_to_used_nodes(phsv[1].PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      phsv[1].export_to_array(cnto);

      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }

    nphs = phsv[2].PHs.size();
    //if (nphs)
    // zone 3 : 
    {
      if (nphs) std::cout << "there are " << phsv[2].PHs.size() << " cells with some uncomputable polygons." << std::endl;
      K_FLD::FloatArray crdtmp(crd);

      ngon_type::compact_to_used_nodes(phsv[2].PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      phsv[2].export_to_array(cnto);

      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    nphs = phsv[3].PHs.size();
    if (nphs) std::cout << "there are " << phsv[3].PHs.size() << " non-centroid-star-shaped cells that need to be split." << std::endl;
    // zone 4 : 
    {
      K_FLD::FloatArray crdtmp(crd);
      //ngon_type::compact_to_used_nodes(phsv[3].PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      phsv[3].export_to_array(cnto);

      //std::cout << "pg : pg " << phsv[3].PGs.size() << " versus " << cnto[0] << std::endl;

      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }

    // zone 5 : neighbors
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(neigh_phs.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      neigh_phs.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }

    // zone 6 : complementary of above selections
    if (remaining_phs.PHs.size())
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(remaining_phs.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      remaining_phs.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
  }
  
  
  delete f; delete cn;
  return l;
}


//=============================================================================
/* Creates 2 zones : 1) outerlayer with firt neighborhoo 2) complementary */
//=============================================================================
PyObject* K_INTERSECTOR::extractOuterLayers(PyObject* self, PyObject* args)
{

  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), outer, remaining;
 
  err = ngon_type::extract_outer_layers(ngi, outer, remaining);
  
  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "extract_outer_layers failed.");
    delete f; delete cn;
    return NULL;
  } 
  
  PyObject *l(PyList_New(0)), *tpl;

  {
    //std::cout << "there are " << outer.PHs.size() << " outer cells detected." << std::endl;
    // zone 1 : outer
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(outer.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      outer.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
    
    // zone 2 : complementary 
    {
      K_FLD::FloatArray crdtmp(crd);
      ngon_type::compact_to_used_nodes(remaining.PGs, crdtmp); //reduce points
      //export to numpy
      K_FLD::IntArray cnto;
      remaining.export_to_array(cnto);
      tpl = K_ARRAY::buildArray(crdtmp, varString, cnto, -1, eltType, false);
      PyList_Append(l, tpl);
      Py_DECREF(tpl);
    }
  }
  
  
  delete f; delete cn;
  return l;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extractNthCell(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int nth(0);

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &nth)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);
 
  
  ngon_unit ph;
  ph.add(ngi.PHs.stride(nth), ngi.PHs.get_facets_ptr(nth));
  
  ngon_type one_ph(ngi.PGs, ph);

  std::vector<E_Int> pgnids, phnids;
  one_ph.remove_unreferenced_pgs(pgnids, phnids);
  ngon_type::compact_to_used_nodes(one_ph.PGs, crd);
  
  K_FLD::IntArray cnto;
  one_ph.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::removeNthCell(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int nth(0);

  if (!PyArg_ParseTuple(args, "Ol", &arr, &nth)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);
  
  ngon_unit phs;
  for (E_Int i = 0; i < ngi.PHs.size(); ++i)
  {
    if (i == nth) continue;
    phs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
  }
  
  ngon_type ng(ngi.PGs, phs);

  std::vector<E_Int> pgnids, phnids;
  ng.remove_unreferenced_pgs(pgnids, phnids);
  ngon_type::compact_to_used_nodes(ng.PGs, crd);
  
  K_FLD::IntArray cnto;
  ng.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extractNthFace(PyObject* self, PyObject* args)
{

  PyObject *arr;
  E_Int nth(0);

  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arr, &nth)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;
    
  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);
 
  
  ngon_unit pg;
  pg.add(ngi.PGs.stride(nth), ngi.PGs.get_facets_ptr(nth));
  
  ngon_type one_ph(pg, true);
  
  ngon_type::compact_to_used_nodes(one_ph.PGs, crd);
  
  K_FLD::IntArray cnto;
  one_ph.export_to_array(cnto);
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  
  delete f; delete cn;
  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::checkCellsClosure(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt), outer, remaining;

  err = ngon_type::check_phs_closure(ngi);

  delete f; delete cn;

#ifdef E_DOUBLEINT
    return Py_BuildValue("l", long(err));
#else
    return Py_BuildValue("i", err);
#endif
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::extrudeUserDefinedBC(PyObject* self, PyObject* args)
{
  PyObject *arr, *pgs;

  if (!PyArg_ParseTuple(args, "OO", &arr, &pgs)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

  std::vector<E_Int> PGlist;
  // Passing the specified wall pgs to the boolean to ignore cells that fall inside bodies
  {

    FldArrayI* inds=NULL;
    E_Int res=0;
    if (pgs != Py_None)
      res = K_NUMPY::getFromNumpyArray(pgs, inds, true);

    std::auto_ptr<FldArrayI> pL(inds); // to avoid to call explicit delete at several places in the code.
  
    //std::cout << "result for NUMPY is : " << res << std::endl;
    if ((res == 1) && (inds != NULL)  && (inds->getSize() != 0))
    {
      E_Int nb_special_pgs = inds->getSize();
      //E_Int minid(INT_MAX), maxid(-1);
      PGlist.resize(nb_special_pgs);
      for (size_t i = 0; i < nb_special_pgs; ++i) 
      {
        PGlist[i]=(*inds)[i]-1;
        //std::cout << pgsList[i] << std::endl;
        //minid = std::min(minid, pgsList[i]);
        //maxid = std::max(maxid, pgsList[i]);
      }

      //std::cout << "min/max : " << minid << "/" << maxid << std::endl;
    }
  }

  ngi.flag_externals(INITIAL_SKIN);
  bool has_been_reversed;
  DELAUNAY::Triangulator dt;
  err = ngon_type::reorient_skins(dt, crd, ngi, has_been_reversed);
  //std::cout << "reversed ? " << has_been_reversed << std::endl;
  if (!err)
  {
    err = ngon_type::append_with_ghost_cells(crd, ngi, PGlist);
    //std::cout << "append_with_ghost_cells status : " << err << std::endl;
  }

  PyObject* tpl = NULL;

  if (!err)
  {
    K_FLD::IntArray cnto;
    ngi.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }
  
  delete f; delete cn;

  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::reorientExternalFaces(PyObject* self, PyObject* args)
{
  PyObject *arr;

  if (!PyArg_ParseTuple(args, "O", &arr)) return NULL;

  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString, *eltType;
  // Check array # 1
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;

  //~ std::cout << "crd : " << crd.cols() << "/" << crd.rows() << std::endl;
  //~ std::cout << "cnt : " << cnt.cols() << "/" << cnt.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngio(cnt), outer, remaining;

  bool has_been_reversed;
  DELAUNAY::Triangulator t;
  err = ngon_type::reorient_skins(t, crd, ngio, has_been_reversed);

  if (has_been_reversed)
    std::cout << "reorientExternalFaces : external faces has been reoriented" << std::endl;
  else
    std::cout << "reorientExternalFaces : external faces orientation is correct" << std::endl;

  PyObject* tpl = NULL;

  if (!err)
  {
    K_FLD::IntArray cnto;
    ngio.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crd, varString, cnto, 8, "NGON", false);
  }
  
  delete f; delete cn;

  return tpl;
}

//=============================================================================
/* XXX */
//=============================================================================
PyObject* K_INTERSECTOR::diffMesh(PyObject* self, PyObject* args)
{
  PyObject *arr1, *arr2;

  if (!PyArg_ParseTuple(args, "OO", &arr1, &arr2)) return NULL;

  K_FLD::FloatArray *f1(0), *f2(0);
  K_FLD::IntArray *cn1(0), *cn2(0);
  char *varString1, *varString2, *eltType1, *eltType2;
  // Check array # 1
  E_Int err = check_is_NGON(arr1, f1, cn1, varString1, eltType1);
  if (err) return NULL;
  // Check array # 2
  err = check_is_NGON(arr2, f2, cn2, varString2, eltType2);
  if (err) return NULL;

  std::auto_ptr<K_FLD::FloatArray> pf1(f1), pf2(f2);   //for memory cleaning
  std::auto_ptr<K_FLD::IntArray> pcn1(cn1), pcn2(cn2); //for memory cleaning

  K_FLD::FloatArray & crd = *f1;
  K_FLD::IntArray & cnt = *cn1;
  K_FLD::FloatArray & crd2 = *f2;
  K_FLD::IntArray & cnt2 = *cn2;

  // std::cout << "crd1 : " << crd.cols() << "/" << crd.rows() << std::endl;
  // std::cout << "cnt1 : " << cnt.cols() << "/" << cnt.rows() << std::endl;
  // std::cout << "crd2 : " << crd2.cols() << "/" << crd2.rows() << std::endl;
  // std::cout << "cnt2 : " << cnt2.cols() << "/" << cnt2.rows() << std::endl;

  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  ngon_type ng(cnt), ng1(cnt), ng2(cnt2);

  E_Int nb_cells1 = ng.PHs.size();
  E_Int nb_cells2 = ng2.PHs.size();

  // concatenate the meshes
  E_Int shft = crd.cols();
  ng2.PGs.shift(shft);
  ng.append(ng2);
  crd.pushBack(crd2);

  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(crd);

  // detect identical cells
  std::vector<E_Int> nids; // cell ids in the concatenated ng
  E_Int nb_match = ng.detect_phs_with_same_centroid<acrd_t> (acrd, nids);

  //std::cout << "detect_phs_with_same_centroid : " << nb_match << std::endl;
  //
  if (nb_match == 0)
    std::cout << "the meshes are totally unmatching" << std::endl;

  std::vector<bool> keep1(nb_cells1, true);
  std::vector<bool> keep2(nb_cells2, true);
  //std::vector<bool> *pKa, *pKb;

  //
  for (size_t i=0; i < nids.size(); ++i)
  {
    E_Int nid = nids[i];
    if (nid == i) continue;

    // 2 cells are matching
    if (i < nb_cells1)
      keep1[i] = false;
    else
      keep2[i-nb_cells1]=false;

    if (nid < nb_cells1)
      keep1[nid] = false;
    else
      keep2[nid-nb_cells1]=false;
  }

  PyObject *l(PyList_New(0)), *tpl;

  {
    ngon_type ng1o;
    std::vector<E_Int> pgnids, phnids;
    ng1.select_phs(ng1, keep1, phnids, ng1o);
    ng1o.remove_unreferenced_pgs(pgnids, phnids);

    K_FLD::FloatArray crdtmp(crd);
    ngon_type::compact_to_used_nodes(ng1o.PGs, crdtmp); //reduce points
    //export to numpy
    K_FLD::IntArray cnto;
    ng1o.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crdtmp, varString1, cnto, 8, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  {
    ngon_type ng2o;
    std::vector<E_Int> pgnids, phnids;
    ng2.select_phs(ng2, keep2, phnids, ng2o);
    ng2o.remove_unreferenced_pgs(pgnids, phnids);

    K_FLD::FloatArray crdtmp(crd);
    ngon_type::compact_to_used_nodes(ng2o.PGs, crdtmp); //reduce points
    //export to numpy
    K_FLD::IntArray cnto;
    ng2o.export_to_array(cnto);
    tpl = K_ARRAY::buildArray(crdtmp, varString1, cnto, 8, "NGON", false);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  return l;
}

//=======================  Intersector/PolyMeshTools/utils.cpp ====================
