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
# include "Nuga/Boolean/Splitter.h"

//#include <iostream>

using namespace std;
using namespace NUGA;


//=============================================================================
/* Agglomerate superfuous faces (overdefined polyhedra) */
//=============================================================================
PyObject* K_INTERSECTOR::splitNonStarCells(PyObject* self, PyObject* args)
{
  PyObject *arr(0);
  E_Float PH_conc_threshold(1./3.);
  E_Float PH_cvx_threshold(0.05);
  E_Float PG_cvx_threshold(1.e-8);

  if (!PYPARSETUPLEF(args, "Oddd", "Offf", &arr, &PH_conc_threshold, &PH_cvx_threshold, &PG_cvx_threshold)) return NULL;

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
  ngon_type ngi(cnt), ngo;
  
  NUGA::Splitter::split_non_star_phs<DELAUNAY::Triangulator>(crd, ngi, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold, ngo);

  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  
  PyObject* tpl = K_ARRAY::buildArray(crd, varString, cnto, -1, eltType, false);;

  delete f; delete cn;
  return tpl;
}


//=======================  Intersector/PolyMeshTools/splitCells.cpp ====================
