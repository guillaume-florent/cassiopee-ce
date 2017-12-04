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

#include "MeshElement/Triangle.h"
#include <limits>

const E_Int K_MESH::Triangle::NB_NODES = 3;
const E_Int K_MESH::Triangle::NB_TRIS = 1;

//=============================================================================
E_Float K_MESH::Triangle::surface(const E_Float* p1, const E_Float* p2, 
                                  const E_Float* p3, size_type dim)
{
  if (dim == 2)
    return surface<2>(p1, p2, p3);
  if (dim == 3)
    return surface<3>(p1, p2, p3);

  return K_CONST::E_MAX_FLOAT;
}

//=============================================================================
K_CONT_DEF::size_type
K_MESH::Triangle::getOppLocalNodeId
(K_CONT_DEF::size_type K, K_CONT_DEF::size_type n,
 const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors)
{
  assert (n < NB_NODES);

  size_type Kadj = neighbors(n, K);

  if (Kadj == E_IDX_NONE) return E_IDX_NONE;

  K_FLD::IntArray::const_iterator pK = connect.col(K);
  size_type N = *(pK+(n+1)%NB_NODES);
  K_FLD::IntArray::const_iterator pKadj = connect.col(Kadj);
  size_type a = getLocalNodeId(pKadj, N);

  return (a+1)%NB_NODES;
}

//=============================================================================
void K_MESH::Triangle::getBoundary
(const Triangle&  T1, const Triangle&  T2, K_MESH::NO_Edge& b)
{
  K_MESH::NO_Edge Ei, Ej;
  for (E_Int i = 0; i < NB_NODES; ++i)
  {
    T1.getBoundary(i, Ei);
    for (E_Int j = 0; j < NB_NODES; ++j)
    {
      T2.getBoundary(j, Ej);
      if (Ei == Ej)
      {
        b = Ei; return;
      }
    }
  }
}

//=============================================================================
void K_MESH::Triangle::getBoundary
(const Triangle&  T1, const Triangle&  T2, E_Int& i1, E_Int& i2)
{
  K_MESH::NO_Edge Ei, Ej;
  for (i1=0; i1 < NB_NODES; ++i1)
  {
    T1.getBoundary(i1, Ei);
    for (i2=0; i2 < NB_NODES; ++i2)
    {
      T2.getBoundary(i2, Ej);
      if (Ei == Ej)
        return;
    }
  }
}

//=============================================================================
E_Int K_MESH::Triangle::getOrientation(const Triangle&  T, 
                                       const E_Int& Ni, const E_Int& Nj, 
                                       E_Bool& same_orient)
{
  same_orient = false;
 
  for (E_Int n = 0; n < NB_NODES; ++n)
  {
    if ((T._nodes[n] == Ni) && (T._nodes[(n+1)%NB_NODES] == Nj))
    {
      same_orient = true;
      return 0;
    }
    if ((T._nodes[n] == Nj) && (T._nodes[(n+1)%NB_NODES] == Ni))
      return 0;
  }

  return -1;
}
