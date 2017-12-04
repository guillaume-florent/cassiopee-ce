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

#ifndef _DELAUNAY_MESHDATA_H_
#define _DELAUNAY_MESHDATA_H_

#include "Def/DefContainers.h"

namespace DELAUNAY
{
  struct MeshData
  {
    typedef E_Int                        size_type;
    typedef K_CONT_DEF::int_vector_type  int_vector_type;
    typedef K_CONT_DEF::bool_vector_type bool_vector_type;
    
    MeshData():pos(0), connectB(0){};

    MeshData(K_FLD::FloatArray& p, const K_FLD::IntArray& cB):pos(&p),connectB(&cB)
    {
    }
    
    void clear()
    {
      hardNodes.clear();
      hardEdges.clear();
      connectM.clear();
      neighbors.clear();
      ancestors.clear();
      colors.clear();
      metrics.clear();
      mask.clear();
    }

    K_FLD::FloatArray*         pos;
    const K_FLD::IntArray*           connectB;
    int_vector_type            hardNodes;
    K_CONT_DEF::non_oriented_edge_set_type hardEdges;
    K_FLD::IntArray            connectM;
    K_FLD::IntArray            neighbors;
    int_vector_type            ancestors;
    int_vector_type            colors;
    K_FLD::FloatArray          metrics;
    bool_vector_type           mask;
  };

  template <typename SurfaceType>
  struct SurfaceMeshData : public MeshData
  {
    SurfaceMeshData(K_FLD::FloatArray& p2D, const K_FLD::FloatArray& p3D, const K_FLD::IntArray& cB, const SurfaceType& s)
      : MeshData(p2D, cB), surface(s), pos3D(p3D)
    {
    }

    const SurfaceType&  surface;
    K_FLD::FloatArray   pos3D;
  };
}

#endif

