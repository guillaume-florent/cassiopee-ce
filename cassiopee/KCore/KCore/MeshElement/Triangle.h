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
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef __K_MESH_TRIANGLE_H__
#define __K_MESH_TRIANGLE_H__

#include <limits>
#include "Def/DefTypes.h"
#include "Def/DefCplusPlusConst.h"
#include "Def/DefContainers.h"
#include "MeshElement/Edge.h"
#include "Fld/ArrayAccessor.h"
#include "Linear/DelaunayMath.h"

#define IS_ZERO(a) (( (a) < E_EPSILON) && ( (a) > -E_EPSILON)) //fixmetol
#define ROUND(x) (((x<E_EPSILON) && (x>-E_EPSILON)) ? 0.: x)

namespace K_MESH
{
  class Triangle
  {
  public:
    typedef E_Int                 node_type;
    typedef Edge                  self_type;
    typedef K_CONT_DEF::size_type size_type;
    typedef K_MESH::NO_Edge       boundary_type;

  public:
    static const size_type NB_NODES;
    static const size_type NB_TRIS;

  public: /* Constructors, Destructor and operators*/

    ///
    Triangle(node_type n0, node_type n1, node_type n2)
    {_nodes[0] = n0; _nodes[1] = n1; _nodes[2] = n2;}
    ///
    template <typename NodeIterator>
    explicit Triangle(const NodeIterator pN)
    {_nodes[0] = *pN; _nodes[1] = *(pN+1); _nodes[2] = *(pN+2);}
    ///
    Triangle(void){_nodes[0] = _nodes[1] = _nodes[2] = E_IDX_NONE;}
    ///
    template <typename AConnect_t>
    inline Triangle (const AConnect_t& cnt, E_Int Ti) { cnt.getEntry(Ti, _nodes);}
    ///
    ~Triangle(void){}

    /// Assignment operator.
    //Triangle& operator=(const Triangle& T){_nodes[0] = T._nodes[0]; _nodes[1] = T._nodes[1]; _nodes[2] = T._nodes[2]; return *this;}

  public: /* Set and Get methods */

    /// Sets the triangle nodes.
    template <typename NodeIterator>
    inline void setNodes(const NodeIterator pN){_nodes[0] = *pN; _nodes[1] = *(pN+1); _nodes[2] = *(pN+2);}

    /// Set the edge nodes.
    inline void setNodes(node_type N0, node_type N1, node_type N2){_nodes[0] = N0; _nodes[1] = N1; _nodes[2] = N2;}
    
    ///
    inline const E_Int* nodes(){return _nodes;}

    /// Gets the i-th node.
    inline const node_type& node(const size_type& i) const {return _nodes[i];}

    /// Gets the i-th edge (required to abstract some algorithms to be valid both in 2D and 3D.
    template <typename BoundaryType>
    inline void getBoundary(E_Int n, BoundaryType& b) const {b.setNodes(_nodes[(n+1)%NB_NODES], _nodes[(n+2)%NB_NODES]);}
    
    template <typename BoundaryType>
    static inline void getBoundary(const E_Int* pS, E_Int n, BoundaryType& b) {b.setNodes(*(pS+(n+1)%NB_NODES), *(pS+(n+2)%NB_NODES));}
    
    /// Reverse the orientation
    inline void reverse() {std::swap(_nodes[1], _nodes[2]);}
    ///
    template <typename NodeIterator>
    static void reverse(NodeIterator pN) {std::swap(*(pN+1), *(pN+2));}

  public : /* Static functions related to edges */  

    ///
    static void getBoundary(const Triangle& T1, const Triangle&  T2, K_MESH::NO_Edge& b); 
    ///
    static void getBoundary(const Triangle& T1, const Triangle&  T2, E_Int& i1, E_Int& i2);

    /// Computes the surface.
    static E_Float surface(const E_Float* p1, const E_Float* p2, const E_Float* p3, size_type dim);

    template <E_Int DIM> inline
      static E_Float surface(const E_Float* p1, const E_Float* p2, const E_Float* p3);
    
    /// Computes the surface vector
    template <E_Int DIM> inline
    static void ndS(const E_Float* p1, const E_Float* p2, const E_Float* p3, E_Float* ndS);

    // Normal
    static void normal(const E_Float* p1, const E_Float* p2, const E_Float* p3, E_Float* normal);
    ///
    static void normal(const K_FLD::FloatArray& coord, const E_Int* pN, E_Float* normal);
    static void normal(const K_FLD::ArrayAccessor<K_FLD::FldArrayF>& coord, const E_Int* pN, E_Float* normal);
    ///
    static void isoG(const K_FLD::FloatArray& coord, const E_Int* pN, E_Float* G);
    static void isoG(const K_FLD::ArrayAccessor<K_FLD::FldArrayF>& coord, const E_Int* pN, E_Float* G);
    static void isoG(const E_Float* P0, const E_Float* P1, const E_Float* P2, E_Float* G);

    /// Returns the rank of N in the storage pointed by pK.
    inline static size_type getLocalNodeId(const size_type* pK, size_type N);

    /// Returns the opposite (i.e not shared) node of the n-th neighbor.
    static size_type getOppLocalNodeId(size_type K, size_type n, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors);

    /// Checks whether the input edge NiNj has the same orientation as one of the triangle's boundaries.
    /// Returns 0 if succeded. Returns -1 if the edge doesn't belong to the triangle.
    static E_Int getOrientation(const Triangle&  T, const E_Int& Ni, const E_Int& Nj, E_Bool& same_orient);

    // Returns the circum disc(ellipse) center in a (an)isotropic uniform field.
    template <typename MetricType> inline
      void
      circumDiskCenter 
      (const K_FLD::FloatArray& pos, const K_FLD::IntArray::const_iterator& pN, double& R2, E_Float* C, const MetricType& m) const;

    /// Computes the parameters of the projected point of P in the base (P0P1, P0P2).
    template <E_Int DIM>
    static void
      project(const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* P, E_Float* UV);

    /// Overload of the above method.
    template <E_Int DIM, typename CoordMatrix, typename NodeIterator>
    static void
      project(const CoordMatrix& pos, NodeIterator pS, const E_Float* P, E_Float* UV);

    ///
    template <typename CoordMatrix, typename NodeIterator>
    static E_Float
      minDistanceToPoint(const CoordMatrix& pos, NodeIterator pS, const E_Float* P, E_Float* UV, bool& inside);

    ///
    template <E_Int DIM>
    static inline void planeLineMinDistance
      (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1,
      E_Float tol, E_Bool tol_is_absolute,
      E_Float& lambda, E_Float * UV,
      E_Bool& parallel, E_Bool& coincident, E_Float& min_distance, bool strict=false);
    
    /// Implementation de l'algo rapide de test d'inetrsection T3-T3 en 3D (INRIA - Rapport de recherche N° 4488)
    template <E_Int DIM>
    static inline bool fast_intersectT3 
      (const E_Float* P1, const E_Float* Q1, const E_Float* R1, const E_Float* P2, const E_Float* Q2, const E_Float* R2);

    ///
    template <E_Int DIM>
    static bool intersect 
      (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1,
      E_Float tol, E_Bool tol_is_absolute, 
      E_Float& u00, E_Float& u01, E_Int& tx, E_Bool& overlap);

    ///
    template <E_Int DIM>
    static bool intersect 
      (const K_FLD::FloatArray& pos, E_Int i0,  E_Int i1,  E_Int i2, E_Int n0, E_Int n1, 
      E_Float tol, E_Bool tol_is_absolute, 
      E_Float& u00, E_Float& u01, E_Int* tx, E_Bool& overlap, E_Bool& coplanar);
    
    ///
    template <E_Int DIM>
    static bool intersectv2 
    (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1, 
    E_Float tol, E_Bool tol_is_absolute, 
    E_Float& u0, E_Float& u1, E_Int& tx, E_Bool& overlap, E_Bool& coincident);

    ///
    template <E_Int DIM>
    static bool intersect2 
      (const K_FLD::FloatArray& pos, E_Int i0,  E_Int i1,  E_Int i2, E_Int n0, E_Int n1, 
      E_Float tol, E_Bool tol_is_absolute, 
      E_Float& u00, E_Float& u01, E_Bool& overlap);

    template <E_Int DIM>
    static E_Float qualityG(const E_Float* P0, const E_Float* P1, const E_Float* P2);
    
    //dummy
    static inline void triangulate(const E_Int* nodes, E_Int* target){target[0]=nodes[0]; target[1]=nodes[1]; target[2]=nodes[2];}
    inline void triangulate(E_Int* target){target[0]=_nodes[0]; target[1]=_nodes[1]; target[2]=_nodes[2];}

   public: /* Operators */ // relevant only for sub class for which nodes are sorted.

     bool operator< (const Triangle& r) const
     {
       if (_nodes[0] < r._nodes[0])
         return true;
       if (!(r._nodes[0] < _nodes[0]) && (_nodes[1] < r._nodes[1])) // (l.node(0) == r.node(0)) AND (l.node(1) < r.node(1)) 
         return true;
       if (!(r._nodes[0] < _nodes[0]) && !(r._nodes[1] < _nodes[1]) && (_nodes[2] < r._nodes[2])) // (l.node(0) == r.node(0)) AND (l.node(1) == r.node(1)) AND (l.node(2) < r.node(2))
         return true;

       return false;
     }

  protected:
    E_Int _nodes[3];

  };

  /////////////////////////////////////////////////////////////////////////////
  
  /**
    Sub class designed for std::set and std::map class storage.
  */

#define MIN_ID(n0, n1, n2) (!(n0 > n1) && !(n0 > n2)) ? n0 : (!(n1 > n0) && !(n1 > n2)) ? n1 : n2; 
#define MAX_ID(n0, n1, n2) (!(n0 < n1) && !(n0 < n2)) ? n0 : (!(n1 < n0) && !(n1 < n2)) ? n1 : n2; 

  /// Non oriented triangle class (nodes are sorted by increasing id).
  class NO_Triangle : public Triangle
  {
     public:
    ///
    NO_Triangle(node_type n0, node_type n1, node_type n2){setNodes(n0, n1, n2);}
    ///
    template <typename NodeIterator>
    explicit NO_Triangle(const NodeIterator pN){setNodes(pN);}
    ///
    NO_Triangle():Triangle(){}
    ///
    ~NO_Triangle(void){}

  public: /* Set and Get methods */

    /// Sets the triangle nodes.
    template <typename NodeIterator>
    inline void setNodes(const NodeIterator pN)
    {
      _nodes[0] = MIN_ID(*pN, *(pN+1), *(pN+2));
      _nodes[2] = MAX_ID(*pN, *(pN+1), *(pN+2));
      _nodes[1] = ((_nodes[0] != *pN) && (_nodes[2] != *pN)) ? *pN : ((_nodes[0] != *(pN+1)) && (_nodes[2] != *(pN+1))) ? *(pN+1) : *(pN+2);
    }

    /// Set the edge nodes.
    inline void setNodes(node_type n0, node_type n1, node_type n2)
    {
      _nodes[0] = MIN_ID(n0, n1, n2);
      _nodes[2] = MAX_ID(n0, n1, n2);
      _nodes[1] = ((_nodes[0] != n0) && (_nodes[2] != n0)) ? n0 : ((_nodes[0] != n1) && (_nodes[2] != n1)) ? n1 : n2;
    }
  };

#define MIN_LOC_ID(n0, n1, n2) (!(n0 > n1) && !(n0 > n2)) ? 0 : (!(n1 > n0) && !(n1 > n2)) ? 1 : 2; 

  /// Oriented triangle class (lowest id is stored first, but cyclic-ordering is preserved).
  class O_Triangle : public Triangle
  {
  public:
    ///
    O_Triangle(node_type n0, node_type n1, node_type n2)
    {
      E_Int id = MIN_LOC_ID(n0, n1, n2);
      if (id == 1)
      {
        _nodes[0] = n1; _nodes[1] = n2; _nodes[2] = n0;
      }
      else if (id == 2)
      {
        _nodes[0] = n2; _nodes[1] = n0; _nodes[2] = n1;
      }
    }
    ///
    template <typename NodeIterator>
    explicit O_Triangle(const NodeIterator pN){setNodes(pN);}
    ///
    O_Triangle():Triangle(){}
    ///
    ~O_Triangle(void){}

  public: /* Set and Get methods */

    /// Sets the triangle nodes.
    template <typename NodeIterator>
    inline void setNodes(const NodeIterator pN)
    {
      E_Int id = MIN_LOC_ID(*pN, *(pN+1), *(pN+2));
      _nodes[0] = *(pN+id); _nodes[1] = *(pN+(id+1)%3); _nodes[2] = *(pN+(id+2)%3);
    }

    /// Set the edge nodes.
    inline void setNodes(node_type n0, node_type n1, node_type n2)
    {
      E_Int id = MIN_LOC_ID(n0, n1, n2);
      if (id == 0)
      {
        _nodes[0] = n0; _nodes[1] = n1; _nodes[2] = n2;
      }
      else if (id == 1)
      {
        _nodes[0] = n1; _nodes[1] = n2; _nodes[2] = n0;
      }
      else if (id == 2)
      {
        _nodes[0] = n2; _nodes[1] = n0; _nodes[2] = n1;
      }
    }
  };

  
  /////////////////////////////////////////////////////////////////////////////

  /// General case : Anisotropic. But approximated by a uniformed field taking m as the field value.
  template <typename MetricValueType> inline
    void
    Triangle::circumDiskCenter 
    (const K_FLD::FloatArray& pos, const K_FLD::IntArray::const_iterator& pS, double& R2,
    E_Float* C, const MetricValueType& m) const 
  {
    const E_Float* const  p1   = pos.col(*pS);
    const E_Float* const  p2   = pos.col(*(pS+1));
    const E_Float* const  p3   = pos.col(*(pS+2));

    E_Float m11 = m[0], m22 = m[2], m12 = m[1];
    R2 = 0.;

    assert (pos.rows() >= 2); // fixme : should be dim instead of rows.

    E_Float a11 = m11 * (p2[0] - p1[0]) + m12 * (p2[1] - p1[1]);
    E_Float a12 = m22 * (p2[1] - p1[1]) + m12 * (p2[0] - p1[0]);
    E_Float a21 = m11 * (p3[0] - p1[0]) + m12 * (p3[1] - p1[1]);
    E_Float a22 = m22 * (p3[1] - p1[1]) + m12 * (p3[0] - p1[0]);

    E_Float n1  = m11 * p1[0] * p1[0] + 2 * m12 * p1[0] * p1[1] + m22 * p1[1] * p1[1];
    E_Float d1  = m11 * p2[0] * p2[0] + 2 * m12 * p2[0] * p2[1] + m22 * p2[1] * p2[1] - n1;
    E_Float d2  = m11 * p3[0] * p3[0] + 2 * m12 * p3[0] * p3[1] + m22 * p3[1] * p3[1] - n1;

    K_LINEAR::DelaunayMath::resoLin(a11, a12, a21, a22, 0.5*d1, 0.5*d2, C[0], C[1]);

    E_Float u = p1[0] - C[0];
    E_Float v = p1[1] - C[1];

    if (C[0] != K_CONST::E_MAX_FLOAT)
      R2 = (m11 * u * u) + (2 * m12 * u * v) + (m22 * v * v);
  }

  // Iso specialization. The metric field value is irrelevant here.
  template <> inline
    void
    Triangle::circumDiskCenter 
    (const K_FLD::FloatArray& pos, const K_FLD::IntArray::const_iterator& pS, double& R2,
    E_Float* C, const E_Float& dummy) const 
  {
    const E_Float* const  p1   = pos.col(*pS);
    const E_Float* const  p2   = pos.col(*(pS+1));
    const E_Float* const  p3   = pos.col(*(pS+2));

    R2 = 0.;

    assert (pos.rows() >= 2); // fixme : should be dim instead of rows.

    E_Float a11 = p2[0] - p1[0];
    E_Float a12 = p2[1] - p1[1];
    E_Float a21 = p3[0] - p1[0];
    E_Float a22 = p3[1] - p1[1];

    E_Float n1  = K_FUNC::sqrNorm<2>(p1);
    E_Float d1  = 0.5 * (K_FUNC::sqrNorm<2>(p2) - n1);
    E_Float d2  = 0.5 * (K_FUNC::sqrNorm<2>(p3) - n1);

    K_LINEAR::DelaunayMath::resoLin(a11, a12, a21, a22, d1, d2, C[0],C[1]);

    if (C[0] != K_CONST::E_MAX_FLOAT)
      R2 = K_FUNC::sqrDistance(p1, C, 2);
  }

  ///
  template <E_Int DIM, typename CoordMatrix, typename NodeIterator>
  void
    Triangle::project(const CoordMatrix& pos, NodeIterator pS, const E_Float* P, E_Float* UV)
  {
    project<DIM>(pos.col(*pS), pos.col(*(pS+1)), pos.col(*(pS+2)), P, UV);
  }

  ///
  template <E_Int DIM>
  void
    Triangle::project(const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* P, E_Float* UV)
  {
    E_Float U[3], V[3], VP[3], Y[3], W[3], d;
    U[2] = V[2] = VP[2] = 0.;

    E_Float& u = UV[0] = K_CONST::E_MAX_FLOAT;
    E_Float& v = UV[1] = K_CONST::E_MAX_FLOAT;

    K_FUNC::diff<DIM>(P1, P0, U);
    K_FUNC::diff<DIM>(P2, P0, V);
    K_FUNC::diff<DIM>(P, P0, VP);

    K_FUNC::crossProduct<3>(U, V, W);
    d = K_FUNC::normalize<3>(W);

    if (d < E_EPSILON) // Degenerated triangle. //fixme error ?
      return;

    d = 1. / d;

    K_FUNC::crossProduct<3>(VP, W, Y); // Y = W ^ X

    // u = (V ^ X) / (U ^ V) = - (V.Y) * d  where X is the projection of VP on the triangle.
    // v = (U.Y) * d;

    u = -d * K_FUNC::dot<DIM>(V, Y);
    v =  d * K_FUNC::dot<DIM>(U, Y);
  }

  ///
  template <typename CoordMatrix, typename NodeIterator>
  E_Float
    Triangle::minDistanceToPoint(const CoordMatrix& pos, NodeIterator pS, const E_Float* P, E_Float* UV, bool& inside)
  {
    E_Float n[3], U[3], V[3], VP[3], d, tol(E_EPSILON), u, v, di;

    assert (pos.rows() >= 3); // fixme : should be dim instead of rows.

    inside = false;

    const E_Float* const  p1   = pos.col(*pS);
    const E_Float* const  p2   = pos.col(*(pS+1));
    const E_Float* const  p3   = pos.col(*(pS+2));

    K_FUNC::diff<3>(p2, p1, U);
    K_FUNC::diff<3>(p3, p1, V);
    K_FUNC::diff<3>(P, p1, VP);

    Triangle::project<3>(pos, pS, P, UV);

    u = UV[0];
    v = UV[1];

    if ((u == K_CONST::E_MAX_FLOAT) || (v == K_CONST::E_MAX_FLOAT))
      // Numerical error : fixme
      return K_CONST::E_MAX_FLOAT;

    // The projection is inside.
    if (((1. - u - v) >= -tol)  && 
      (u >= -tol) && (v >= -tol))
    {
      K_FUNC::crossProduct<3>(U,V,n);
      K_FUNC::normalize<3>(n);
      d = K_FUNC::dot<3>(VP, n);
      inside = true;
      return ::fabs(d);
    }

    // Find d as the minimum edge distance to P.

    // d(P, A)
    UV[1] = 0.;
    d = Edge::edgePointMinDistance<3>(p1, p2, P, UV[0]);

    // d(P, B)
    u = 0.;
    di = Edge::edgePointMinDistance<3>(p1, p3, P, v);

    if (di < d)
    {
      UV[0] = u;
      UV[1] = v;
      d = di;
    }

    //
    di = Edge::edgePointMinDistance<3>(p2, p3, P, v);
    u = 1. - v;

    if (di < d)
    {
      UV[0] = u;
      UV[1] = v;
      d = di;
    }
    return d;
  }

  ///
  template <E_Int DIM>
  inline void
    K_MESH::Triangle::planeLineMinDistance
    (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1,
    E_Float tol, E_Bool tol_is_absolute,
    E_Float& lambda, E_Float * UV,
    E_Bool& parallel, E_Bool& coincident, E_Float& min_distance, bool strict)
  {
    E_Float E[DIM], U[DIM], V[DIM], W[DIM], P0Q0[DIM], P0Q1[DIM], tol2(tol*tol);
    E_Float L2, x0, x1, dx, d;

    lambda = UV[0] = UV[1] = min_distance = K_CONST::E_MAX_FLOAT;
    coincident = parallel = false;

    K_FUNC::diff<DIM>(Q1, Q0, E);
    L2 = K_FUNC::sqrNorm<DIM>(E);

    if (L2 < E_EPSILON*E_EPSILON)
      return;

    if (tol_is_absolute)
      tol2 /= L2;

    K_FUNC::diff<DIM>(P1, P0, U);
    K_FUNC::diff<DIM>(P2, P0, V);
    K_FUNC::diff<DIM>(Q0, P0, P0Q0);
    K_FUNC::diff<DIM>(Q1, P0, P0Q1);

    K_FUNC::crossProduct<DIM>(U,V,W);
    d = K_FUNC::normalize<DIM>(W);

    x0 = K_FUNC::dot<DIM>(W, P0Q0);
    x1 = K_FUNC::dot<DIM>(W, P0Q1);
    x1 = (IS_ZERO(x1)) ? 0. : x1; //for robustness : if Q1 is lying on the plane, lambda = -x0/dx should be 1
    dx = x1-x0;

    parallel = (dx*dx < (tol2 * L2)); //checkme

    if (!parallel) //(1 intersection)
    {
      lambda = -x0/dx;
      if (strict)
      {
        E_Float tolR = ::sqrt(tol2);
        if ( (lambda <= -tolR) || (lambda >= 1. + tolR) )
        {
          return ;
        }
      }
      min_distance = 0.;

      E_Float IP[DIM], VP[DIM], Y[DIM];
      for (E_Int i = 0; i < DIM; ++i)
        IP[i] = Q0[i] + lambda * E[i];

      d = 1. / d;

      K_FUNC::diff<DIM>(IP, P0, VP);
      K_FUNC::crossProduct<DIM>(VP, W, Y); // Y = W ^ X
      // u = (V ^ X) / (U ^ V) = - (V.Y) * d  where X is the projection of VP on the triangle.
      // v = (U.Y) * d;
      UV[0] = -d * K_FUNC::dot<DIM>(V, Y);
      UV[1] =  d * K_FUNC::dot<DIM>(U, Y);
    }
    else //parallel.
    {
      x0 = (x0 < 0.) ? -x0 : x0;
      x1 = (x1 < 0.) ? -x1 : x1;
      min_distance = (x0 < x1) ? x0 : x1;
    }

    coincident = (parallel && (min_distance*min_distance < tol2 * L2));
  }


  ///
  template <>
  inline void
    K_MESH::Triangle::planeLineMinDistance<2>
    (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1,
    E_Float tol, E_Bool tol_is_absolute,
    E_Float& lambda, E_Float * UV,
    E_Bool& parallel, E_Bool& coincident, E_Float& min_distance, bool strict)
  {
    lambda = UV[0] = UV[1] = K_CONST::E_MAX_FLOAT;
    parallel = coincident = true;
    min_distance = 0.;
  }
  
#define SIGN(a) ((a < -E_EPSILON) ? -1 : ((a > E_EPSILON) ? 1 : 0))  
  
  inline E_Int get_region (const E_Float* M, const E_Float* p, const E_Float* q, const E_Float* r)
  {
    /*
                         \ R3(16)/
                          \    /
                           \ /
                           /\ r
                         /   \
             R13(32)   /      \
                     /         \  R23(8)
                   /       0    \
                 /               \
  - - - - - - -/------------------\ - - - - -
     R1(1)   / p      R12(2)      q\  R2(4)
           /                        \
    */
#define CROSS_PRD(a,b,c) ( (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]) )
#define POTENTIAL_X(a) ((a == 9)||(a==10)||(a==18)||(a==34)||(a==36)||(a==40)) //sum of the edge ends' region value
    
    E_Float a1 = CROSS_PRD(M, p, q);
    E_Float a2 = CROSS_PRD(M, q, r);
    E_Float a3 = CROSS_PRD(M, r, p);
    
    if ((a1 > -E_EPSILON) && (a2 > -E_EPSILON) && (a3 > -E_EPSILON)) // is inside
      return 0;
    
    E_Int s = SIGN(a1);
    if (s == SIGN(a2))
    {
      if (s == 1)
        return 32;
      else // -1
        return 4; 
    }
    else if (s == SIGN(a3))
    {
      if (s == 1)
        return 8;
      else // -1
        return 1; 
    }
    else // SIGN(a2) == SIGN(a3)
    {
      if (s == 1)
        return 16;
      else // -1
        return 2; 
    }
  }
  
  inline void solid_angle(E_Int reg, const E_Float* p, const E_Float* q, const E_Float* r, const E_Float*& s, const E_Float*&t)
  {
    /*
                           /\ 
                         / r \
                       /      \
                     /         \ 
           \       /       0    \    /
             \   /p             q\  /
               /-------------------/ 
               t\      R12(2)     /s
                  \              /   
                    \           /
                      \        /
                        \     /
                          \  /
                            \ P
    */
    switch (reg)
    {
      case 1/*R1*/  : s=q;t=r;break;
      case 2/*R12*/ : s=p;t=q;break;
      case 4/*R2*/  : s=r;t=p;break;
      case 8/*R23*/ : s=r;t=q;break;
      case 16/*R3*/ : s=p;t=q;break;
      case 32/*R13*/: s=p;t=r;break;
      default       : break;
    }
  }
  
  /// 2D
  template <>
  inline bool K_MESH::Triangle::fast_intersectT3<2>
  (const E_Float* P1, const E_Float* Q1, const E_Float* R1, const E_Float* P2, const E_Float* Q2, const E_Float* R2)
  {
    E_Int region[6];
    
    region[0] = get_region(P2, P1, Q1, R1);
    if (region[0] == 0) return true;
    region[1] = get_region(Q2, P1, Q1, R1);
    if (region[1] == 0) return true;
    region[2] = get_region(R2, P1, Q1, R1);
    if (region[2] == 0) return true;
    region[3] = get_region(P1, P2, Q2, R2);
    if (region[3] == 0) return true;
    region[4] = get_region(Q1, P2, Q2, R2);
    if (region[4] == 0) return true;
    region[5] = get_region(R1, P2, Q2, R2);
    if (region[5] == 0) return true;
 
    const E_Float *s=NULL, *t=NULL;
    if (POTENTIAL_X(region[0]+region[1])) //P2Q2
    {
      solid_angle(region[0], P1, Q1, R1, s, t);
      if ((CROSS_PRD(P2,s,Q2) > E_EPSILON) && (CROSS_PRD(P2,Q2,t) > E_EPSILON))
        return true;
    }
    if (POTENTIAL_X(region[1]+region[2])) //Q2R2
    {
      solid_angle(region[1], P1, Q1, R1, s, t);
      if ((CROSS_PRD(Q2,s,R2) > E_EPSILON) && (CROSS_PRD(Q2,R2,t) > E_EPSILON))
        return true;
    }
    if (POTENTIAL_X(region[2]+region[0])) //R2P2
    {
      solid_angle(region[2], P1, Q1, R1, s, t);
      if ((CROSS_PRD(R2,s,P2) > E_EPSILON) && (CROSS_PRD(R2,P2,t) > E_EPSILON))
        return true;
    }
    return false;
  }
  
  template<typename T>
  inline void permut1(T&a, T&b, T&c)
  { T tmp(b); b=a;a=c;c=tmp; }
  
  template<typename T>
  inline void permut2(T&a, T&b, T&c)
  { T tmp(b); b=c;c=a;a=tmp; }
  
  ///
  template <>
  inline bool K_MESH::Triangle::fast_intersectT3<3>
  (const E_Float* P1, const E_Float* Q1, const E_Float* R1, const E_Float* P2, const E_Float* Q2, const E_Float* R2)
  {
    
    // Does T1 intersect Pi2 ?
    E_Float d11 = K_FUNC::zzdet4(P2,Q2,R2,P1);
    E_Float d21 = K_FUNC::zzdet4(P2,Q2,R2,Q1);
    E_Float d31 = K_FUNC::zzdet4(P2,Q2,R2,R1);
        
    if ( (d11 < -E_EPSILON) && (d21 < -E_EPSILON) && (d31 < -E_EPSILON)) // P1,Q1,and R1 are strictly on the same side of Pi2 -> no intersection
      return false;
    if ( (d11 > E_EPSILON) && (d21 > E_EPSILON) && (d31 > E_EPSILON))    // P1,Q1,and R1 are strictly on the same side of Pi2 -> no intersection
      return false;
    
    if (IS_ZERO(d11) && IS_ZERO(d21) && IS_ZERO(d31)) // Pi1 == Pi2 => 2D intersection test
    {
      // Go to the Plane
      K_FLD::FloatArray P(3,3);
      E_Float *u(P.col(0)), *v(P.col(1)), *w(P.col(2));
      
      K_FUNC::diff<3>(Q1, P1, u);
      K_FUNC::diff<3>(R1, P1, v);
      K_FUNC::crossProduct<3>(u, v, w);
      K_FUNC::normalize<3>(w);
      K_FUNC::normalize<3>(u);
      K_FUNC::crossProduct<3>(w, u, v);
      K_FUNC::normalize<3>(v);
      
      K_FLD::FloatArray::inverse3(P);
      //do the transform
      E_Float p1[2], q1[2], r1[2], p2[2], q2[2], r2[2];
      for (size_t i = 0; i < 2; ++i)
      {
        p1[i] = q1[i] = r1[i] = p2[i] = q2[i] = r2[i] = 0.;
        for (size_t j = 0; j < 3; ++j)
        {
          p1[i] += P(i,j) * P1[j];
          q1[i] += P(i,j) * Q1[j];
          r1[i] += P(i,j) * R1[j];
          p2[i] += P(i,j) * P2[j];
          q2[i] += P(i,j) * Q2[j];
          r2[i] += P(i,j) * R2[j];
        }
      }
      //
      return fast_intersectT3<2>(p1, q1, r1, p2, q2, r2);
    }
    
    // T1 intersects Pi2. Check now if T2 intersects Pi1
    
    E_Float d12 = K_FUNC::zzdet4(P1,Q1,R1,P2);
    E_Float d22 = K_FUNC::zzdet4(P1,Q1,R1,Q2);
    E_Float d32 = K_FUNC::zzdet4(P1,Q1,R1,R2);
    
    if ( (d12 < -E_EPSILON) && (d22 < -E_EPSILON) && (d32 < -E_EPSILON)) // P2,Q2,and R2 are strictly on the same side of Pi1 -> no intersection
      return false;
    if ( (d12 > E_EPSILON) && (d22 > E_EPSILON) && (d32 > E_EPSILON))    // P2,Q2,and R2 are strictly on the same side of Pi1 -> no intersection
      return false;
    
    // T1 intersects Pi2 and T2 intersects Pi1
    
    //permut to have p1 (reps. p2) on positive top of Pi2(resp. Pi1)
    const E_Float *p1(P1),*q1(Q1),*r1(R1),*p2(P2),*q2(Q2),*r2(R2);
    // P1, Q1, R1
    E_Int s11(SIGN(d11)), s21(SIGN(d21)), s31(SIGN(d31));
        
    if (s11 == s21)
    {
      permut1(p1,q1,r1);//p1 = R1; q1 = P1; r1 = Q1;
      if (s31 == -1)
      {
        std::swap(q2,r2);
        std::swap(d22, d32);
      }
    }
    else if (s11 == s31)
    {
      permut2(p1,q1,r1);//p1 = Q1; q1 = R1; r1 = P1;
      if (s21 == -1)
      {
        std::swap(q2,r2);
        std::swap(d22, d32);
      }
    }
    else if (s21 == s31)
    {
      if (s11 == -1)
      {
        std::swap(q2,r2);
        std::swap(d22, d32);
      }
    }
    else // 3 different signs
    {
      // if s11 == 1, nothing to do
      
      if (s21 == 1)
        permut2(p1,q1,r1);//p1 = Q1; q1 = R1; r1 = P1;
      else if (s31 == 1)
        permut1(p1,q1,r1);//p1 = R1; q1 = P1; r1 = Q1;
    }
    
    // p2, q2, r2
    E_Int s12(SIGN(d12)), s22(SIGN(d22)), s32(SIGN(d32));
    if (s12 == s22)
    {
      permut1(p2,q2,r2); //p2 = r2; q2 = p2; r2 = q2;
      if (s32 == -1) std::swap(q1,r1);
    }
    else if (s12 == s32)
    {
      permut2(p2,q2,r2); //p2 = q2; q2 = r2; r2 = p2;
      if (s22 == -1) std::swap(q1,r1);
    }
    else if (s22 == s32)
    {
      if (s12 == -1) std::swap(q1,r1);
    }
    else // 3 different signs
    {
      // if s12 == 1, nothing to do
      
      if (s22 == 1)
        permut2(p2,q2,r2); //p2 = q2; q2 = r2; r2 = p2;
      else if (s32 == 1)
        permut1(p2,q2,r2); //p2 = r2; q2 = p2; r2 = q2;
    }
    
    return ((K_FUNC::zzdet4(p1,q1,p2,q2) < E_EPSILON) && (K_FUNC::zzdet4(p1,r1,r2,p2) < E_EPSILON));
  }

  // Triangle-Edge intersection. ==> used only in CompGeom/trianglesIntersection.cpp
  template <E_Int DIM>
  bool
    K_MESH::Triangle::intersect 
    (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1, 
    E_Float tol, E_Bool tol_is_absolute, 
    E_Float& u0, E_Float& u1, E_Int& tx, E_Bool& overlap) 
  {
    E_Float lambda, UV[2], min_d, eps(E_EPSILON), tol1(tol);
    E_Bool  parallel, coincident;

    tx = 0;

    overlap = false;
    u0 = u1 = K_CONST::E_MAX_FLOAT;

    planeLineMinDistance<DIM>(P0, P1, P2, Q0, Q1, tol, tol_is_absolute, lambda, UV, parallel, coincident, min_d);

    E_Float L = ::sqrt(K_FUNC::sqrDistance(Q0, Q1, DIM));
    if (tol_is_absolute)
    {
      if (L > E_EPSILON)
        tol1 /= L;
    }

    if (min_d > tol1 * L)
      return false;
    else if (!coincident) // The line intersect the plane at a unique point (3D case only)
    {
      u0  = lambda;
      return ( (u0 > -tol1) && (u0 < 1. + tol1) &&
        ((1. - UV[0] -  UV[1]) >= -tol1)  && 
        (UV[0] >= -tol1) && 
        (UV[1] >= -tol1) );
    }
    else                  // The line is lying on the plane (up to 2 intersections).
    {
      E_Int Xcount = 0;
      E_Float Xu[] = {K_CONST::E_MAX_FLOAT, K_CONST::E_MAX_FLOAT};

      E_Float        dum0, dum1;
      E_Bool         ovlap;
      const E_Float* pP[] = {P0, P1, P2};

      for (E_Int i = 0; (i < NB_NODES) && (Xcount < 2); ++i)
      {
        if (Edge::intersect<DIM>(pP[i], pP[(i+1)%3], Q0, Q1, tol, tol_is_absolute, dum0, dum1, u0, u1, ovlap))
        {
          if (ovlap)
          {
            Xcount = 2;
            Xu[0] = u0;
            Xu[1] = u1;
            tx = i<<2;
            break;
          }
          if ((u0 == 0.) || (u0 == 1.))                    // Just sharing a node.
            continue;
          if ((Xcount == 1) && (::fabs(Xu[0] - u0) < eps)) // Already taken into account.
            continue;

          Xu[Xcount++] = u0;
          tx += i<<2;
        }
      }

      if (Xcount < 2)
      {
        E_Float UV1[2];
        project<DIM>(P0, P1, P2, Q0, UV1);

        bool Q0_is_inside_T = (((1. - UV1[0] -  UV1[1]) >= -eps)  && 
          (UV1[0] >= -eps) && (UV1[1] >= -eps));

        if (Q0_is_inside_T)
          Xu[Xcount++] = 0.;
      }

      if (Xcount < 2)
      {
        E_Float UV2[2];
        project<DIM>(P0, P1, P2, Q1, UV2);

        bool Q1_is_inside_T = (((1. - UV2[0] -  UV2[1]) >= -eps)  && 
          (UV2[0] >= -eps) && (UV2[1] >= -eps));

        if (Q1_is_inside_T)
          Xu[Xcount++] = 1.;
      }

      overlap = (Xcount == 2);

      u0 = (Xu[0] < Xu[1]) ? Xu[0] : Xu[1];
      u1 = (Xu[1] > Xu[0]) ? Xu[1] : Xu[0];

      return (Xcount > 0);
    }
  }
  
  // Triangle-Edge intersection.  ==> used in maskGen
  template <E_Int DIM>
  bool
    K_MESH::Triangle::intersectv2 
    (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1, 
    E_Float tol, E_Bool tol_is_absolute, 
    E_Float& u0, E_Float& u1, E_Int& tx, E_Bool& overlap, E_Bool& coincident) 
  {
    E_Float lambda, UV[2], min_d, eps(E_EPSILON), tol12(tol*tol);
    E_Bool  parallel;

    tx = 0;

    coincident=overlap = 0;
    u0 = u1 = K_CONST::E_MAX_FLOAT;

    planeLineMinDistance<DIM>(P0, P1, P2, Q0, Q1, tol, tol_is_absolute, lambda, UV, parallel, coincident, min_d);

    E_Float L2 = K_FUNC::sqrDistance(Q0, Q1, DIM);
    if (tol_is_absolute)
    {
      if (min_d*min_d > tol12)
        return false;
      
      if (L2 > E_EPSILON*E_EPSILON)
        tol12 /= L2; //make it relative
    }
    else
    {
      if (min_d*min_d > tol12 * L2)
        return false;
    }

    if (!coincident) // The line intersect the plane at a unique point (3D case only)
    {
      E_Float tol1=::sqrt(tol12);
      u0  = lambda;
      return ( (u0 > -tol1) && (u0 < 1. + tol1) &&
        ((1. - UV[0] -  UV[1]) >= -tol1)  && 
        (UV[0] >= -tol1) && 
        (UV[1] >= -tol1) );
    }
    else                  // The line is lying on the plane (up to 2 intersections).
    {
      E_Int Xcount = 0;
      E_Float Xu[] = {K_CONST::E_MAX_FLOAT, K_CONST::E_MAX_FLOAT};

      E_Float        dum0, dum1;
      E_Bool         ovlap;
      const E_Float* pP[] = {P0, P1, P2};

      for (E_Int i = 0; (i < NB_NODES) && (Xcount < 2); ++i)
      {
        if (Edge::intersect<DIM>(pP[i], pP[(i+1)%3], Q0, Q1, tol, tol_is_absolute, dum0, dum1, u0, u1, ovlap))
        {
          if (ovlap)
          {
            Xcount = 2;
            Xu[0] = u0;
            Xu[1] = u1;
            tx = i<<2;
            break;
          }
          if ((::fabs(u0) < eps) || (::fabs(u0 - 1.) < eps))    // Just sharing a node.
            continue;
          if ((Xcount == 1) && (::fabs(Xu[0] - u0) < eps)) // Already taken into account.
            continue;

          Xu[Xcount++] = u0;
          tx += i<<2;
        }
      }

      if (Xcount < 2)
      {
        E_Float UV1[2];
        project<DIM>(P0, P1, P2, Q0, UV1);

        // it should be tol1 instead of eps but not necessary to compute tol1 here : edge interferences are caught
        // by calls to Edge::intersect, so at this point if it is inside, it is "really" inside (with no interf.)
        bool Q0_is_inside_T = (((1. - UV1[0] -  UV1[1]) > eps)  && 
          (UV1[0] > eps) && (UV1[1] > eps));

        if (Q0_is_inside_T)
          Xu[Xcount++] = 0.;
      }

      if (Xcount < 2)
      {
        E_Float UV2[2];
        project<DIM>(P0, P1, P2, Q1, UV2);

        // it should be tol1 instead of eps but not necessary to compute tol1 here : edge interferences are caught
        // by calls to Edge::intersect, so at this point if it is inside, it is "really" inside (with no interf.)
        bool Q1_is_inside_T = (((1. - UV2[0] -  UV2[1]) > eps)  && 
          (UV2[0] > eps) && (UV2[1] > eps));

        if (Q1_is_inside_T)
          Xu[Xcount++] = 1.;
      }

      overlap = (Xcount == 2);

      u0 = (Xu[0] < Xu[1]) ? Xu[0] : Xu[1];
      u1 = (Xu[1] > Xu[0]) ? Xu[1] : Xu[0];

      return (Xcount > 0);
    }
  }

  // Triangle-Edge intersection. // USED IN TRI_CONFORMIZER and therefore Booleans
  template <E_Int DIM>
  bool
    K_MESH::Triangle::intersect 
    (const K_FLD::FloatArray& pos, E_Int i0,  E_Int i1,  E_Int i2, E_Int n0, E_Int n1, 
    E_Float tol, E_Bool tol_is_absolute, 
    E_Float& u0, E_Float& u1, E_Int* tx, E_Bool& overlap, E_Bool& coplanar) 
  {
    E_Float lambda, UV[2], min_d, eps(E_EPSILON), tolR2(tol*tol);
    E_Bool  parallel;

    const E_Float* P0 = pos.col(i0);
    const E_Float* P1 = pos.col(i1);
    const E_Float* P2 = pos.col(i2);
    const E_Float* Q0 = pos.col(n0);
    const E_Float* Q1 = pos.col(n1);

    overlap = 0;
    u0 = u1 = K_CONST::E_MAX_FLOAT;
    tx[0] = tx[1] = 0;

    planeLineMinDistance<DIM>(P0, P1, P2, Q0, Q1, tol, tol_is_absolute, lambda, UV, parallel, coplanar, min_d, true/*strict*/);

    E_Float L2 = K_FUNC::sqrDistance(Q0, Q1, DIM);
    if (tol_is_absolute)
    {
      if (min_d*min_d > tolR2)
        return false;
      
      if (L2 > eps*eps)
        tolR2 /= L2; //make it relative
    }
    else
    {
      if (min_d*min_d > tolR2 * L2)
        return false;
    }

    if (!coplanar) // The line intersect the plane at a unique point (3D case only)
    {
      E_Float tolR=::sqrt(tolR2);
//      fixme : Probablement a reactiver pour faire notamment passer les calculs en mode soft (C3_T3)
//      E_Float IP[DIM];
//      for (size_t k=0; k < DIM; ++k)
//        IP[k]=(1.-lambda)*Q0[k] + lambda*Q1[k];
//        
////      E_Float s0 = K_MESH::Triangle::surface(IP, P0, P1, DIM);
////      E_Float s1 = K_MESH::Triangle::surface(IP, P1, P2, DIM);
////      E_Float s2 = K_MESH::Triangle::surface(IP, P2, P0, DIM);
////      
////      E_Float Li2 = K_FUNC::sqrDistance(P0, P1, DIM);
////      E_Float h02 = 4.*s0*s0/Li2;
//      E_Float l1;
//      E_Float h02bis = K_MESH::Edge::linePointMinDistance2<DIM>(P0,P1, IP, l1);
//      
////      Li2 = K_FUNC::sqrDistance(P1, P2, DIM);
////      E_Float h12 = 4.*s1*s1/Li2;
//      E_Float h12bis = K_MESH::Edge::linePointMinDistance2<DIM>(P1,P2, IP, l1);
//      
////      Li2 = K_FUNC::sqrDistance(P2, P0, DIM);
////      E_Float h22 = 4.*s2*s2/Li2;
//      E_Float h22bis = K_MESH::Edge::linePointMinDistance2<DIM>(P2,P0, IP, l1);
//
//      if (h02bis < tolR2*L2)
//        tx[0] += 1;
//      if (h12bis < tolR2*L2)
//        tx[0] += 2;
//      if (h22bis < tolR2*L2)
//        tx[0] += 4;
      if (::fabs(UV[0]) <= tolR) // ie. along P0P2, hence 4
        tx[0] += 4;
      if (::fabs(UV[1]) <= tolR) // ie. along P0P1, hence 1
        tx[0] += 1;
      if (::fabs(1. - UV[0] -  UV[1]) <= tolR)// ie. along P1P2, hence 2
        tx[0] += 2;
      
      if (tx[0] == 3 && ((i1==n0)|| (i1==n1))) //false X : they just share anode
        return false;
      if (tx[0] == 5 && ((i0==n0)|| (i0==n1))) //false X : they just share anode
        return false;
      if (tx[0] == 6 && ((i2==n0)|| (i2==n1))) //false X : they just share anode
        return false;
      
      u0  = lambda;
      bool is_x = ( (u0 > -tolR) && (u0 < 1. + tolR) &&
        ((1. - UV[0] -  UV[1]) >= -tolR)  && 
        (UV[0] >= -tolR) && 
        (UV[1] >= -tolR) );
      
//      //ROBUSTNESS : consistency with Mesher<T, MetricType>::__forceEdge
//      //ensure that the IP is indeed far from the edges otherwise the mesher would fail to restore a boundary
//      if (is_x && tx[0]==0)
//      {
////        E_Float IP[DIM];
////        for (size_t k=0; k < DIM; ++k)
////          IP[k]=(1.-lambda)*Q0[k] + lambda*Q1[k];
////        
////        E_Float U[3], V[3], W[3];
////        K_MESH::Triangle::normal(P0, P1, P2, W);
////        K_FUNC::diff<3>(P1, P0, U);
////        K_FUNC::normalize<3>(U);
////        K_FUNC::crossProduct<3>(W, U, V);
////
////        // Build the transformation matrix.
////        K_FLD::FloatArray P(3,3);
////        for (E_Int i = 0; i < 3; ++i)
////        {
////          P(i, 0) = U[i];
////          P(i, 1) = V[i];
////          P(i, 2) = W[i];
////        }
////        // Transform the working coordinates.
////        K_FLD::FloatArray::inverse3(P);
////        E_Float pP0[3], pP1[3], pP2[3], pIP[3];
////        
////        for (size_t j = 0; j < 3; ++j)
////          pP0[j]=pP1[j]=pP2[j]=pIP[j]=0.;
////        
////        for (size_t j = 0; j < 2; ++j)
////        {
////          for (size_t n = 0; n < 3; ++n)
////          {
////            pP0[j] += P(j, n) * P0[n] ;
////            pP1[j] += P(j, n) * P1[n] ;
////            pP2[j] += P(j, n) * P2[n] ;
////            pIP[j] += P(j, n) * IP[n] ;
////          }
////        }
////        
////        E_Float s = K_MESH::Triangle::surface(pP0, pP1, pIP, 2);
////        if (s <= tolR) //fixmetol < or <= ?
////          tx[0] +=1;
////        s = K_MESH::Triangle::surface(pP1, pP2, pIP, 2);
////        if (s <= tolR) //fixmetol < or <= ?
////          tx[0] +=2;
////        s = K_MESH::Triangle::surface(pP2, pP0, pIP, 2);
////        if (s <= tolR) //fixmetol < or <= ?
////          tx[0] +=4;    
//      }
      
      return is_x;
    }
    else  // The line is lying on the plane (up to 2 intersections). => WE THEN ENFORCE COPLANARITY FOR ROBUSTNESS
    {
      E_Int Xcount = 0;
      E_Float Xu[] = {K_CONST::E_MAX_FLOAT, K_CONST::E_MAX_FLOAT};

      E_Float        dum0, dum1;
      E_Bool         ovlap;
      E_Int          pI[] = {i0, i1, i2};

      for (E_Int i = 0; (i < NB_NODES) && (Xcount < 2); ++i)
      {
        if (Edge::intersect<DIM>(pos, pI[i], pI[(i+1)%3], n0, n1, tol, tol_is_absolute, dum0, dum1, u0, u1, ovlap, true/*enforce coplanarity*/))
        {
          if (ovlap)
          {
            Xcount = 2;
            Xu[0] = u0;
            Xu[1] = u1;
            tx[0]=tx[1] = (int)::pow(2.,(int)i);
            break;
          }

          if ((Xcount == 1) && (::fabs(Xu[0] - u0) < eps)) /*fixmetol*/     // Already taken into account.
            continue;

          Xu[Xcount] = u0;
          
          tx[Xcount++] += (int)::pow(2.,(int)i);
        }
      }
      
      E_Float tolR = ::sqrt(tolR2);

      if (Xcount < 2 && Xu[0] != 0.) //do not try Q0 again if it's taken into account
      {
        E_Float UV1[2];
        project<DIM>(P0, P1, P2, Q0, UV1);

        // it should be tol1 instead of eps but not necessary to compute tol1 here : edge interferences are caught
        // by calls to Edge::intersect, so at this point if it is inside, it is "really" inside (with no interf.)
        //tolR=::sqrt(tolR2);
        
        bool Q0_interfers_T = (((1. - UV1[0] -  UV1[1]) > -tolR)  && 
          (UV1[0] > -tolR) && (UV1[1] > -tolR));
        
        if (Q0_interfers_T)
        {
        
          bool Q0_is_inside_T = (((1. - UV1[0] -  UV1[1]) > tolR)  && 
                                (UV1[0] > tolR) && (UV1[1] > tolR));

          if (Q0_is_inside_T)
          {
            assert (tx[Xcount]==0);
            Xu[Xcount++] = 0.;
          }
          else if ((::fabs(1. - UV1[0] -  UV1[1]) < tolR) && (UV1[0] < 1.) && (UV1[1] < 1.)) //last 2 conditions means "not sharing a node"
          {
            tx[Xcount] +=2;
            Xu[Xcount++] = 0.;
          }
          else if ((::fabs(UV1[0]) < tolR) && (UV1[1] > 0.) && (UV1[1] < 1.)) //last 2 conditions means "not sharing a node"
          {
            tx[Xcount] +=4;
            Xu[Xcount++] = 0.;
          }
          else if ((::fabs(UV1[1]) < tolR) && (UV1[0] > 0.) && (UV1[0] < 1.)) //last 2 conditions means "not sharing a node"
          {
            tx[Xcount] +=1;
            Xu[Xcount++] = 0.;
          }
        }
      }

      if (Xcount < 2)
      {
        E_Float UV2[2];
        project<DIM>(P0, P1, P2, Q1, UV2);
        
        bool Q1_interfers_T = (((1. - UV2[0] -  UV2[1]) > -tolR)  && 
                              (UV2[0] > -tolR) && (UV2[1] > -tolR));

        if (Q1_interfers_T)
        {
          // it should be tol1 instead of eps but not necessary to compute tol1 here : edge interferences are caught
          // by calls to Edge::intersect, so at this point if it is inside, it is "really" inside (with no interf.)
          bool Q1_is_inside_T = (((1. - UV2[0] -  UV2[1]) > tolR)  && 
                                (UV2[0] > tolR) && (UV2[1] > tolR));

          if (Q1_is_inside_T)
          {
            assert (tx[Xcount]==0);
            Xu[Xcount++] = 1.;
          }
          else if ((::fabs(1. - UV2[0] -  UV2[1]) < tolR) && (UV2[0] < 1.) && (UV2[1] < 1.)) //last 2 conditions means "not sharing a node"
          {
            tx[Xcount] +=2;
            Xu[Xcount++] = 1.;
          }
          else if ((::fabs(UV2[0]) < tolR) && (UV2[1] > 0.) && (UV2[1] < 1.)) //last 2 conditions means "not sharing a node"
          {
            tx[Xcount] +=4;
            Xu[Xcount++] = 1.;
          }
          else if ((::fabs(UV2[1]) < tolR) && (UV2[0] > 0.) && (UV2[0] < 1.)) //last 2 conditions means "not sharing a node"
          {
            tx[Xcount] +=1;
            Xu[Xcount++] = 1.;
          }
        }
      }

      overlap = (Xcount == 2);
      
      if (Xu[0] < Xu[1])
      {
        u0=Xu[0];u1=Xu[1];
      }
      else
      {
        u0=Xu[1];u1=Xu[0];
        std::swap(tx[0], tx[1]);
      }

      return (Xcount > 0);
    }
  }

  //=============================================================================
  template <E_Int DIM>
  E_Float
    K_MESH::Triangle::qualityG(const E_Float* P0, const E_Float* P1, const E_Float* P2)
  {
    // Q = 4 * sqrt(3) * S / (L * p), S suface, L longest edge, p perimeter.

    E_Float S, p = 0., d, L = 0.;

    S = surface<DIM>(P0, P1, P2);

    d = ::sqrt(K_FUNC::sqrDistance(P1, P0, DIM));
    p += d;
    L = (L < d) ? d : L;

    d = ::sqrt(K_FUNC::sqrDistance(P2, P0, DIM));
    p += d;
    L = (L < d) ? d : L;

    d = ::sqrt(K_FUNC::sqrDistance(P2, P1, DIM));
    p += d;
    L = (L < d) ? d : L;

    d = 4. * K_CONST::E_SQRT3 * (S / (L * p));
    d = (d > 1.) ? 1. : d;

    return d;
  }

  //=============================================================================
  template <> inline
  E_Float
  K_MESH::Triangle::surface<3>
  (const E_Float* P0, const E_Float* P1, const E_Float* P2)
  {
    E_Float U[3], V[3], W[3];

    K_FUNC::diff<3>(P1, P0, U);
    K_FUNC::diff<3>(P2, P0, V);
    K_FUNC::crossProduct<3>(U, V, W);

    return 0.5 * ::sqrt(K_FUNC::sqrNorm<3>(W));
  }

  //=============================================================================
  template <> inline 
  E_Float
  K_MESH::Triangle::surface<2>
  (const E_Float* p1, const E_Float* p2, const E_Float* p3)
  {
    return 0.5 * ((p1[0]-p3[0])*(p2[1]-p3[1]) - (p1[1]-p3[1])*(p2[0]-p3[0]));
  }

  //=============================================================================
  inline void
  K_MESH::Triangle::normal
  (const E_Float* P0, const E_Float* P1, const E_Float* P2, E_Float* normal)
  {
    E_Float U[3], V[3];

    K_FUNC::diff<3>(P1, P0, U);
    K_FUNC::diff<3>(P2, P0, V);
    K_FUNC::crossProduct<3>(U, V, normal);
    K_FUNC::normalize<3>(normal);
  }
  
  //=============================================================================
  inline void K_MESH::Triangle::normal(const K_FLD::FloatArray& coord, const E_Int* pN, E_Float* normal)
  {
    E_Float U[3], V[3];
    const E_Float* P0 = coord.col(*pN);
    const E_Float* P1 = coord.col(*(pN+1));
    const E_Float* P2 = coord.col(*(pN+2));

    K_FUNC::diff<3>(P1, P0, U);
    K_FUNC::diff<3>(P2, P0, V);
    K_FUNC::crossProduct<3>(U, V, normal);
    K_FUNC::normalize<3>(normal);
  }
  
  inline
  void K_MESH::Triangle::normal
  (const K_FLD::ArrayAccessor<K_FLD::FldArrayF>& coord, const E_Int* pN, E_Float* normal)
  {        
    const E_Float* P = coord.array().begin(coord.posX(0));
    E_Int stride = coord.row_stride();
    const E_Float* P0 = &P[*pN];     // points to the x coord of first node.
    const E_Float* P1 = &P[*(pN+1)]; // points to the x coord of second node.
    const E_Float* P2 = &P[*(pN+2)]; // points to the x coord of third node.
    
    E_Float U[3], V[3];
    K_FUNC::diff<3>(P1, P0, stride, U);
    K_FUNC::diff<3>(P2, P0, stride, V);
    
    K_FUNC::crossProduct<3>(U, V, normal);
    K_FUNC::normalize<3>(normal);
  }
  
  //=============================================================================
  template <> inline 
  void K_MESH::Triangle::ndS<3>(const E_Float* P0, const E_Float* P1, const E_Float* P2, E_Float* ndS)
  {
    E_Float U[3], V[3];
    
    ndS[0]=ndS[1]=ndS[2]=0.;

    K_FUNC::diff<3>(P1, P0, U);
    K_FUNC::diff<3>(P2, P0, V);
    
    // prevent numerical error when computing cross product. fixme : should be done evrywhere a cross product or determinant is done ?
    for (size_t i=0; i < 3; ++i)
    {
      U[i]=ROUND(U[i]);
      V[i]=ROUND(V[i]);
    }
    
    K_FUNC::crossProduct<3>(U, V, ndS);
    
    ndS[0] *= 0.5;
    ndS[1] *= 0.5;
    ndS[2] *= 0.5;
  }
  
  //=============================================================================
  inline void K_MESH::Triangle::isoG(const K_FLD::FloatArray& coord, const E_Int* pN, E_Float* G)
  {
    const E_Float* P0 = coord.col(*pN);
    const E_Float* P1 = coord.col(*(pN+1));
    const E_Float* P2 = coord.col(*(pN+2));
    
    for (size_t i = 0; i < 3; ++i)
      G[i]=(P0[i]+P1[i]+P2[i])/3.;
  }
  
  inline void K_MESH::Triangle::isoG
  (const K_FLD::ArrayAccessor<K_FLD::FldArrayF>& coord, const E_Int* pN, E_Float* G)
  {
    const K_FLD::FldArrayF& array = coord.array();
    const E_Float* X[3];
    
    X[0] = array.begin(coord.posX(0));
    X[1] = array.begin(coord.posX(1));
    X[2] = array.begin(coord.posX(2));
        
    for (size_t i = 0; i < 3; ++i)
      G[i]=(X[i][*pN]+X[i][*(pN+1)]+X[i][*(pN+2)])/3.;
  }
  
  inline void K_MESH::Triangle::isoG(const E_Float* P0, const E_Float* P1, const E_Float* P2, E_Float* G)
  {
    for (size_t i = 0; i < 3; ++i)
      G[i]=(P0[i]+P1[i]+P2[i])/3.;
  }
  
  //=============================================================================
  K_CONT_DEF::size_type
  K_MESH::Triangle::getLocalNodeId(const size_type* pK, size_type N)
  {
    for (size_type n = 0; n < K_MESH::Triangle::NB_NODES; ++n)
    {
      if (*(pK+n) == N) return n;
    }
    assert (false);// should never get here.
    return E_IDX_NONE;
  }

} // End Namespace K_MESH
#endif

