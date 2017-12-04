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

#ifndef __KCORE_SEARCH_BBTREE_H__
#define __KCORE_SEARCH_BBTREE_H__

#include "Fld/DynArray.h"
#include "Fld/ArrayAccessor.h"
#include "Def/DefContainers.h"
#include "MeshElement/Triangle.h"

namespace K_SEARCH
{

///
template <E_Int DIM>
class BoundingBox{

public:
  /// Constructors
  BoundingBox(const E_Float *mB, const E_Float* MB)
  {for (E_Int i = 0; i < DIM; ++i){minB[i] = *(mB++); maxB[i] = *(MB++);}}

  ///
  BoundingBox(const K_FLD::FloatArray& pos, const E_Int* nodes, size_t n)
  {
    const E_Float* Pi;
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = K_CONST::E_MAX_FLOAT; maxB[i] = -K_CONST::E_MAX_FLOAT;}

    for (size_t i = 0; i < n; ++i)
    {
      Pi = pos.col(nodes[i]);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t>
  BoundingBox(const K_FLD::ArrayAccessor<CoordinateArray_t>& pos, const E_Int* nodes, size_t n)
  {
    this->compute(pos, nodes, n);
  }

  inline const E_Float& getMinB(E_Int axis) const {return minB[axis];}
  inline const E_Float& getMaxB(E_Int axis) const {return maxB[axis];}

  // Warning : strict comparison
  bool operator==(const BoundingBox& rhs){

    for (E_Int i = 0; i < DIM; ++i){
      if (minB[i] != rhs.minB[i]) return false;
      if (maxB[i] != rhs.maxB[i]) return false;
    }
    return true;
  }
  
  BoundingBox& operator=(const BoundingBox& rhs){

    for (E_Int i = 0; i < DIM; ++i){
      minB[i] = rhs.minB[i];
      maxB[i] = rhs.maxB[i];
    }
    return *this;
  }
  
  bool is_included(const BoundingBox& rhs){
    for (E_Int i = 0; i < DIM; ++i)
    {
      if (minB[i] < rhs.minB[i]) return false;
      if (maxB[i] > rhs.maxB[i]) return false;
    }
    
    return true;
  }

  static bool intersection(const BoundingBox& b1, const BoundingBox& b2, BoundingBox& b)
  {
    std::vector<E_Float> palma(4);
        
    for (E_Int i = 0; i < DIM; ++i)
    {
      if (b1.minB[i] > b2.maxB[i]) return false;
      if (b2.minB[i] > b1.maxB[i]) return false;
    }  

    for (E_Int i = 0; i < DIM; ++i)
    {
      palma[0] = b1.minB[i]; palma[1] = b1.maxB[i]; palma[2] = b2.minB[i]; palma[3] = b2.maxB[i];
      std::sort(palma.begin(), palma.end());
      b.minB[i] = palma[1];
      b.maxB[i] = palma[2];
    }

    return true;
  }

  ~BoundingBox()
  {
  }
  
  template <template <typename ELEM, typename ALLOC = std::allocator<ELEM> > class Vector>
  void compute(const K_FLD::FloatArray& pos, const Vector<E_Int> & indices)
  {
    const E_Float* Pi;
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = K_CONST::E_MAX_FLOAT; maxB[i] = -K_CONST::E_MAX_FLOAT;}

    for (size_t i = 0; i < indices.size(); ++i)
    {
      Pi = pos.col(indices[i]);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t, template <typename ELEM, typename ALLOC = std::allocator<ELEM> > class Vector>
  void compute
  (const K_FLD::ArrayAccessor<CoordinateArray_t>& pos, const Vector<E_Int> & indices)
  {
    E_Float Pi[DIM];
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = K_CONST::E_MAX_FLOAT; maxB[i] = -K_CONST::E_MAX_FLOAT;}

    for (size_t i = 0; i < indices.size(); ++i)
    {
      pos.getEntry(indices[i], Pi);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t>
  void compute
  (const K_FLD::ArrayAccessor<CoordinateArray_t>& pos)
  {
    E_Float Pi[DIM];
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = K_CONST::E_MAX_FLOAT; maxB[i] = -K_CONST::E_MAX_FLOAT;}

    for (E_Int i = 0; i < pos.size(); ++i)
    {
      pos.getEntry(i, Pi);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t>
  void compute
  (const CoordinateArray_t& pos, const E_Int* nodes, size_t n, E_Int index_start = 0)
  {
    const E_Float* Pi;
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = K_CONST::E_MAX_FLOAT; maxB[i] = -K_CONST::E_MAX_FLOAT;}

    for (size_t i = 0; i < n; ++i)
    {
      Pi = pos.col(nodes[i] - index_start);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t>
  void compute
  (const K_FLD::ArrayAccessor<CoordinateArray_t>& pos, const E_Int* nodes, size_t n)
  {
    E_Float Pi[DIM];
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = K_CONST::E_MAX_FLOAT; maxB[i] = -K_CONST::E_MAX_FLOAT;}

    for (size_t i = 0; i < n; ++i)
    {
      pos.getEntry(nodes[i], Pi);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }

public://fixme //private:
  BoundingBox(){for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = K_CONST::E_MAX_FLOAT; maxB[i] = -K_CONST::E_MAX_FLOAT;}}

  
public://fixme //private:
  E_Float minB[DIM];
  E_Float maxB[DIM];
};

typedef BoundingBox<2> BBox2D;
typedef BoundingBox<3> BBox3D;

///
template <E_Int DIM, typename BBoxType = BoundingBox<DIM> >
class BbTree {

  public: /** Typedefs */

    typedef           BbTree                                self_type;
    typedef           K_CONT_DEF::size_type                 size_type;
    typedef           K_FLD::IntArray                       tree_array_type;
   

#define BBTREE_ROWS 2 // the first row contains the node id, the second(third) contains the left(right) child column id.

  public: /** Constructors and Destructor */

    /// Builds a tree and inserts the boxes from begin to end.
    BbTree(const std::vector<BBoxType*>& boxes, E_Float tolerance=E_EPSILON);

    /// Destructor.
    ~BbTree()
    {
      if (_root_id != 0)
        for (size_t i = _root_id; i < _boxes.size(); ++i) delete _boxes[i];
    };

  public: /** Query methods */

    /// Returns all the boxes which overlap the input box.
    /** Warning: out is not cleared upon entry.*/
    void getOverlappingBoxes(const E_Float* minB, const E_Float* maxB, std::vector<size_type>& out) const;
    
    /// Returns true if an overlapping box is found (flase otherwise).
    bool hasAnOverlappingBox(const E_Float* minB, const E_Float* maxB) const;

    /// Returns all the boxes which intersects the input ray defined by P0 and P1.
    /** Warning: out is not cleared upon entry.*/
    void getIntersectingBoxes(const E_Float* P0, const E_Float* P1, std::vector<size_type>& out, E_Float tolerance=E_EPSILON, bool strict=false) const;
    
    /// Returns the global bounding box : the tree's bounding box
    const BBoxType* getGlobalBox(){return _boxes[_root_id];}
    
    ///
    static E_Bool boxesAreOverlapping(const BBoxType * bb1, const BBoxType * bb2, const E_Float& tol);
    
    ///
    static E_Bool box1IsIncludedinbox2(const BBoxType * bb1, const BBoxType * bb2, const E_Float& tol);

  private: /** Insertion method */

    /** Inserts the boxes in the tree when constructing the BbTree.
     *  Returns the id of the boxes bounding box.*/
    size_type __insert(const std::vector<size_type>::const_iterator& begin,
                       const std::vector<size_type>::const_iterator& end);

    /// Computes the bounding box of the input boxes.
    void __getBoxesBox(const std::vector<BBoxType*>&  boxes, const std::vector<size_type>::const_iterator& begin,
                       const std::vector<size_type>::const_iterator& end, BBoxType& box) const;

    ///
    inline E_Int __getLongestSideAxis(const BBoxType& box) const;

    ///
    void __getOverlappingBoxes(BBoxType* box, E_Int node, std::vector<size_type>& out) const;
    
    ///
    bool __hasAnOverlappingBox(BBoxType* box, E_Int node) const;

    ///
    void __getIntersectingBoxes(const E_Float* P0, const E_Float* P1, E_Int node, const E_Float& abstol, std::vector<size_type>& out) const;
    ///
    void __getIntersectingBoxesStrict(const E_Float* P0, const E_Float* P1, E_Int node, const E_Float& abstol, std::vector<size_type>& out, const BBoxType& boxSeg) const;

    ///
    inline E_Bool __boxIntersectRay(const BBoxType& box, const E_Float* P0, const E_Float* P1, const E_Float& abstol) const;
    ///
    inline E_Bool __boxIntersectSeg(const BBoxType& box, const E_Float* P0, const E_Float* P1, const E_Float& abstol, const BBoxType& boxSeg) const;

#ifdef E_TIME1
public:
  static double _append_tree ;
  static double _get_box_boxes;
  static double _get_longest_axis;
  static double _build_vectors;
  static double _fill_vectors;
#endif

 // private:
  public:
    /// BB tree
    tree_array_type          _tree;

    /// boxes
    std::vector<BBoxType*>   _boxes;

    /// tolerance
    E_Float                  _tolerance;

    /// Empty column for the tree.
    size_type                _newL[BBTREE_ROWS];

    /// 
    size_type                _root_id;
    
}; // End class BbTree

typedef BbTree<2, BBox2D> BbTree2D;
typedef BbTree<3, BBox3D> BbTree3D;

} // end namespace


#include "Search/BbTree.cxx"




#endif /* BBTREE_H_ */
