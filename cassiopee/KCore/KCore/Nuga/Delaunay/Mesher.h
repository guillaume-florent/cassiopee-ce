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

#ifndef __DELAUNAY_MESHER_H__
#define __DELAUNAY_MESHER_H__

#include "Kernel.h"
#include "Predicates.h"
#include "MeshData.h"
#include "MesherMode.h"
#include "Metric.h"
#include "Connect/MeshTool.h"
#include "Refiner.h"
#include "Search/KdTree.h"
#include "macros.h"
#include "MeshElement/Edge.h"
#include "chrono.h"
#include "Linear/DelaunayMath.h"
#include <list>
#include <deque>

#ifdef DEBUG_MESHER
#include "iodata.h"
#include <sstream>
#include "IO/io.h"
#endif
#include <iostream>

namespace DELAUNAY
{

  template <typename T, typename MetricType>
  class Mesher
  {

  public:
    typedef K_CONT_DEF::size_type                                        size_type;
    typedef K_CONT_DEF::int_set_type                                     int_set_type;
    typedef K_CONT_DEF::int_vector_type                                  int_vector_type;
    typedef K_CONT_DEF::bool_vector_type                                 bool_vector_type;
    typedef K_CONT_DEF::int_pair_type                                    int_pair_type;
    typedef K_CONT_DEF::int_pair_vector_type                             int_pair_vector_type;
    typedef K_CONT_DEF::non_oriented_edge_set_type                       non_oriented_edge_set_type;
    typedef K_MESH::Triangle                                             element_type;
    typedef K_MESH::Edge                                                 edge_type;
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray>                      coord_access_type;
    typedef K_SEARCH::KdTree<>                                           tree_type;
    typedef typename DELAUNAY::Kernel<T>                                 kernel_type;


  public:
    Mesher(MetricType& metric, const MesherMode& mode);
    ~Mesher(void);

    E_Int run (MeshData& data);

  private:

    E_Int initialize();

    E_Int triangulate();

    E_Int restoreBoundaries(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                            K_FLD::IntArray& neighbors, int_vector_type& ancestors);

    E_Int setColors(size_type Nbox, MeshData& data);

    E_Int finalize(MeshData& data, E_Int N0);

    E_Int refine();

    E_Int clean_data(MeshData& data, const bool_vector_type& mask);

  private:

    E_Int __compute_bounding_box(const int_vector_type& cloud, E_Float& minX, E_Float& minY, E_Float& maxX, E_Float& maxY);

    E_Int __getPipe(size_type N0, size_type N1, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
      const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, int_set_type& pipe,
      int_pair_vector_type& X_edges);
    
    E_Int __get_xedge_on_shell
    (size_type N0, size_type N1, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
     const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, E_Int& S, E_Int& n, E_Float &tolerance);

    E_Int __forceEdge(size_type N0, size_type N1, /*int_set_type& pipe,*/ int_pair_vector_type& Xedges, 
      const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, K_FLD::IntArray& neighbors,
      int_vector_type& ancestors);    

  private:

    MeshData*                   _data;

    const MesherMode            _mode;

    MetricType&                 _metric;

    K_CONNECT::MeshTool*        _tool;

    coord_access_type*          _posAcc;
    tree_type*                  _tree;

    kernel_type*                _kernel;

    int_set_type                _box_nodes;

    E_Int                       _err;

    E_Int                       _N0;
    
#ifdef DEBUG_MESHER
    public:
      bool dbg_flag;
#endif
  };


  template <typename T, typename MetricType>
  Mesher<T, MetricType>::Mesher(MetricType& metric, const MesherMode& mode)
    :_data(0), _mode(mode), _metric(metric),_tool(0), _posAcc(0), _tree(0)
  {
    std::srand(1);//to initialize std::rand and ensure repeatability
  }

  template <typename T, typename MetricType>
  Mesher<T, MetricType>::~Mesher(void)
  {
    if (_tool)
    {
      delete _tool;
      _tool = 0;
    }

    if (_posAcc)
    {
      delete _posAcc;
      _posAcc = 0;
    }

    if (_tree)
    {
      delete _tree;
      _tree = 0;
    }

    if (_kernel)
    {
      delete _kernel;
      _kernel = 0;
    }
  }

  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::run (MeshData& data)
  {
    _err = 0;
    _data = &data;

#ifdef E_TIME
    chrono c;
    c.start();
#endif

#ifdef E_TIME
    std::cout << "INIT" << std::endl;
#endif

    // Check data etc...
    _err = initialize();

#ifdef E_TIME
    std::cout << c.elapsed() << std::endl;
    c.start();
    std::cout << "TRIANGULATE" << std::endl;
#endif

    // Triangulate.
    _err = triangulate();
    
    if (_err)
    {
      std::cout << "error triangulating" << std::endl;
#ifdef DEBUG_MESHER
      MIO::write("err_tria.mesh", *_data->pos, *_data->connectB, "BAR");
#endif
      return _err;
    }

#ifdef E_TIME
    std::cout << c.elapsed() << std::endl;
    c.start();
    std::cout << "RESTORE BOUNDARIES" << std::endl;
#endif

    //DELAUNAY::iodata::read("out.mdat", data);
    // Restore missing hard edges.
    _err = restoreBoundaries(*data.pos, data.connectM, data.neighbors, data.ancestors);
    if (_err)
    {
      std::cout << "error restoring boundaries" << std::endl;
      return _err;
    }
      

#ifdef E_TIME
    std::cout << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_MESHER
    if (dbg_flag)
      MIO::write("triangulationC.mesh", *_data->pos, _data->connectM, "TRI", &_data->mask);
#endif

#ifdef E_TIME
    std::cout << "SET COLORS" << std::endl;
    c.start();
#endif

    // Find out the subdomains.
    _err = setColors(data.pos->cols()-1, data);//fixme : Nbox
    if (_err)
    {
      std::cout << "error setting colors" << std::endl;
      return _err;
    }

    if (_mode.mesh_mode == MesherMode::REFINE_MODE)
    {
#ifdef E_TIME
    c.start();
    std::cout << "REFINE" << std::endl;
#endif

    // Refine
    _err = refine();
    if (_err)
      return _err;
    }

#ifdef E_TIME
    std::cout << c.elapsed() << std::endl;
    c.start();
    std::cout << "FINALIZE MESH" << std::endl;
#endif

    _err = finalize(data, _N0);

#ifdef E_TIME
    std::cout << " " << c.elapsed() << std::endl;
    c.start();
#endif

    if (_err)
      std::cout << "error finalizing" << std::endl;
    
    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::initialize()
  {
    // Fast returns.
    if (_err)   return _err;

    const K_FLD::IntArray& connectB = *_data->connectB;
    size_type       COLS(connectB.cols());
    
    // Check that IDs in connect are referring to column in pos.

    // Check that in 2D, connectB and pos have 2 rows

    //etc...

    // Keep the minimum index in pos that can be changed.
    _N0 = _data->pos->cols();

    // Store the hard edges in a non oriented set.
    E_Int Ni, Nj;
    K_FLD::IntArray::const_iterator pS = connectB.begin();
    std::set<E_Int> hNodes;
    for (size_type i = 0; i < COLS; ++i, pS = pS+2)
    {
      Ni = *pS;
      Nj = *(pS+1);
      _data->hardEdges.insert(K_MESH::NO_Edge(Ni, Nj));
      hNodes.insert(Ni);
      hNodes.insert(Nj);
    }

    // Reset hard nodes to be consistent with hard edges
    hNodes.insert(ALL(_data->hardNodes)); // Append with the input hard nodes.
    _data->hardNodes.clear();
    for (int_set_type::const_iterator it = hNodes.begin(); it != hNodes.end(); ++it)
      _data->hardNodes.push_back(*it);

    // Build the initial mesh.
    E_Float minX, minY, maxX, maxY, L;
    __compute_bounding_box (_data->hardNodes, minX, minY, maxX, maxY);

    L = std::max(maxY-minY, maxX-minX);

    double factor = 0.1;

    minX -= factor*L;
    minY -= factor*L;
    maxX += factor*L;
    maxY += factor*L;

    E_Float c1[2] = {minX,minY};
    E_Float c2[2] = {maxX,minY};
    E_Float c3[2] = {maxX,maxY};
    E_Float c4[2] = {minX,maxY};

    size_type C1,C2,C3,C4;
    _data->pos->pushBack(c1, c1+2); C1 = _data->pos->cols()-1;
    _data->pos->pushBack(c2, c2+2); C2 = _data->pos->cols()-1;
    _data->pos->pushBack(c3, c3+2); C3 = _data->pos->cols()-1;
    _data->pos->pushBack(c4, c4+2); C4 = _data->pos->cols()-1;

    _box_nodes.insert(C1);
    _box_nodes.insert(C2);
    _box_nodes.insert(C3);
    _box_nodes.insert(C4);

    E_Int T1[3] = {C1,C2,C4};
    E_Int T2[3] = {C2,C3,C4};

    _data->connectM.pushBack (T1, T1+3);
    _data->connectM.pushBack (T2, T2+3);

    _data->ancestors.resize(_data->pos->cols(), E_IDX_NONE);
    _data->ancestors[C1] = _data->ancestors[C2] = _data->ancestors[C4] = 0;
    _data->ancestors[C3] = 1;


    size_type def = E_IDX_NONE;
    _data->neighbors.resize(3, 2, &def);
    _data->neighbors(0,0) = 1;
    _data->neighbors(1,1) = 0;

    std::vector<E_Int> indices;
    indices.push_back(C1);
    indices.push_back(C2);
    indices.push_back(C3);
    indices.push_back(C4);

    // KdTree initialisation
    _posAcc = new K_FLD::ArrayAccessor<K_FLD::FloatArray>(*_data->pos);
    _tree = new tree_type(*_posAcc, indices); //add the box nodes only.

    /*
    _data->ancestors.reserve(10*pos.cols());
    _data->connectM.reserve(3,10*pos.cols());
    _data->neighbors.reserve(3,10*pos.cols());
    _tree->reserve(10*pos.cols());// fixme : in triangulation mode (only har nodes)/ in refine mode (rough estimatio of the total)
    */

    _tool = new K_CONNECT::MeshTool(*_tree);

    _kernel = new kernel_type(*_data, *_tool);

    //fixme
    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::triangulate()
  {
    // Fast returns.
    if (_err) return _err;

    if (!_mode.do_not_shuffle)
      std::random_shuffle (ALL(_data->hardNodes));

    size_type nb_nodes(_data->hardNodes.size()), Ni;

    std::vector<E_Int> newIds;
    int unconstrained;
    for (size_type i = 0; (i < nb_nodes) && (_err == 0); ++i)
    {
      Ni = _data->hardNodes[i];

      _err = _kernel->insertNode(Ni, _metric[Ni], unconstrained);
      _tree->insert(Ni);  
#ifdef DEBUG_MESHER
      if (dbg_flag)
      {
        std::ostringstream o;
        o << "t_" << i << ".mesh";
        //K_CONVERTER::DynArrayIO::write(o.str().c_str(), _data->pos, _data->connectM, "TRI", &_data->mask);
        MIO::write(o.str().c_str(), *_data->pos, _data->connectM, "TRI", &_data->mask);
      }
#endif
    }

    if (_err) return _err;

#ifdef DEBUG_MESHER
    if (dbg_flag)
      MIO::write("triangulation0.mesh", *_data->pos, _data->connectM, "TRI", &_data->mask);
#endif

    clean_data(*_data, _data->mask); // Clean by removing invalidated elements and update the data structure. 


#ifdef E_TIME
    std::cout << "cav time        : "           << _kernel->cavity_time << std::endl;
    std::cout << "    init cav time   : "       << _kernel->init_cavity_time << std::endl;
    std::cout << "            base time   : "   << _kernel->_base_time << std::endl;
    std::cout << "            append time : "   << _kernel->_append_time << std::endl;
    std::cout << "    fix cav time    : "       << _kernel->fix_cavity_time << std::endl;
    std::cout << "    sort bound time : "       << _kernel->sorting_bound_time << std::endl;
    std::cout << "remesh time     : "           << _kernel->remesh_time << std::endl;
    std::cout << "inval time      : "           << _kernel->inval_time << std::endl;
    std::cout << std::endl;
#endif

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::restoreBoundaries
    (const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
    K_FLD::IntArray& neighbors, int_vector_type& ancestors)
  {

    std::vector<K_MESH::NO_Edge> missing_edges;
    int_set_type pipe;
    int_pair_vector_type Xedges;

    size_type cols(connect.cols());
    K_FLD::IntArray::const_iterator pS;

    // Get all the triangulation edges.
    non_oriented_edge_set_type all_edges;
    for (size_type j = 0; j < cols; ++j)
    {
      pS = connect.col(j);
      for (size_type i = 0; i < element_type::NB_NODES; ++i)
        all_edges.insert(K_MESH::NO_Edge(*(pS+i), *(pS + (i+1)%element_type::NB_NODES)));
    }

    // Get the missing hard edges.
    std::set_difference (ALL(_data->hardEdges), ALL(all_edges), std::back_inserter(missing_edges));

    size_type sz = (size_type)missing_edges.size();
    for (size_type i = 0; (i < sz) && !_err; ++i)
    {
#ifdef DEBUG_MESHER
      if (dbg_flag)
        MIO::write("triangulationi.mesh", pos, connect, "TRI", &(_data->mask));
      E_Int ni,nj;
      ni = missing_edges[i].node(0);
      nj = missing_edges[i].node(1);
#endif
      _err = __getPipe(missing_edges[i].node(0), missing_edges[i].node(1), pos, connect, neighbors,
        ancestors, pipe, Xedges);
      if (_err){
        std::cout << "error getting pipe : " << _err << std::endl;
      }

      if (!_err)
        _err = __forceEdge(missing_edges[i].node(0), missing_edges[i].node(1), /*pipe,*/ Xedges,
        pos, connect, neighbors, ancestors);
      if (_err){
        std::cout << "error forcing edge" << std::endl;
        //fixme : store this edge as not forceable (and carry on if there is a "non strict" mode...)
      }
    }

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::setColors(size_type Nbox, MeshData& data)
  {
    size_type cols(data.connectM.cols()), Sseed(E_IDX_NONE), S, Sn, Ni, Nj;
    K_FLD::IntArray::const_iterator pS;

    data.colors.clear();
    data.colors.resize(cols, E_IDX_NONE);

    // Get an element connected to the box.
    for (size_type i = 0; (i < cols) && (Sseed == E_IDX_NONE); ++i)
    {
      pS = data.connectM.col(i);
      for (size_type j = 0; (j < element_type::NB_NODES) && (Sseed == E_IDX_NONE); ++j)
      {
        if (*(pS+j) == Nbox)
          Sseed = i;
      }
    }

    if (Sseed == E_IDX_NONE) // Error
      return 1;

    std::vector<E_Int> cpool;
    size_type color = 0;
    size_type colored = 0;
    while (colored != cols)
    {
      cpool.push_back(Sseed);
      data.colors[Sseed] = color;
      ++colored;

      while (!cpool.empty())
      {
        S = cpool.back();
        cpool.pop_back();
        pS = data.connectM.col(S);

        for (size_type i = 0; i < element_type::NB_NODES; ++i)
        {
          Sn = data.neighbors(i, S);

          if (Sn == E_IDX_NONE)
            continue;

          if (data.colors[Sn] != E_IDX_NONE)
            continue;

          Ni = *(pS + (i+1) % element_type::NB_NODES);
          Nj = *(pS + (i+2) % element_type::NB_NODES);

          if (_data->hardEdges.find(K_MESH::NO_Edge(Ni, Nj)) != _data->hardEdges.end())
            continue;

          data.colors[Sn] = color;
          ++colored;
          cpool.push_back(Sn);
        }
      }

      ++color;
      if (colored != cols)
      {
        bool found = false;
        size_type i = 0;
        for (; (i < cols) && !found; ++i)
          found = (data.colors[i] == E_IDX_NONE);
        Sseed = i-1;
      }
    }

    if ((color > 2) && (_mode.remove_holes == true)) // detect eventual interior
    {
      std::vector<K_FLD::IntArray> connects;
      connects.resize(color);
      E_Int c;
      K_FLD::IntArray bound;
      // Split by color.
      for (E_Int Si = 0; Si < data.connectM.cols(); ++Si)
      {
        c = data.colors[Si];
        pS = data.connectM.col(Si);
        if (c >= 1)
          connects[c].pushBack(pS, pS+3);
      }

      // Store the hard edges in an oriented set.
      E_Int Ni, Nj;
      K_FLD::IntArray::const_iterator pS = data.connectB->begin();
      K_CONT_DEF::oriented_edge_set_type hard_edges;
      for (size_type i = 0; i < data.connectB->cols(); ++i, pS = pS+2)
      {
        Ni = *pS;
        Nj = *(pS+1);
        hard_edges.insert(K_MESH::Edge(Ni, Nj));
      }

      // Reset the color to zero for interior parts.
      E_Int nbc = (E_Int)connects.size(), c1, Si, nbound;
      bool invalid_color;
      for (c1 = 1; c1 < nbc; ++c1)
      {
        //DynArrayIO::write("part.mesh", data.pos, connects[c]);
        _tool->getBoundary(connects[c1], bound);
        //KDynArrayIO::write("bound.mesh", data.pos, bound);
        invalid_color = true;
        nbound = bound.cols();

        for (Si = 0; (Si < nbound) && invalid_color; ++Si)
          invalid_color &= (hard_edges.find(K_MESH::Edge(bound.col(Si))) == hard_edges.end());

        if (invalid_color)
          for (Si = 0; Si < data.connectM.cols(); ++Si)
            if (data.colors[Si] == c1)
              data.colors[Si] = 0;
      }
    }

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::refine()
  {
    K_FLD::FloatArray& pos = *_data->pos;

    _kernel->setConstraint(_data->hardEdges);

    std::vector<size_type> refine_nodes;

    coord_access_type posAcc(pos);
    tree_type filter_tree(posAcc);// fixme : should be only triangulation nodes.

    size_type Ni, nb_refine_nodes;

    Refiner<MetricType> saturator(_metric);

#ifdef E_TIME
    chrono c;
#endif

    float contrained;
    bool carry_on = false;
    do
    {
      carry_on = false;
#ifdef E_TIME
      c.start();
#endif
      saturator.computeRefinePoints(*_data, _box_nodes, _data->hardEdges, refine_nodes);

#ifdef E_TIME
      std::cout << "__compute_refine_points : " << c.elapsed() << std::endl;
      c.start();
#endif

      saturator.filterRefinePoints(*_data, _box_nodes, refine_nodes, filter_tree);

#ifdef E_TIME
      std::cout << "__filter_refine_points : " << c.elapsed() << std::endl;
      c.start();
#endif

      nb_refine_nodes = refine_nodes.size();
      carry_on = (nb_refine_nodes != 0);

      std::random_shuffle (ALL(refine_nodes));

      _data->ancestors.resize(_data->pos->cols(), E_IDX_NONE);

      for (size_type i = 0; (i < nb_refine_nodes) && !_err; ++i)
      {
        Ni = refine_nodes[i];
        _err = _kernel->insertNode(Ni, _metric[Ni], contrained);
        _tree->insert(Ni);
      }

      if (_err)
        return _err;

#ifdef E_TIME
      std::cout << "insertion : " << c.elapsed() << std::endl;
      c.start();
#endif

      this->clean_data(*_data, _data->mask);// Clean the current mesh by removing invalidated elements.

#ifdef E_TIME
      std::cout << "compacting : " << c.elapsed() << std::endl;
      std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
#endif
    }
    while (carry_on);

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::__compute_bounding_box
    (const int_vector_type& cloud, E_Float& minX, E_Float& minY, E_Float& maxX, E_Float& maxY)
  {
    minX = minY = K_CONST::E_MAX_FLOAT;
    maxX = maxY = - K_CONST::E_MAX_FLOAT;

    K_FLD::FloatArray::const_iterator pK;
    size_type nb_nodes = (size_type)cloud.size();

    for (size_type i = 0; i < nb_nodes; ++i)
    {
      pK = _data->pos->col(cloud[i]);
      minX = std::min(minX, *pK);
      maxX = std::max(maxX, *pK);
      minY = std::min(minY, *(pK+1));
      maxY = std::max(maxY, *(pK+1));
    }

    return _err;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::__getPipe
    (size_type N0, size_type N1, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
    const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, int_set_type& pipe, int_pair_vector_type& X_edges)
  {
    int_vector_type Ancs;
    size_type ni, S(E_IDX_NONE), n(E_IDX_NONE), count(0), nopp, nstart, Ni, Nj;
    K_FLD::IntArray::const_iterator pS;
    E_Bool done(false), intersect;
    E_Float tolerance(_tool->getTolerance())/*fixmetol*/, tol2(tolerance*tolerance)/*fixmetol*//*, u00, u01, u10, u11, bestu=-1.;*//*used in case of poor ancestors to pick the right one*/;
    size_type dim = pos.rows();

    pipe.clear();
    X_edges.clear();
    
    E_Int err=__get_xedge_on_shell(N0, N1, pos, connect, neighbors, ancestors, S, n, tolerance);
    if (err || S==E_IDX_NONE)
      return err;
    
    X_edges.push_back(int_pair_type(S,n));
    pipe.insert(S);

    while (!done)
    {
      count  = 0;
      nopp   = element_type::getOppLocalNodeId(S, n, connect, neighbors);
      nstart = (nopp+1)%element_type::NB_NODES;
        
      S = neighbors(n, S);
      pS = connect.col(S);
        
      for (size_type i = 0; (i < edge_type::NB_NODES) && !done; ++i)
      {
        ni  = (nstart+i) % element_type::NB_NODES;
        Ni = *(pS + (ni+1) % element_type::NB_NODES);
        Nj = *(pS + (ni+2) % element_type::NB_NODES);

        done = (Nj == N1); // End of pipe reached.

        intersect  = !done && (element_type::surface(pos.col(N0), pos.col(Ni), pos.col(N1), dim) > tol2);
        intersect &= (element_type::surface(pos.col(N0), pos.col(N1), pos.col(Nj), dim) > tol2);
        // intersect = !done && edge_type::intersect<2>(pos.col(N0), pos.col(N1), pos.col(Ni), pos.col(Nj), tolerance, true/*absolute*/,u00, u01, u10, u11, overlap);
        if (intersect)
        {
          ++count;
          n = ni;
          X_edges.push_back(int_pair_type(S,n));
        }
      }

      if (done && (count != 0))  // Error : we must have Nj = N1 at i = 0 when reaching the end of pipe.
      {err=2; break;}
      if (!done && (count != 1)) // Error : one node is on the edge N0N1.
      {err=3; break;}

      pipe.insert(S);
    }
    
    return err;
  }
  
  inline E_Int sign2D(const E_Float*P0, const E_Float*P1, const E_Float*P)
  {
    E_Float perpdot = (P0[1]-P1[1])*(P[0]-P0[0]) + (P1[0]-P0[0])*(P[1]-P0[1]);
    //return (perpdot < -E_EPSILON) ? -1 : (perpdot > E_EPSILON) ? 1 : 0;
    return (perpdot < 0.) ? -1 : (perpdot > 0.) ? 1 : 0;
  }
  
  inline E_Float proj(const E_Float*P0, const E_Float*P1, const E_Float*P)
  {
    return (P1[0]-P0[0])*(P[0]-P0[0]) + (P1[1]-P0[1])*(P[1]-P0[1]);
  }
  
  ///
  template <typename T, typename MetricType>
  E_Int
  Mesher<T, MetricType>::__get_xedge_on_shell
  (size_type N0, size_type N1, const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
   const K_FLD::IntArray& neighbors, const int_vector_type& ancestors, E_Int& Sx, E_Int& nx, E_Float &tolerance)
  {
    int_vector_type Ancs;
    size_type Si, ni, Ni, Nj;
    K_FLD::IntArray::const_iterator pS;
    //E_Bool done(false), intersect, overlap;
    E_Float dum1,dum2,dum3;
    E_Bool dum4;
            
    if (ancestors[N0] == E_IDX_NONE)
      return 4;
    
    const E_Float* P0=pos.col(N0);
    const E_Float* P1=pos.col(N1);

    pS = connect.col(ancestors[N0]);
    _tool->getAncestors(N0, connect, ancestors, neighbors, std::back_inserter(Ancs));
    
    // Get the ancestor(s) intersecting N0N1. Should be only one, numerical error otherwise.
    size_type AncNb = (size_type)Ancs.size();
    E_Int SIGNi(-2), SIGNip1(-2);
    E_Float ux0=-K_CONST::E_MAX_FLOAT, ux;
    for (size_type i = 0; i < AncNb; ++i)
    {
      Si = Ancs[i];
      pS = connect.col(Si);
      ni = element_type::getLocalNodeId(pS, N0);

      Ni = *(pS + (ni+1)%element_type::NB_NODES);
      Nj = *(pS + (ni+2)%element_type::NB_NODES);

      if ((Ni == N1) || (Nj == N1))
      {
        Sx=E_IDX_NONE;
        return 0; // N0N1 is already in, so return OK
      }
      
      if (neighbors(ni, Si) == E_IDX_NONE)
        continue;
        
      SIGNi=sign2D(P0,P1, pos.col(Ni));
      SIGNip1=sign2D(P0,P1, pos.col(Nj));
      
      if (SIGNi*SIGNip1 <= 0)
      {
        //ux = proj(P0,P1, pos.col(Nj));
        //bool intersect = 
        edge_type::intersect<2>(P0, P1, pos.col(Ni), pos.col(Nj), tolerance, true/*absolute*/,ux, dum1, dum2, dum3, dum4);
        
        if (SIGNi*SIGNip1 < 0 && ux > 0. && ux < 1.)
        {
          Sx=Si;
          nx=ni;
          return 0;
        }
      
        if (ux0 < ux && ux < 1.1 && ux > -1.1)
        {
          Sx=Si;
          nx=ni;
          ux0=ux;
        } 
      }
    }
    return 0;
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::__forceEdge
    (size_type N00, size_type N11, /*int_set_type& pipe,*/ int_pair_vector_type& Xedges,
    const K_FLD::FloatArray& pos, K_FLD::IntArray& connect, K_FLD::IntArray& neighbors,
    int_vector_type& ancestors)
  {
    size_type r, S, b, Sn, bn, Ni, Nj, Nk, Nl, S1, S2, dim(pos.rows()), b1, b2, Xnb;
    size_t nb_try(0);
    K_FLD::IntArray::iterator pS, pSn;
    E_Bool isConvex;
    E_Float tolerance(_tool->getTolerance()), s1, s2, tol2(tolerance*tolerance);
    const size_type& NB_NODES = element_type::NB_NODES;

    while (!Xedges.empty()  && (nb_try < 10*Xedges.size()))
    {
      // Choose randomly an intersecting edge.

      Xnb = (size_type)Xedges.size();
      r   = std::rand() % Xnb;//fixme : not a good random...

      int_pair_type& E = Xedges[r];

      S = E.first;
      b = E.second;

      pS = connect.col(S);
      Ni = *(pS + b);
      Nj = *(pS + (b+1) % NB_NODES);
      Nl = *(pS + (b+2) % NB_NODES);

      Sn = neighbors(b, S);

      pSn = connect.col(Sn);
      bn = element_type::getOppLocalNodeId(S, b, connect, neighbors);

      Nk = *(pSn + bn);

      // Convexity test : skip if the pair of element (S, Sn) is not convex (i.e. edge swapping is not doable).
      // Test done by visibility criterion.

      s1 = element_type::surface(pos.col(Ni), pos.col(Nj), pos.col(Nk), dim);
      s2 = element_type::surface(pos.col(Ni), pos.col(Nk), pos.col(Nl), dim);
      isConvex  = (s1 > tol2) && (s2 > tol2);

      if (!isConvex)
      {
        ++nb_try;
        continue;
      }
      nb_try = 0;

      // Neighbors to modify
      S1 = neighbors((b+1) % NB_NODES, S);
      S2 = neighbors((bn+1) % NB_NODES, Sn);
      b1 = element_type::getOppLocalNodeId(S, (b+1) % NB_NODES, connect, neighbors);
      b2 = element_type::getOppLocalNodeId(Sn, (bn+1) % NB_NODES, connect, neighbors);

      // Update elements S and Sn (connect)
      *(pS  + (b+2)  % NB_NODES) = Nk;
      *(pSn + (bn+2) % NB_NODES) = Ni;

      // Update the neighboring (neighbors)
      neighbors((b+1)  % NB_NODES, S)  = Sn;
      neighbors((bn+1) % NB_NODES, Sn) = S;
      if ((S1 != E_IDX_NONE) && (b1 != E_IDX_NONE))
        neighbors(b1, S1)              = Sn;
      neighbors(bn, Sn)                = S1;
      if ((S2 != E_IDX_NONE) && (b2 != E_IDX_NONE))
        neighbors(b2, S2)              = S;
      neighbors(b, S)                  = S2;

      // Update the ancestors.
      ancestors[Ni] = ancestors[Nj] = ancestors[Nk] = S;
      ancestors[Nl] = Sn;

      // Remove the edge.
      E = Xedges.back(); // compacting : put the last in the current for the next pass
      Xedges.pop_back(); // and remove the last.
      if (--Xnb == 0)
        return 0;

      // Update (if they exist) the intersecting edge references:
      // (Sn, (bn+1)) -> (S,b)
      // (S, (b+1))   -> (Sn, bn)
      for (size_type i = 0; i < Xnb; ++i)
      {
        int_pair_type& E = Xedges[i];
        if ((E.first == Sn) && (E.second == (bn+1) % NB_NODES))
        {
          E.first = S;
          E.second = b;
        }
        else if ((E.first == S) && (E.second == (b+1) % NB_NODES))
        {
          E.first = Sn;
          E.second = bn;
        }
      }

      // If the three nodes of the new S (or Sn) are in the same side of N0N1, remove S(or Sn).
      // Otherwise insert swapE in Xedges. 
      size_type count = 0;
      E_Float N0N1[2], V[2], v;
      const E_Float* P0 = pos.col(N00);
      K_FUNC::diff<2> (P0, pos.col(N11), N0N1);

      for (size_type i = 0; i < NB_NODES; ++i)
      {
        K_FUNC::diff<2> (P0, pos.col(*(pS+i)), V);
        K_FUNC::crossProduct<2> (N0N1, V, &v);

        if (v > 0.)
          ++count;
        else if (v < 0.)
          --count;
      }

      if (::abs(count) > 1)
      {
        //pipe.erase(S);
        continue;
      }

      count = 0;
      for (size_type i = 0; i < NB_NODES; ++i)
      {
        K_FUNC::diff<2> (P0, pos.col(*(pSn+i)), V);
        K_FUNC::crossProduct<2> (N0N1, V, &v);
        if (v > 0.)
          ++count;
        else if (v < 0.)
          --count;
      }

      if (::abs(count) > 1)
      {/*pipe.erase(Sn);*/}
      else
        Xedges.push_back(int_pair_type(S, (b+1) % NB_NODES));// Swapped edge.
    }// While loop

    if (Xedges.empty())
      return 0;
    else
    {
      return 1;//fixme : handle errors in a better way (nb_try...)
    }
  }

  ///
  template <typename T, typename MetricType>
  E_Int
    Mesher<T, MetricType>::clean_data(MeshData& data, const bool_vector_type& mask)
  {
    std::vector<size_type> newIds;

    // Remove invalidated elements
    K_FLD::IntArray::compact (_data->connectM, mask, newIds);
    //Update neighbors
    K_FLD::IntArray::compact (_data->neighbors, newIds);
    K_FLD::IntArray::changeIndices (_data->neighbors, newIds);
    
    //Update colors
    size_type cols(data.connectM.cols()), nb_nodes(data.ancestors.size());
    if (!_data->colors.empty())
    {
      size_type sz = (size_type)newIds.size(), newId;
      for (size_type i = 0; i < sz; ++i)
      {
        newId = newIds[i];
        if (newId != E_IDX_NONE)
          _data->colors[newId] = _data->colors[i];
      }
      _data->colors.resize(cols);
    }

    //Update mask
    _data->mask.clear();
    _data->mask.resize(_data->connectM.cols(), true);

    //Update ancestors
    for (size_type i = 0; i < nb_nodes; ++i)
    {
      if (_data->ancestors[i] != E_IDX_NONE)
        _data->ancestors[i] = newIds[_data->ancestors[i]];
    }
    
#ifdef E_DEBUG
    //std::cout << _data->connectM;
    //std::cout << _data->neighbors;
#endif

    return (cols - data.connectM.cols());
  }

  ///
  template <typename T, typename MetricType>
  E_Int
  Mesher<T, MetricType>::finalize(MeshData& data, E_Int N0)
  {
    size_type cols(data.connectM.cols());

    data.mask.resize(cols);

    for (size_type i = 0; i < cols; ++i) // mask box elements in top of invalidated ones.
    {
      if (data.mask[i])
        data.mask[i] = (data.colors[i] != 0);
    }

    clean_data(data, data.mask);

    if (cols && data.connectM.cols() == 0)// all the elements are masked
      return 1;

    std::vector<E_Int> new_IDs;
    _tool->compact_to_mesh(*data.pos, data.connectM, new_IDs, &N0); // Remove unused nodes above N0.

    // update ancestors.
    for (size_t i = 0; i < data.ancestors.size(); ++i)
      if (new_IDs[i] != E_IDX_NONE)
      data.ancestors[new_IDs[i]] = data.ancestors[i];
    data.ancestors.resize(data.pos->cols());

    return 0;
  }
}
#endif
