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

#ifndef __NGON_BOOLEANOPERATOR_H__
#define	__NGON_BOOLEANOPERATOR_H__

//#define DEBUG_EXTRACT
//#define DEBUG_W_PYTHON_LAYER

#include "Fld/ngon_t.hxx"
#include "Search/BbTree.h"
#include "Search/ZoneExtractor.h"
#include "TRI_Conformizer.h"
#include "Connect/BARSplitter.h"
#include "Connect/merge.h"
#include "MeshElement/Polyhedron.h"
#include "Connect/MeshTool.h"
#ifdef FLAG_STEP
#include "chrono.h"
#endif

// For the T3Mesher
#include "Nuga/Delaunay/Triangulator.h"
//////////////////

#ifdef DEBUG_EXTRACT
#include "debug.h"
#endif
#ifdef DEBUG_W_PYTHON_LAYER
#include <fstream>
#include <sstream>
#endif
//#include "medit.hxx"

#if defined(DEBUG_BOOLEAN) && !defined(MANUAL)
#include "TRI_debug.h"
#include "NGON_debug.h"
#include <fstream>
//#include "TRI_BooleanOperator.h"
#include <iostream>
#endif

#define Vector_t std::vector
#define TEMPLATE_COORD_CONNECT template <typename Coordinate_t, typename Connectivity_t>
#define NGON_BOOLEAN_CLASS NGON_BooleanOperator<Coordinate_t,Connectivity_t>
#define NGON_BOOLEAN_FLD_SPECIALIZED NGON_BooleanOperator<K_FLD::FldArrayF,K_FLD::FldArrayI>
#define NGON_BOOLEAN_DYN_SPECIALIZED NGON_BooleanOperator<K_FLD::FloatArray,K_FLD::IntArray>
#define NGON_DBG NGON_debug<Coordinate_t, Connectivity_t>
#define TRI_DBG TRI_debug

#define ROUND(x) (((x<E_EPSILON) && (x>-E_EPSILON)) ? 0.: x)
#define robust_det4(a,b,c,d) K_FUNC::zzdet3(ROUND(d[0]-a[0]), ROUND(d[1]-a[1]), ROUND(d[2]-a[2]), ROUND(b[0]-d[0]), ROUND(b[1]-d[1]), ROUND(b[2]-d[2]), ROUND(c[0]-d[0]), ROUND(c[1]-d[1]), ROUND(c[2]-d[2]))

#define zSIGN(a) ( (a < 0.)? -1 : (a==0.) ? 0 : 1 )
#define zABS(a) ((a < 0) ? -a : a)

#define BLANKED 0.
#define VISIBLE 1.

namespace NUGA
{

/// Generalized Polyhedra Meshes Boolean Operator
TEMPLATE_COORD_CONNECT
class NGON_BooleanOperator {
public : 
    enum eInterPolicy { SOLID_NONE, SOLID_RIGHT, SURFACE_RIGHT, BOTH_SURFACE};
    enum eMergePolicy {PRESERVE_LEFT, PRESERVE_RIGHT/*, PRESERVE_X*/};
    enum eAggregation {NONE=0, CONVEX=1, FULL=2};
    enum eOperation   {INTER, DIFF, UNION, MODIFIED_SOLID};
    
    enum eRetCode {EMPTY_X=-7, OK=0,ERROR=1, IMMERSED_1};
    
    typedef ngon_t<Connectivity_t>               ngon_type;    
    typedef K_FLD::ArrayAccessor<Coordinate_t>   ACoordinate_t;
    typedef K_FLD::ArrayAccessor<Connectivity_t> AConnectivity_t;
    
private:
      enum eZone {Z_1=0, Z_2=1, Z_IN=2, /*Z_ALL=3,*/ Z_NONE=E_IDX_NONE};

public:
  /// Constructor with the 2 input GPM M1 & M2.
  //For FldArrays
  NGON_BooleanOperator(const K_FLD::FldArrayF& pos1, E_Int px, E_Int py, E_Int pz, const K_FLD::FldArrayI& cNGON1,
                       const K_FLD::FldArrayF& pos2, E_Int px2, E_Int py2, E_Int pz2, const K_FLD::FldArrayI& cNGON2,
                       E_Float tolerance, eAggregation aggtype);
  //For DynArrays
  NGON_BooleanOperator(const K_FLD::FloatArray& pos1, const K_FLD::IntArray& cNGON1,
                       const K_FLD::FloatArray& pos2, const K_FLD::IntArray& cNGON2,
                       E_Float tolerance, eAggregation aggtype);
  /// Destructor.
  ~NGON_BooleanOperator(void){}
  
public:
  /// 
  E_Int Intersection(Coordinate_t& coord, Connectivity_t& connect, eInterPolicy XPol = SOLID_RIGHT, eMergePolicy MergePol = PRESERVE_RIGHT);
  ///
  E_Int Diff(Coordinate_t& coord, Connectivity_t& connect, eInterPolicy XPol = SOLID_RIGHT, eMergePolicy MergePol = PRESERVE_RIGHT);
  /// 
  E_Int Union(Coordinate_t& coord, Connectivity_t& connect, eInterPolicy XPol = SOLID_RIGHT, eMergePolicy MergePol = PRESERVE_RIGHT);
  ///
  E_Int Modified_Solid(Coordinate_t& coord, Connectivity_t& connect, eMergePolicy MergePol = PRESERVE_RIGHT);
  ///
  E_Int Diffsurf(Coordinate_t& coord, Connectivity_t& connect);
  
public:
  ///
  void passPGs(E_Int rk, const Vector_t<E_Int>& PGlist){ //rk=0 for walls (xcellN), rk=1 for ghosts

    _pglist2[rk]=PGlist;
    
#ifdef DEBUG_W_PYTHON_LAYER
    std::ostringstream o;
    o << "PG_list_" << rk << ".txt";
  std::ofstream of(o.str().c_str());
  E_Int nb_pgs = PGlist.size();
  for (E_Int i=0; i < nb_pgs; ++i)
  {
    const E_Int& PGi = PGlist[i];
    of << PGi << std::endl;
  }
  of.close();
#endif
    
  }
  
  void nb_ghost(const ngon_type& ng)
  {
    if (_nb_cells2 == 0) return;
    
    E_Int count=0;
    for (E_Int i=0; i < ng.PHs.size(); ++i)
    {
      E_Int Anci = ng.PHs._ancEs(1,i);
      //if (ng.PHs._type[i] == GHOST)++count;   
      if (Anci >= _nb_cells2 && Anci != E_IDX_NONE) ++count;
    }
    std::cout << count << " ghosts exist" << std::endl;
  }
  
  void setTriangulatorParams(bool do_no_shuff, bool improve_qual)
  {
    _do_not_shuffle = do_no_shuff;
    _improve_quality_by_swapping = improve_qual;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Functions based on mapping between input and output meshes //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /// Conservative Interpolation for Chimera
  ///rec_ids[i] gives the i-th receiver, its donnors range from don_ids[xdon[i]] to don_ids[xdon[1+1]-1].
  // Same for coefs : from don_coefs[xdon[i]] to don_coefs[xdon[1+1]-1].
  E_Int volume_coefficients
       (E_Int rec_op /*0 or 1*/, std::vector<E_Int>& rec_ids, std::vector<E_Int>& xdon, std::vector<E_Int>& don_ids, std::vector<E_Float>& don_coefs);
  ///
  E_Int volume_and_centroid_coefficients
       (E_Int rec_op, std::vector<E_Int>& rec_ids, std::vector<E_Int>& xdon, std::vector<E_Int>& don_ids,
        std::vector<E_Float>& don_coefs, K_FLD::FloatArray& piece_centroids, K_FLD::FloatArray& crd, ngon_type& ngoper, std::vector<E_Int>& piece_ids);
  ///
  E_Int conservative_transfer_pwc
       (E_Int rec_op /*0 or 1*/, const K_FLD::FloatArray& fields_don, K_FLD::FloatArray& fields_rec, E_Int field = -1 /*-1 means all*/);
  ///
  E_Int conservative_transfer_pwl
       (E_Int rec_op /*0 or 1*/, const K_FLD::FloatArray& fields_don, const K_FLD::FloatArray& grads_don, K_FLD::FloatArray& fields_rec, E_Int field = -1 /*-1 means all*/);
  ///
  E_Int total_mass(const K_FLD::FloatArray& crd, const K_FLD::IntArray& ng_cnt, const K_FLD::FloatArray& fields, const std::vector<E_Float>& vols, std::vector<E_Float>& masses);
  
  ///
  E_Int XcellN (std::vector<E_Float>& cellN);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  public:
  ///
  static E_Int agglomerate_by_collapsing_PGs(K_FLD::FloatArray& crd, ngon_type& ng, E_Float vmin);
  static E_Int agglomerate_by_collapsing_PGs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids);
 
private:
  ///
  void __init(eInterPolicy XPol, eMergePolicy MPol, eOperation Op);
  ///
  eRetCode __compute();
  ///
  E_Int __compute_solid_modified_layer();
  ///
  eRetCode __process_intersections(ngon_unit& wPGs, E_Int& nb_pgs1, ngon_unit& extrawPGs, K_FLD::IntArray& connectT3, Vector_t<E_Int>& nT3_to_PG, Vector_t<E_Int> & priority);
  ///
  E_Int __process_duplicates(const ngon_unit&wPGs, K_FLD::IntArray& connectT3, Vector_t<E_Int>& nT3_to_PG, Vector_t<E_Int>& priority);
  ///
  void __refine_open_PGs(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& colors, ngon_unit& extrawPGs);
  ///
  eRetCode __get_working_PGs(eInterPolicy XPol, eMergePolicy MPol, ngon_unit& wPGs, E_Int& nb_pgs1, ngon_unit& extrawPGs);
  ///
  E_Int __conformize(K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, Vector_t<E_Int>& colors, E_Int X0, Vector_t<E_Int>& priority);
  ///
  eRetCode __focus_on_intersection_zone(eInterPolicy XPol, ngon_type& wNG1, ngon_type& wNG2, ngon_type& rNG1, ngon_type& rNG2);
  ///
  eRetCode __reorient_externals(eInterPolicy XPol, ngon_type& wNG1, ngon_type& wNG2);
  ///
  void __flag_PHs_sharing_nodes_with_selected_PGs(const ngon_type& wNG, Vector_t<bool>& keep, const Vector_t<E_Int>& pg_oids);
  
  ///
  void __refine_working_area(ngon_type& wNG, const Vector_t<E_Int>& indices, ngon_type& remainingNG);
  ///
  void __refine_working_area(ngon_type& wNG, const Vector_t<bool>& flag, ngon_type& remainingNG);
 
  ///
  void __create_ph_boxes(const ngon_type& wNG, const ACoordinate_t& coord, Vector_t<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BBox3D*& pool, K_SEARCH::BBox3D& GBbox);
  ///
  void __create_pg_boxes(const ngon_type& wNG, const ACoordinate_t& coord, Vector_t<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BBox3D*& pool, K_SEARCH::BBox3D& GBbox, Vector_t<E_Int>*oids=0, bool only_skin=false);
  ///
  void __compact_and_join(ngon_type& ngio, K_FLD::FloatArray& coord);
  ///
  void __destroy_tree(Vector_t<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BbTree3D *tree, K_SEARCH::BBox3D* pool);
  
  //// PH contruction methods ////
  ///
  E_Int __build_PHT3s(const K_FLD::FloatArray& coord,  const K_FLD::IntArray& connectT3, const K_FLD::IntArray& connectT3o, const Vector_t<E_Int>& is_skin,
                      std::map<E_Int, Vector_t<E_Int> >& PHT3s);
  ///
  E_Int __build_PHs(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& is_skinT3, const Vector_t<E_Int>& nT3_to_PG, ngon_type& ngout);
  
  ///
  E_Int __build_connect_hard(const K_FLD::FloatArray& coord, ngon_unit& extrawPGs, E_Int nb_pgs1, const K_FLD::IntArray& connectT3,
                            K_FLD::IntArray& connectHard, Vector_t<E_Int>& nT3_to_PG, Vector_t<E_Int>& is_skin);
  
  //// Zoning methods ////
  ///
  E_Int __classify_skin_PHT3s(std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift
#ifdef DEBUG_BOOLEAN
, const K_FLD::IntArray& connectT3o
#endif
);
  ///
  E_Int __classify_soft();
  ///
  void __set_skin_PHT3s_zones(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift, Vector_t<eZone>& zones
#ifdef DEBUG_BOOLEAN
, const K_FLD::IntArray& connectT3
#endif
);// T3 -> PHT3
  ///
  void __flag_border_OUT_zones(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift, eZone Z, E_Int col, Vector_t<eZone>& zones, Vector_t<eZone>& T3z);
  ///
  void __flag_border_IN_zones(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift, const Vector_t<eZone>& T3z, Vector_t<eZone>& zones);
  ///
  void __set_PH_history(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift, E_Int nb_pgs1, const K_FLD::IntArray& F2E, K_FLD::IntArray& ancPH, bool soft);

#ifndef DEBUG_BOOLEAN
 private:
#else
   public:
#endif
  ///
  E_Int __assemble_PGT3s(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& T3_to_PG, std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s);
  ///
  E_Int __split_non_connex_PGT3s(const K_FLD::IntArray& connectT3, std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s);
  ///
  E_Int __split_multiply_shared_PGT3s(const K_FLD::IntArray& connectT3, std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s);
  ///
  void __build_T3_to_PH(const K_FLD::IntArray& connectT3, std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s, K_FLD::IntArray& neighbors);
  ///
  template <eAggregation AGG>
  E_Int __aggregate_PHs(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& is_skinT3, const std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s, ngon_type& ngout)
  {
#ifdef FLAG_STEP
    chrono c;
    c.start();
    std::cout << "NGON Boolean : __aggregate_PHs (AGG MODE :" << AGG << ") ..." << std::endl;
#endif
#ifdef DEBUG_EXTRACT
  std::ostringstream xfname;
  xfname << /*MIO::wdir <<*/ "_extract_boxes.txt" ;
  std::ofstream extract_file (xfname.str().c_str());
#endif

      ngon_unit pgsi, pgs, phs;

      phs.reserve(PH_to_PGT3s.size());
      pgs.reserve(connectT3.cols());//fixme : is a better approx possible ?

      std::deque<E_Int> molec_ph;
      E_Int face(0), sz, key;
      bool skin_ph;

      Vector_t<E_Int> ph_bit;
      std::map <E_Int, Vector_t<E_Int> > key_to_phbits;
      std::map <E_Int, Vector_t<E_Int> >::const_iterator itPHbits;
      std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >::const_iterator it;
      std::map<E_Int, Vector_t<E_Int> >::const_iterator itPG;
      
#ifdef FLAG_STEP
  E_Int zcount=0;
#endif
#ifdef DEBUG_BOOLEAN
  E_Int faultyId = 86644;
#endif

      E_Int aphi = 0;
      for (it = PH_to_PGT3s.begin(); it != PH_to_PGT3s.end(); ++it, ++aphi)
      {
#ifdef FLAG_STEP
        ++zcount;
#endif
#ifdef DEBUG_BOOLEAN
        const E_Int& PHi = it->first;
        //std::cout << "PHi : " << PHi << std::endl;
#endif

#ifdef DEBUG_EXTRACT
        if (it->first==19571)
        {
          K_FLD::FloatArray crd;
          Gen_debug::get_PHT3_points(it->first, connectT3, PH_to_PGT3s, _coord, crd); // get the coordinates for this PH  
          Gen_debug::add_box_to_extract_file(crd, extract_file);
          extract_file.close();//just here because of the infinite loop : should be closed at the end
        }
#endif

        const std::map<E_Int, Vector_t<E_Int> >& Aggs = it->second;
        itPG = Aggs.begin();

        molec_ph.clear();

#ifdef DEBUG_BOOLEAN
          //K_FLD::IntArray cT3;
          //Vector_t<E_Int> colors;
#endif

          skin_ph = false;
          for (; itPG != Aggs.end(); ++itPG)
          {

#ifdef DEBUG_BOOLEAN
            const E_Int& PGi = itPG->first;
            //std::cout << "PGi : " << PGi << std::endl;
            _enabled = false;
            //if (phs.size()==faultyId && (PGi==170301 || PGi==170299))
              //_enabled=true;
#endif      

              const Vector_t<E_Int> & T3s = itPG->second;
              key = *std::min_element(T3s.begin(), T3s.end());

#ifdef DEBUG_BOOLEAN
              //if (key==10991)
              //_enabled=true;
              /*if (_enabled)
              {
                K_FLD::IntArray connectM;
                connectM.reserve(3, T3s.size());
                for (size_t i = 0; i < T3s.size(); ++i)
                  connectM.pushBack(connectT3.col(T3s[i]), connectT3.col(T3s[i])+3);  
                std::ostringstream o;
                o << "pgt3_" << PGi << ".mesh";
                MIO::write(o.str().c_str(), _coord, connectM, "TRI");
                cT3.pushBack(connectM);
                colors.resize(cT3.cols(), PGi);
              }*/
#endif

              const E_Int & is_skin = (is_skinT3[key] != INNER) && (is_skinT3[key] != CONNEXION_SKIN);
              skin_ph |= is_skin;

              itPHbits = key_to_phbits.find(key);

              if (itPHbits == key_to_phbits.end()) // Not already aggregated
              {
                pgsi.clear();
                __aggregate<AGG>(connectT3, T3s, pgsi);

#ifdef DEBUG_BOOLEAN
                /*if (_enabled)
                {
                  ACoordinate_t ac(_coord);
                  std::ostringstream o;
                  o << "pgagg_" << PGi;
                  NGON_DBG::write(o.str().c_str(), ac, pgsi);
                }*/
#endif

                // FLAGS AND HISTORY
                pgs.append(pgsi);
                pgs._type.resize(pgs.size(), is_skin ? INITIAL_SKIN : INNER);//IMPORTANT : use flag to mark skin PGs (for soft)

                E_Int hPG = _anc_PG[_nT3_to_oPG[key]];
                E_Int I = (hPG < _nb_pgs1) ? 0 : 1;
                E_Int anc[] = {E_IDX_NONE, E_IDX_NONE};
                anc[I]=hPG;
                for (E_Int kk = 0; kk < pgsi.size(); ++kk)
                  pgs._ancEs.pushBack(anc, anc+2);

                // Build the corresponding PH bit.
                sz = pgsi.size();
                ph_bit.resize(sz);
                for (E_Int k = 0; k < sz; ++k) ph_bit[k] = ++face;
                key_to_phbits[key] = ph_bit;
                itPHbits = key_to_phbits.find(key);
              }

              const Vector_t<E_Int>& phb = itPHbits->second;

#ifdef DEBUG_BOOLEAN
              size_t T3sz = T3s.size(); 
              size_t phbsz = phb.size(); 
              assert(!phb.empty());
              if (AGG == NONE)
                assert(T3sz == phbsz);
              else
                assert(phbsz <= T3sz);
#endif

              // Create the PH progressively.
              molec_ph.insert(molec_ph.end(), phb.begin(), phb.end());
          }

#ifdef DEBUG_BOOLEAN
         /*if (phs.size() == faultyId)
         {
           std::ostringstream o;
           o<< "pht3_" << PHi << ".mesh";
		   MIO::write(o.str().c_str(), _coord, cT3, "TRI", 0, &colors);
		 }*/
#endif

          const E_Int* pAnc = _anc_PH_for_PHT3s.col(aphi);
          molec_ph.push_front(molec_ph.size());// this is now a real PH molecule.
          phs.add(molec_ph);
          phs._type.push_back(skin_ph ? INITIAL_SKIN : INNER);
          phs._ancEs.pushBack(pAnc, pAnc + 2);
      }

      assert(phs.size() <= PH_to_PGT3s.size());

#ifdef DEBUG_BOOLEAN
      ///*E_Int countin(0), countex(0);
      //for (size_t i = 0; i < phs.size(); ++i)
      //{
      //  if (phs._external[i])
      //    ++countex;
      //  else
      //    ++countin;
      //}
      //std::cout << "in/ex : " << countin << "/" <<countex<< std::endl;*/
#endif

      ngout.PGs = pgs;
      ngout.PHs = phs;

#ifdef FLAG_STEP
      std::cout << "NGON Boolean : __aggregate_PHs : " << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_BOOLEAN
      if (ngout.PGs.size()*ngout.PHs.size())
      {
		  ngout.PGs.updateFacets();//fixme : required ?
		  ngout.PHs.updateFacets();

		  std::set<K_MESH::NO_Edge> free_edges;
		  for (size_t i = 0; i < ngout.PHs.size(); ++i)
		  {
			//std::cout << i << std::endl;
            bool ok = ngon_type::is_closed(ngout, i, free_edges);
#ifndef DEBUG_BOOLEAN
        assert(ok);
#else
        if (ok)
          continue;

        //NGON_DBG::draw_PH("pho.plt", _coord, ngout, i);

        std::cout << i << "-th PH is not closed !" << std::endl;
			  NGON_DBG::draw_PH("pho.tp", _coord, ngout, i);
        //if (i == 16)
        //{
        //  for (size_t j = 0; j < ngout.PHs.stride(i); ++j)
        //    NGON_DBG::draw_PGT3(_coord, ngout.PGs, ngout.PHs.get_facet(i, j) - 1);
        //}
        E_Int z = ngout.PHs.get_facet(16, 6);
        E_Int z1 = ngout.PHs.get_facet(16, 7);
        {
          K_FLD::IntArray connectE;
          std::set<K_MESH::NO_Edge>::const_iterator it;
          for (it = free_edges.begin(); it != free_edges.end(); ++it)
          {
            K_MESH::Edge e(it->node(0), it->node(1));
            connectE.pushBack(e.begin(), e.end());
          }
          connectE.shift(-1);
          MIO::write("free_edges.mesh", _coord, connectE, "BAR");
        }
#endif
		  }
	  }
	  ///*
	  //#ifdef FLDARR
	  //  K_FLD::FldArrayF crd;
	  //  coord.convert(crd);
	  //  K_FLD::ArrayAccessor<Coordinate_t> ac(crd, 1, 2, 3);
	  //#else
	  //  K_FLD::FloatArray & crd = _coord;
	  //  K_FLD::ArrayAccessor<Coordinate_t> ac(_coord);
	  //#endif    
	  //  NGON_DBG::write("aggregated", ac, ngout);  
	  //*/
#endif
      
#ifdef DEBUG_EXTRACT
  extract_file.close();
#endif

    return 0;
  }

  //// Neighbouring methods ////
  ///
  void __build_neighbors_table(const K_FLD::FloatArray& coord,  const K_FLD::IntArray& connectT3, const K_FLD::IntArray& connectT3o, K_FLD::IntArray& neighbors);
  ///
  void __sort_T3_sharing_an_edge(E_Int E0, E_Int E1, E_Int shift, 
                                 const K_FLD::FloatArray& normals, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, Vector_t<E_Int>& T3indices);
  ///
  void __update_neigboring(E_Int E0, E_Int E1, E_Int shift, const K_FLD::IntArray& connectT3,
                           const Vector_t<E_Int>& T3s, K_FLD::IntArray& neighbors);  
  //// PHT3 methods ////
  ///
  E_Int __assemble_PHT3s(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3,
                       const K_FLD::IntArray& neighbors, std::map<E_Int, Vector_t<E_Int> >& PHT3s);
  ///
  E_Int __remove_parasite_PHT3s(const Vector_t<E_Int>& is_skin, std::map<E_Int, Vector_t<E_Int> >& PHT3s, const K_FLD::IntArray& connectT3);
  ///
  E_Int __check_PHT3s_closure(bool and_manifoldness, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, std::map<E_Int, Vector_t<E_Int> >& PHT3s);
  
  //// Aggregation methods ////
  ///
  template <eAggregation AGG> inline 
  E_Int __aggregate(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs);
  ///
  inline E_Int __aggregate_none (const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs);
  ///
  inline E_Int __aggregate_full(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs);
  ///
  inline E_Int __aggregate_convex(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs);
  ///
  bool __is_convex(const K_CONT_DEF::int_pair_vector_type &boundaries, const std::map< E_Int, std::pair<E_Int, E_Int> >& node_to_nodes, 
                   const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const K_FLD::FloatArray& normals,
                   const Vector_t<E_Int>& globalId, std::deque<E_Int>& sorted_contour, E_Int& K0, E_Int& n0, E_Int& Eip1);
  ///
  E_Int __cut_mesh_convex_line(const K_FLD::IntArray& connectT3, const K_FLD::FloatArray& normals, const Vector_t<E_Int>& globalId,
                               const K_FLD::IntArray& connectB, E_Int K0, E_Int n0, E_Int Eip1,
                                     K_FLD::IntArray& neighbors, K_FLD::IntArray& connectT31, K_FLD::IntArray& connectT32);
  //// Orientation methods ////
  ///
  void __duplicate_orient(K_FLD::IntArray& connectT3);
  ///
  void __compute_normals(const K_FLD::ArrayAccessor<K_FLD::FloatArray>& acrd, const K_FLD::ArrayAccessor<K_FLD::IntArray>& acnt,  const ngon_unit& PGs, const Vector_t<E_Int> T3_to_PG, K_FLD::FloatArray& normals);
  /// to have the reversed T3s.
  void __swap(K_FLD::IntArray& arr);
  ///
  void __remove_orientations(std::map<E_Int, Vector_t<E_Int> >& PHT3s, E_Int shift);
  ///
  ///
  E_Int __discard_holes_by_box(const K_FLD::FloatArray& coord, ngon_type& wNG);
  E_Int __discard_prescribed_polygons(const K_FLD::FloatArray& coord, ngon_type& wNG, const Vector_t<E_Int>& PGlist);
  ///
  
  
  inline bool __is_ghost(E_Int i){ return (i >= _nb_cells2) ;} 
  
  void __remove_ghosts(ngon_type& ng)
  {
    if (_nb_cells2 == 0) return;
    
    ngon_unit realPHs;
    
    E_Int nb_phs = ng.PHs.size(), count(-1);
    Vector_t<E_Int> nids(nb_phs, E_IDX_NONE);
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      if (ng.PHs._ancEs(1,i) < _nb_cells2)
      {
        realPHs.add(ng.PHs.stride(i), ng.PHs.get_facets_ptr(i));
        nids[i] = ++count;
      }  
    }
    
    if (count == nb_phs) return;
    
    std::cout << "removing  " << nb_phs - count - 1 << " ghosts" << std::endl;
    
    realPHs.compact_attributes(ng.PHs, nids);
    ng.PHs = realPHs;
    
    Vector_t<E_Int> pgnids, phnids;
    ng.remove_unreferenced_pgs(pgnids, phnids);
    
  }
  
 private:  
  NGON_BooleanOperator();
  NGON_BooleanOperator(const NGON_BooleanOperator& orig);

 private:
  
  eInterPolicy _XPol;
  eMergePolicy _MPol;
  eAggregation _AggPol;
  eOperation   _Op;
  bool _build_hard;
  E_Float _tolerance;
  E_Float _convexity_tol;
  
private:private:
  bool _processed;// tell whether __create_boolean_zones() has been run or not.
  const Connectivity_t& _cNGON1;
  const Connectivity_t& _cNGON2;
  ACoordinate_t _aCoords1;
  Coordinate_t _crd2;
public:
  ngon_type _ngXs, _ngXh, _ng1, _ng2;
  ngon_type *_ngoper; // gather upon exit the result ngon (concatenation of some of _ngXs, _ngXh, _ng1, _ng2)
  
  //global output coordinates.
  K_FLD::FloatArray _coord;
  
  Vector_t<eZone> _zones;
  K_FLD::IntArray _anc_PH_for_PHT3s; //for soft
  K_FLD::IntArray _anc_PG;
  K_FLD::FloatArray _normals;
  
  Vector_t<E_Int> _nodes_history; //for modified solid's layer
  
  bool _do_not_shuffle, _improve_quality_by_swapping;
  
  Vector_t<std::pair<E_Float, E_Int> > _palmares;
  std::set<E_Int> _tmp_set_int;
  Vector_t<E_Int> _tmp_vec;
  std::deque<E_Int> _tmp_deq;
  K_FLD::IntArray _tmp_IntArray, _tmp_IntArray2;
  K_FLD::FloatArray _tmp_FloatArray;
  std::set<K_MESH::Edge> _tmp_set_oedge;
  std::map<E_Int, E_Int> _tmp_map_nn;
  
  Vector_t<E_Int> _nT3_to_oPG;
  E_Int _nb_pgs1;
  K_FLD::IntArray _F2E, _extraF2E;
  
  Vector_t<E_Int> _pglist2[2]; // 0 : WALL / 1 : GHOST EXTRUSION
  
  E_Int _nb_cells2;

#ifdef DEBUG_BOOLEAN
  bool _enabled;
#endif

};

typedef NUGA::NGON_BooleanOperator<K_FLD::FldArrayF, K_FLD::FldArrayI> Fld_NGON_BooleanOperator; 
typedef NUGA::NGON_BooleanOperator<K_FLD::FloatArray, K_FLD::IntArray> Dyn_NGON_BooleanOperator; 

/// Constructor for FldArrays
TEMPLATE_COORD_CONNECT
NGON_BOOLEAN_CLASS::NGON_BooleanOperator
(const K_FLD::FldArrayF& pos1, E_Int px, E_Int py, E_Int pz, const K_FLD::FldArrayI& cNGON1,
 const K_FLD::FldArrayF& pos2, E_Int px2, E_Int py2, E_Int pz2, const K_FLD::FldArrayI& cNGON2, E_Float tolerance, eAggregation aggtype)
: _tolerance(tolerance), _convexity_tol(1.e-2), _AggPol(aggtype), _processed(false),
 _cNGON1(cNGON1), _cNGON2(cNGON2), _aCoords1(pos1, px, py, pz), _crd2(pos2, px2, py2, pz2), _do_not_shuffle(true), _improve_quality_by_swapping(false)
{}

/// Constructor for DynArrays
TEMPLATE_COORD_CONNECT
NGON_BOOLEAN_CLASS::NGON_BooleanOperator
(const K_FLD::FloatArray& pos1, const K_FLD::IntArray& cNGON1,
 const K_FLD::FloatArray& pos2, const K_FLD::IntArray& cNGON2, E_Float tolerance, eAggregation aggtype)
  : _AggPol(aggtype), _tolerance(tolerance),  _convexity_tol(1.e-2), _processed(false),
 _cNGON1(cNGON1), _cNGON2(cNGON2), _aCoords1(pos1), _crd2(pos2), _do_not_shuffle(true), _improve_quality_by_swapping(false)
{}

///
TEMPLATE_COORD_CONNECT
NGON_BOOLEAN_CLASS::NGON_BooleanOperator() {}
///
TEMPLATE_COORD_CONNECT
NGON_BOOLEAN_CLASS::NGON_BooleanOperator(const NGON_BooleanOperator& orig) {}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::Intersection
(Coordinate_t& coord, Connectivity_t& connect,eInterPolicy XPol, eMergePolicy MPol)
{
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif

  eRetCode ret(OK);

  coord.clear();
  connect.clear();
  
  __init(XPol, MPol, INTER);

  if (!_processed) // process the transform.
    ret = __compute();

  if (ret == OK)
  {
    // color and split : 1, 2 or IN (coloring algo based on "external" info).
    E_Int err = __classify_soft();
    if (err)
      return err;

    coord = _coord;
    __compact_and_join(_ngXs, coord);
    _ngXs.export_to_array(connect);
    if (connect.cols() == 0)ret =ERROR;
    //connect should not be empty upon exit to silent the python, but in reality it can be (when it has 4 cols)
    //in this case we add this test for remapping use
    if (_ngXs.PGs.size() != 0)
      _ngoper=&_ngXs; // for mapping upon exit
  }

#ifdef DEBUG_BOOLEAN
  /*{
    E_Int nb_phs = _ngXs.PHs.size();
    K_FLD::IntArray connect;
    for (size_t i=0; i < nb_phs; ++i)
    {
      ngon_type ngo;
      Vector_t<bool> keep(nb_phs, false);
      keep[i]=true;
      ngon_type::select_phs(_ngXs, keep, ngo);
      connect.clear();
      ngo.export_to_array(connect);
      std::ostringstream o;
      o << "e_" << i << ".plt";
      MIO::write(o.str().c_str(), coord, connect, "NGON");
    }
  }*/
  
  /*std::cout << _ngXs.PHs.size() << std::endl;
  std::cout << _ngXs.PGs.size() << std::endl;
  NGON_DBG::draw_PGT3s(coord, _ngXs.PGs);*/
  
  if (ret != ERROR)
  {
    std::ostringstream o;
    o << "I";
    if (XPol == SOLID_RIGHT) o << "1";
    o << ".plt";
    MIO::write(o.str().c_str(), coord, connect, "NGON");
    std::cout << "total number of elements " << connect(0, 0) << std::endl;
  }
  else
    std::cout << "Intersection failed." << std::endl;
#endif

#ifdef FLAG_STEP
  if (connect.cols()) std::cout << "NGON Boolean : nb final PHs : "  << connect[2+connect[1]] << std::endl;
  std::cout << "NGON Booelan : total CPU : " << c.elapsed() << std::endl;
#endif

  return ret;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::Diff
(Coordinate_t& coord, Connectivity_t& connect, eInterPolicy XPol, eMergePolicy MPol)
{
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif

  eRetCode ret(OK);

  coord.clear();
  connect.clear();

  __init(XPol, MPol, DIFF);

  if (!_processed) // process the transform.
    ret = __compute();

  if (ret == OK)
  {
    // color and split : 1, 2 or IN (coloring algo based on "external" info).
    E_Int err = __classify_soft();
    if (err)
      return err;

    coord = _coord;
    __compact_and_join(_ng1, coord);
    _ng1.export_to_array(connect);
    _ngoper=&_ng1; // for mapping upon exit
  }
  else if (ret == EMPTY_X)
  {
    coord = _aCoords1.array();
    _ng1 = _cNGON1;
    _ng1.export_to_array(connect);
    _ngoper=&_ng1; // for mapping upon exit
  }
//  else if (ret == IMMERSED_1)
//  {
//    //nothing to return
//  }

#ifdef DEBUG_BOOLEAN
  if (ret == EMPTY_X || ret == OK)
  {
    MIO::write("D12.plt", coord, connect, "NGON");
    std::cout << "total number of elements " << connect(0, 0) << std::endl;
  }
  else if (ret == ERROR)
    std::cout << "Diff failed." << std::endl;
  else if (ret == IMMERSED_1)
    std::cout << "Void answer : operand #1 is fully inside operand #2" << std::endl;
#endif

#ifdef FLAG_STEP
  if (connect.cols()) std::cout << "NGON Boolean : nb final PHs : "  << connect[2+connect[1]] << std::endl;
  std::cout << "NGON Boolean : total CPU : " << c.elapsed() << std::endl;
#endif

  return ret;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::Diffsurf
(Coordinate_t& coord, Connectivity_t& connect)
{ 
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif
  
  E_Int err(0);
  
  coord.clear();
  connect.clear();
  
  __init(SURFACE_RIGHT, PRESERVE_RIGHT, DIFF);
 
  if (!_processed)
  {
    // process the transform.
    err = __compute();
    if (err)
      return err;
    
    // color and split : 1, 2 or IN (coloring algo based on "external" info).
    err = __classify_soft();
    if (err)
      return err;
  }
  
  coord=_coord;
  __compact_and_join(_ng1, coord);
  _ng1.export_to_array(connect);
  _ngoper=&_ng1; // for mapping upon exit
 
  
#ifdef DEBUG_BOOLEAN
  MIO::write("D12.plt", coord, connect, "NGON");
#endif
  
#ifdef FLAG_STEP
  if (connect.cols()) std::cout << "NGON Boolean : nb final PHs : "  << connect[2+connect[1]] << std::endl;
  std::cout << "NGON Booelan : total CPU : " << c.elapsed() << std::endl;
#endif

  return err;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::Union
(Coordinate_t& coord, Connectivity_t& connect, eInterPolicy XPol, eMergePolicy MPol)
{
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif

  eRetCode ret(OK);

  coord.clear();
  connect.clear();

  __init(XPol, MPol, UNION);
  
  if (!_processed) // process the transform.
      ret = __compute();
  
  if (ret == OK) // valid computation and common part exists
  {
    switch (XPol)
    {
      case BOTH_SURFACE:
      case SOLID_NONE:
      {
        coord = _coord;
        __compact_and_join(_ngXs, coord);
        _ngXs.export_to_array(connect);
        _ngoper=&_ngXs; // for mapping upon exit

        break;
      }
      case SOLID_RIGHT:
      {
        // color and split : 1, 2 or IN (coloring algo based on "external" info).
        E_Int err = __classify_soft();
        if (err)
          return err;
        
#ifdef MANUAL
        {
          K_FLD::FloatArray coord=_coord;
          __compact_and_join(_ngXh, coord);
          K_FLD::IntArray connect;
          _ngXh.export_to_array(connect);
          MIO::write("subL.plt", coord, connect, "NGON");
        }
#endif
          
#ifdef DEBUG_BOOLEAN 
        std::cout << "nb hards : " << _ngXh.PHs.size() << std::endl;
        std::cout << "nb types : " << _ngXh.PHs._type.size() << std::endl;
        nb_ghost(_ngXh);
#endif

        //remove ghosts
        __remove_ghosts(_ngXh);
        
        _ngXh.append(_ng1); // D12
          
        coord = _coord;
        
        _ngXh.export_to_array(connect);
        _ngoper=&_ngXh; // for mapping upon exit
        
        break;
      }
      default: break;
    }
  }
  else if (ret == EMPTY_X)
  {
    coord.pushBack(_aCoords1.array());
    E_Int nb_pts1 = coord.getSize();
    coord.pushBack(_crd2);

    ngon_type wNG1(_cNGON1), wNG2(_cNGON2);
    wNG2.PGs.shift(nb_pts1);
    wNG1.append(wNG2);
    wNG1.export_to_array(connect);
    
    // fixme : not handle because anc_PH/PG structure is not appropriate :
    // here we should fill the first row of anc_PH/PG for wNG1 and the second for wNG2.
    _ngoper = 0; 
  }

  if (ret != ERROR)
  {
    ngon_type uni(connect);
    ngon_type::clean_connectivity(uni, coord);
    ngon_type::simplify_pgs(uni, coord);
    __compact_and_join(uni, coord);
    uni.export_to_array(connect);

#ifdef DEBUG_BOOLEAN 
#ifndef WIN32
    MIO::write("U.plt", coord, connect, "NGON");
#else
    MIO::write("U.tp", coord, connect, "NGON");
#endif
    std::cout << "total number of elements " << connect(0, 0) << std::endl;
#endif
  }
#ifdef DEBUG_BOOLEAN 
  else
    std::cout << "Union failed." << std::endl;
#endif

#ifdef FLAG_STEP
  E_Int sz = (connect.cols() == 0)  ? 0 : connect[2+connect[1]];
  std::cout << "NGON Boolean : nb final PHs : "  << sz << std::endl;
  std::cout << "NGON Boolean : total CPU : " << c.elapsed() << std::endl;
#endif

  return ret;
}

/*
     * {
      // CONSERVATIVE REMAPPING
      ngon_unit neighbors;
      std::vector<bool> wall(_ngXs.PGs.size(), false);
  
      for (size_t i=0; i < _ngXs.PGs.size(); ++i)
      wall[i]=(_nPG_to_oPG[i] < _nb_pgs1 ) ? true : false;
 
      _ngXs.build_ph_neighborhood(neighbors, wall);
      
      // Color.
      Vector_t<E_Int> colors;
      //fixme : ElementType is not relevant for this algo, should be moved somewhere else
      K_CONNECT::EltAlgo<K_MESH::Triangle>::coloring (neighbors, colors);
      E_Int colmax = *std::max_element(colors.begin(), colors.end())+1;
  
#ifdef DEBUG_BOOLEAN
     {
       Vector_t<bool> flag(colors.size());  
      for (size_t c=0; c < colmax; ++c)
      {
        for (size_t i=0; i < colors.size(); ++i)
          flag[i]=(colors[i]==c);
    
       std::ostringstream o;
       o << "map_" << c;
    
       NGON_DBG::write(o.str().c_str(), ACoordinate_t(coord), _ngXs, &flag);
      }
     }
#endif
    }
     */

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::Modified_Solid
(Coordinate_t& coord, Connectivity_t& connect, eMergePolicy MergePol)
{ 
  E_Int err(0);
  
  __init(SOLID_RIGHT, MergePol, MODIFIED_SOLID);
  
  if (!_processed)
  {
    // process the transform.
    err = __compute();
    if (err) return err;
  }
  
  coord=_coord;
  __compact_and_join(_ngXh, coord);
  _ngXh.export_to_array(connect);
  
#ifdef DEBUG_BOOLEAN
  MIO::write("Solid.plt", coord, connect, "NGON");
#endif
  
  return err;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::volume_coefficients
(E_Int rec_op, std::vector<E_Int>& rec_ids, std::vector<E_Int>& xdon, std::vector<E_Int>& don_ids, std::vector<E_Float>& don_coefs)
{
  // Mapping (rec_op+1)%2 to rec_op : 0->1 or 1->0
  
  E_Int err(0), don_op((rec_op+1)%2);

  rec_ids.clear();
  xdon.clear();
  don_ids.clear();
  don_coefs.clear();

  if (rec_op != 0 && rec_op != 1)
    return 1;
  
  Coordinate_t coord;
  Connectivity_t connect;
  Intersection(coord, connect, SOLID_NONE, PRESERVE_RIGHT);
  
  if (!_ngoper)
    return 1;
  
  ngon_type ng_don = (rec_op==0) ? _cNGON2 : _cNGON1;
  ngon_type ng_rec = (rec_op==0) ? _cNGON1 : _cNGON2;
  
  ng_don.PGs.updateFacets();
  ng_don.PHs.updateFacets();
  ng_rec.PGs.updateFacets();
  ng_rec.PHs.updateFacets();
  
  const K_FLD::FloatArray & crd_don = (rec_op==0) ? _crd2 : _aCoords1.array();
  const K_FLD::FloatArray & crd_rec = (rec_op==0) ? _aCoords1.array() : _crd2;
  
  std::map<E_Int, std::vector<E_Int> > recId_to_bitIds;
  std::map<E_Int, std::vector<E_Int> >::iterator it;
  
  E_Int nb_phs = _ngoper->PHs.size();
  if (nb_phs == 0)  // Empty intersection.
    return 0;

  //std::cout <<  _ngoper->anc_PH << std::endl;
  //
  for (size_t PHi=0; PHi < nb_phs; ++PHi)
  {
    const E_Int& rec_PHi = _ngoper->anc_PH(rec_op, PHi); // receiver ancestor id.
    const E_Int& don_PHi = _ngoper->anc_PH(don_op, PHi); // donnor ancestor id.
    
    if (rec_PHi == E_IDX_NONE || don_PHi == E_IDX_NONE) // only interested in PH in intersection zone.
      continue;
    
    recId_to_bitIds[rec_PHi].push_back(PHi);
  }
  
#ifdef DEBUG_BOOLEAN
//  E_Int PHi=0;
//  E_Int rid = _ngoper->anc_PH(rec_op, PHi);
//  E_Int did = _ngoper->anc_PH(don_op, PHi);
//  
//  NGON_DBG::draw_PH("receiver.tp", crd_rec, ng_rec, rid);
//  NGON_DBG::draw_PH("donnor.tp", crd_don, ng_don, did);
//  
//  for (size_t i=0; i < recId_to_bitIds[rid].size(); ++i)
//  {
//    std::ostringstream o;
//    o << "don_" << recId_to_bitIds[rid][i] << ".tp";
//    NGON_DBG::draw_PH(o.str().c_str(), coord, *_ngoper, recId_to_bitIds[rid][i]);
//  }
#endif
  
  E_Float v_rec;
  
  //init xdon : recPH delimiter for don_ids and don_coefs
  xdon.resize(recId_to_bitIds.size()+1);
  xdon[0]=0;
  E_Int count(1);
  E_Float Gdum[3];
  DELAUNAY::Triangulator dt;
  //
  for (it = recId_to_bitIds.begin(); it != recId_to_bitIds.end(); ++it, ++count)
  {
    const E_Int& rec_PHi = it->first;
    
    E_Int nb_sources = it->second.size();
    
    rec_ids.push_back(rec_PHi);
    xdon[count]=nb_sources + xdon[count-1];
    
    K_MESH::Polyhedron<STAR_SHAPED>::metrics<DELAUNAY::Triangulator>(dt, crd_rec, ng_rec, rec_PHi, v_rec, Gdum); //volume of donnor ancestor to compute the fraction
    
    if (v_rec <= 0.)
    {
#ifdef DEBUG_BOOLEAN
      std::cout << "sliver-like" << std::endl;
      //NGON_DBG::draw_PH("nullvol.tp", crd_rec, ng_rec, recPHi);
#endif
      v_rec=0.;
    }
    else
      v_rec = 1/v_rec;

    E_Float vcumul=0;
    //
    for (size_t i=0; i < nb_sources; ++i)
    {
      const E_Int& PHibit = it->second[i];
      const E_Int& don_PHi = _ngoper->anc_PH(don_op, PHibit);
      
      don_ids.push_back(don_PHi);

      E_Float v;
      K_MESH::Polyhedron<STAR_SHAPED>::metrics<DELAUNAY::Triangulator>(dt, coord, *_ngoper, PHibit, v, Gdum); //volume of the piece of donnor PH
      vcumul += v;
      
      don_coefs.push_back(v*v_rec);
    }
    
#ifdef DEBUG_BOOLEAN
    E_Float err = vcumul * v_rec;
    err = (err - 1.)/ err; //make it relative
    if (::fabs(err) > E_EPSILON)
      std::cout << "erreur relative entre v_rec vs vcumul : " << ::fabs(err) << std::endl;
#endif
    
  }
  
  return err;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::volume_and_centroid_coefficients
(E_Int rec_op, std::vector<E_Int>& rec_ids, std::vector<E_Int>& xdon, std::vector<E_Int>& don_ids, std::vector<E_Float>& don_coefs, K_FLD::FloatArray& piece_centroids, K_FLD::FloatArray& crd, ngon_type& ngoper, std::vector<E_Int>& piece_ids)
{
  // Mapping (rec_op+1)%2 to rec_op : 0->1 or 1->0
  
  DELAUNAY::Triangulator dt;
  
  E_Int err(0), don_op((rec_op+1)%2);

  rec_ids.clear();
  xdon.clear();
  don_ids.clear();
  don_coefs.clear();
  piece_centroids.clear();
  piece_ids.clear();

  if (rec_op != 0 && rec_op != 1)
    return 1;

  Coordinate_t coord;
  Connectivity_t connect;
  Intersection(coord, connect, SOLID_NONE, PRESERVE_RIGHT);

  if (!_ngoper)
    return 1;
  
  ngoper = *_ngoper;
  ngoper.PGs.updateFacets();
  ngoper.PHs.updateFacets();
  crd=coord;

  ngon_type ng_don = (rec_op==0) ? _cNGON2 : _cNGON1;
  ngon_type ng_rec = (rec_op==0) ? _cNGON1 : _cNGON2;
  
  ng_don.PGs.updateFacets();
  ng_don.PHs.updateFacets();
  ng_rec.PGs.updateFacets();
  ng_rec.PHs.updateFacets();

  const K_FLD::FloatArray & crd_don = (rec_op==0) ? _crd2 : _aCoords1.array();
  const K_FLD::FloatArray & crd_rec = (rec_op==0) ? _aCoords1.array() : _crd2;

  std::map<E_Int, std::vector<E_Int> > recId_to_bitIds;
  std::map<E_Int, std::vector<E_Int> >::iterator it;

  E_Int nb_phs = _ngoper->PHs.size();
  if (nb_phs == 0)  // Empty intersection.
    return 0;

  //std::cout <<  _ngoper->anc_PH << std::endl;
  //
  for (size_t PHi=0; PHi < nb_phs; ++PHi)
  {
    const E_Int& rec_PHi = _ngoper->anc_PH(rec_op, PHi); // receiver ancestor id.
    const E_Int& don_PHi = _ngoper->anc_PH(don_op, PHi); // donnor ancestor id.

    if (rec_PHi == E_IDX_NONE || don_PHi == E_IDX_NONE) // only interested in PH in intersection zone. should not happen when same domain
    {
#ifdef DEBUG_BOOLEAN
      std::cout << "MISS MISS MISS : " << PHi << " has been discareded" << std::endl;
      std::ostringstream o;
      o << "don_" << PHi << ".tp";
      //foo_draw_PH(o.str().c_str(), coord, *_ngoper, PHi);
#endif
      //continue;
      return 1;
    }

    recId_to_bitIds[rec_PHi].push_back(PHi);
  }

#ifdef DEBUG_BOOLEAN
  //E_Int PHi=0;
  E_Int rid = 609;//_ngoper->anc_PH(rec_op, PHi);
  //E_Int did = _ngoper->anc_PH(don_op, PHi);
  //NGON_DBG::draw_PH("donnor.tp", crd_don, ng_don, did);
  
  NGON_DBG::draw_PH("receiver.tp", crd_rec, ng_rec, rid);
  
  for (size_t i=0; i < recId_to_bitIds[rid].size(); ++i)
  {
    std::ostringstream o;
    o << "don_" << recId_to_bitIds[rid][i] << ".tp";
    NGON_DBG::draw_PH(o.str().c_str(), coord, *_ngoper, recId_to_bitIds[rid][i]);
  }
#endif

  // Compute the centroids of the donnor pieces
  K_FLD::FloatArray pcentroids;
  //gcc error: expected primary-expression before '>' token
  ngon_type::template centroids<DELAUNAY::Triangulator>(*_ngoper, coord, pcentroids);
  
  /*K_FLD::FloatArray pcentroids2;
  centers_of_mass(*_ngoper, coord, pcentroids2);
  
  for (size_t i=0; i < pcentroids.cols(); ++i)
  {
    std::cout << pcentroids(0,i) << "/" << pcentroids(1,i) << "/" << pcentroids(2,i) << std::endl;
    std::cout << pcentroids2(0,i) << "/" << pcentroids2(1,i) << "/" << pcentroids2(2,i) << std::endl;
    
    E_Float d = ::sqrt(K_FUNC::sqrDistance(pcentroids.col(i), pcentroids2.col(i), 3));
    
  }*/
  
#ifdef DEBUG_BOOLEAN
  K_FLD::FloatArray rec_centroids;
  K_MESH::Polyhedron<STAR_SHAPED>::centroids<DELAUNAY::Triangulator>(ng_rec, crd_rec, rec_centroids);
#endif

  E_Float v_rec, Gdum[3];

  //init xdon : recPH delimiter for don_ids and don_coefs
  xdon.resize(recId_to_bitIds.size()+1);
  xdon[0]=0;
  E_Int count(1);
  //
  for (it = recId_to_bitIds.begin(); it != recId_to_bitIds.end(); ++it, ++count)
  {
    const E_Int& rec_PHi = it->first;
    
    E_Int nb_sources = it->second.size();
    
    rec_ids.push_back(rec_PHi);
    xdon[count]=nb_sources + xdon[count-1];

    K_MESH::Polyhedron<STAR_SHAPED>::metrics<DELAUNAY::Triangulator>(dt, crd_rec, ng_rec, rec_PHi, v_rec, Gdum); //volume of receptor to compute the fraction

    if (v_rec <= 0.)
    {
#ifdef DEBUG_BOOLEAN
      std::cout << "sliver-like" << std::endl;
      //NGON_DBG::draw_PH("nullvol.tp", crd_rec, ng_rec, recPHi);
#endif
      v_rec=0.;
    }
    else
      v_rec = 1./v_rec;

#ifdef DEBUG_BOOLEAN
    E_Float vcumul=0;
    E_Float acuG[3];
    acuG[0]=acuG[1]=acuG[2]=0.;
#endif
    //
    for (size_t i=0; i < nb_sources; ++i)
    {
      const E_Int& PHibit = it->second[i];
      const E_Int& don_PHi = _ngoper->anc_PH(don_op, PHibit);

      don_ids.push_back(don_PHi);

      E_Float v;
      K_MESH::Polyhedron<STAR_SHAPED>::metrics<DELAUNAY::Triangulator>(dt, coord, *_ngoper, PHibit, v, Gdum); //volume of the piece of donnor PH

#ifdef DEBUG_BOOLEAN
      vcumul += v;
      acuG[0] += pcentroids(0,PHibit)*v;
      acuG[1] += pcentroids(1,PHibit)*v;
      acuG[2] += pcentroids(2,PHibit)*v;
#endif

      don_coefs.push_back(v*v_rec);
      piece_centroids.pushBack(pcentroids.col(PHibit), pcentroids.col(PHibit)+3);
      piece_ids.push_back(PHibit);
    }
    
#ifdef DEBUG_BOOLEAN
    E_Float err = vcumul * v_rec;
    err = (err - 1.)/ err; //make it relative
    if (::fabs(err) > E_EPSILON)
      std::cout << "erreur relative entre v_rec vs vcumul : " << ::fabs(err) << std::endl;
    
    acuG[0] /= vcumul;
    acuG[1] /= vcumul;
    acuG[2] /= vcumul;
    
    E_Float d = ::sqrt(K_FUNC::sqrDistance(acuG, rec_centroids.col(rec_PHi), 3));
    if (d > E_EPSILON)
      std::cout << "barycenter calulation is inconsistent" << std::endl;
#endif
    
  }
  
  return err;
}

/// Conservative transfer of a piecewise constant field
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::conservative_transfer_pwc
(E_Int rec_op /*0 or 1*/, const K_FLD::FloatArray& fields_don, K_FLD::FloatArray& fields_rec, E_Int field/*-1 means all*/)
{
  //
  ngon_type ngd = (rec_op==0) ? _cNGON2 : _cNGON1; //donnor ngon
  ngon_type ngr = (rec_op==0) ? _cNGON1 : _cNGON2; //receptor ngon
  const K_FLD::FloatArray & dcrd = (rec_op==0) ? _crd2 : _aCoords1.array(); // donnor coords
  const K_FLD::FloatArray & rcrd = (rec_op==0) ? _aCoords1.array() : _crd2; //receptor coords

  E_Int nb_dphs = ngd.PHs.size();
  if (fields_don.cols() != nb_dphs)
  {
    std::cout << "ERROR : inconsistency between donnor fields/meshes sizes : " << fields_don.cols() << "/" << nb_dphs << std::endl;
    return 1;
  }
  
  ngd.PGs.updateFacets();
  ngd.PHs.updateFacets();
  ngr.PGs.updateFacets();
  ngr.PHs.updateFacets();

#ifdef FLAG_STEP
  std::cout << "volume coefs ..." << std::endl;
#endif
  std::vector<E_Int> rids/*recpetor ids*/, xdon /*delimiter indirection for dids and vcoefs*/, dids;
  std::vector<E_Float> vcoefs;
  E_Int err = volume_coefficients(1, rids, xdon, dids, vcoefs);

  if (err==1)
    return 1;
  
  E_Int nb_sols = 1;
  E_Int nb_fields = fields_don.rows();
  // if out of sacope value, consider all, only one otherwise
  if (field <= -1 || field >= nb_fields)
  {
    field=0;
    nb_sols = nb_fields;
  }
  
  E_Int nb_rphs = ngr.PHs.size();

  //Redimension in case it is required
  fields_rec.resize(nb_sols, nb_rphs, 0.);
 
#ifdef FLAG_STEP
  std::cout << "transfer " << nb_sols << " PWC fields on " << nb_rphs << " elements..." << std::endl;
#endif
  // Now tranfer
  for (E_Int i=0; i < rids.size(); ++i)
  {
    E_Int& reci = rids[i];
    for (E_Int f=0; f < nb_sols; ++f)
    {
      E_Int fld = (field+f)%nb_fields;
      fields_rec(fld, reci)=0.;
      for (E_Int k=xdon[i]; k<xdon[i+1]; ++k)
        fields_rec(fld, reci) += fields_don(fld, dids[k])*vcoefs[k];
    }
  }
  
  return 0;
}

/// Conservative transfer of a piecewise linear field
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::conservative_transfer_pwl
(E_Int rec_op /*0 or 1*/, const K_FLD::FloatArray& fields_don, const K_FLD::FloatArray& grads_don, K_FLD::FloatArray& fields_rec, E_Int field/*-1 means all*/)
{
  //
  ngon_type ngd = (rec_op==0) ? _cNGON2 : _cNGON1; //donnor ngon
  ngon_type ngr = (rec_op==0) ? _cNGON1 : _cNGON2; //receptor ngon
  const K_FLD::FloatArray & dcrd = (rec_op==0) ? _crd2 : _aCoords1.array(); // donnor coords
  const K_FLD::FloatArray & rcrd = (rec_op==0) ? _aCoords1.array() : _crd2; //receptor coords
  
  DELAUNAY::Triangulator dt;

  E_Int nb_dphs = ngd.PHs.size();
  if (fields_don.cols() != nb_dphs)
  {
    std::cout << "ERROR : inconsistency between donnor fields/meshes sizes : " << fields_don.cols() << "/" << nb_dphs << std::endl;
    return 1;
  }
  
  ngd.PGs.updateFacets();
  ngd.PHs.updateFacets();
  ngr.PGs.updateFacets();
  ngr.PHs.updateFacets();

#ifdef FLAG_STEP
  std::cout << "volume coefs ..." << std::endl;
#endif

  std::vector<E_Int> rids/*recpetor ids*/, xdon /*delimiter indirection for dids and vcoefs*/, dids, pids;
  std::vector<E_Float> vcoefs;
  K_FLD::FloatArray piece_centroids, pcrd;
  ngon_type ngoper;
  E_Int err = volume_and_centroid_coefficients(1, rids, xdon, dids, vcoefs, piece_centroids, pcrd, ngoper, pids);
  if (err==1)
    return 1;
  
  K_FLD::FloatArray don_centroids;
  //gcc error: expected primary-expression before '>' token
  ngon_type::template centroids<DELAUNAY::Triangulator>(ngd, dcrd, don_centroids);

  E_Int nb_sols = 1;
  E_Int nb_fields = fields_don.rows();
  // if out of sacope value, consider all, only one otherwise
  if (field <= -1 || field >= nb_fields)
  {
    field=0;
    nb_sols = nb_fields;
  }
  
  E_Int nb_rphs = ngr.PHs.size();

  //Redimension in case it is required
  fields_rec.resize(nb_sols, nb_rphs, 0.);

#ifdef FLAG_STEP
  std::cout << "transfer " << nb_sols << " PWL fields on " << nb_rphs << " elements..." << std::endl;
#endif

#ifdef DEBUG_BOOLEAN
  K_FLD::FloatArray mass_don_by_pieces(nb_sols, fields_don.cols(), 0.);
  std::vector<E_Float> accumulated_piece_vols_for_donnor(fields_don.cols(), 0.);
  std::vector<E_Float> accumulated_piece_vols_for_receptor(fields_rec.cols(), 0.);
  
  std::vector<E_Float> dvols;
  //gcc error: expected primary-expression before '>' token
  ngon_type::template volumes<DELAUNAY::Triangulator>(dcrd, ngd, dvols);
  std::vector<E_Float> rvols;
  //gcc error: expected primary-expression before '>' token
  ngon_type::template volumes<DELAUNAY::Triangulator>(rcrd, ngr, rvols);
  std::vector<E_Float> pvols;
  //gcc error: expected primary-expression before '>' token
  ngon_type::template volumes<DELAUNAY::Triangulator>(pcrd, ngoper, pvols);
  E_Float V1(0.), V2(0.), VI(0.);
  for (size_t i=0; i < dvols.size(); ++i) V1 += dvols[i];
  for (size_t i=0; i < rvols.size(); ++i) V2 += rvols[i];
  for (size_t i=0; i < pvols.size(); ++i) VI += pvols[i];
  
  std::cout << " VOLUME COMPARISON : " << std::endl;
  std::cout << "V1 : " << V1 << std::endl;
  std::cout << "VI : " << VI << std::endl;
  std::cout << "V2 : " << V2 << std::endl;
  std::cout << "err(V1,VI) : " << ::fabs(V1-VI)/V1 << std::endl;
  std::cout << "err(V2,VI) : " << ::fabs(V2-VI)/V2 << std::endl;
  std::cout << "err(V1,V2) : " << ::fabs(V1-V2)/V1 << std::endl;

  K_FLD::FloatArray rec_centroids;
  //gcc error: expected primary-expression before '>' token
  ngon_type::template centroids<DELAUNAY::Triangulator>(ngr, rcrd, rec_centroids);
  E_Float recG[3], acuG[3], dmax(0.);
  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(rcrd);
#endif
  
  E_Float GcGcp[3];
  
  // Now tranfer
  for (E_Int i=0; i < rids.size(); ++i)
  {
    E_Int& reci = rids[i];
    if (reci == 0)
    {
      std::cout << "dealing with " << i << "-th receptor : " << reci << std::endl;
      std::cout << "donnor fieces from " <<  xdon[i] << " to " << xdon[i+1] << std::endl;
    }
    for (E_Int f=0; f < nb_sols; ++f)
    {
      E_Int fld = (field+f)%nb_fields;
      fields_rec(fld, reci) = 0.;
      
#ifdef DEBUG_BOOLEAN
      acuG[0]=acuG[1]=acuG[2]=0.;
#endif
      
      for (E_Int k=xdon[i]; k<xdon[i+1]; ++k)
      {
        const E_Float* Gc = don_centroids.col(dids[k]); // centroid of the donnor cell
        const E_Float* Gcp = piece_centroids.col(k); // centroid of the donnor cell's current piece
        K_FUNC::diff<3>(Gcp, Gc, GcGcp); // vector between the centroid of the donno cell and its current piece
        
        E_Float gradf = K_FUNC::dot<3>(GcGcp, grads_don.col(dids[k]));
        
        //std::cout << "grads : " << grads_don(0, dids[k]) << "/" << grads_don(1, dids[k]) << "/" << grads_don(2, dids[k]) << std::endl;
        fields_rec(fld, reci) += (fields_don(fld, dids[k]) + gradf) * vcoefs[k];
        //std::cout << "gradf/vcoeef : " << gradf << "/" << vcoefs[k] << std::endl;
         
#ifdef DEBUG_BOOLEAN
        E_Float vv, GG[3];
        K_MESH::Polyhedron<STAR_SHAPED>::metrics<DELAUNAY::Triangulator>(dt, pcrd, ngoper, pids[k], vv, GG);
        E_Float vvv= vcoefs[k] * rvols[reci];
        bool vok = (::fabs(vv - vvv) < E_EPSILON);
        bool xok = (GG[0] == Gcp[0]);
        bool yok = (GG[1] == Gcp[1]);
        bool zok = (GG[2] == Gcp[2]);
        
        if (!vok || !xok || !yok || !zok)
        {
          std::cout << "error" << std::endl;
        }
        
//        if (reci == 0)
//        {
//          std::cout << k << "-th donnor is : " << dids[k] << " with vol coeff & masss : " << vcoefs[k] << " and " << fields_don(fld, dids[k]) << std::endl;
//          std::cout << "its piece volume is : " << vv << " over donnor vol : " << dvols[dids[k]] << std::endl;
//          std::cout << "the receptor mass is : " << rvols[reci] << std::endl;
//          std::cout << "so we should have : " << vcoefs[k] << "equal " << vv/rvols[reci] << std::endl;
//        }
        
        E_Float vpiece = vcoefs[k] * rvols[reci];
        
        mass_don_by_pieces(fld, dids[k]) += (fields_don(fld, dids[k]) + gradf) * vpiece;
        
        accumulated_piece_vols_for_donnor[dids[k]] += vpiece;
        accumulated_piece_vols_for_receptor[reci] += vpiece;
        
        acuG[0] += Gcp[0]*vcoefs[k];
        acuG[1] += Gcp[1]*vcoefs[k];
        acuG[2] += Gcp[2]*vcoefs[k];
#endif
      }
      
#ifdef STAR_SHAPE_WAY
      fields_rec(fld, reci) /= rvols[reci];
#endif

    }
#ifdef DEBUG_BOOLEAN
    E_Float recGx = rec_centroids(0, reci);
    E_Float recGy = rec_centroids(1, reci);
    E_Float recGz = rec_centroids(2, reci);
    E_Float d = ::sqrt(K_FUNC::sqrDistance(acuG, rec_centroids.col(reci), 3));
    dmax = (d > dmax) ? d: dmax;
    assert (d < E_EPSILON);
#endif
  }
  
#ifdef DEBUG_BOOLEAN
  //std::cout << "DMAX is " << dmax << std::endl;
  for (size_t i=0; i < fields_rec.cols(); ++i)
  {
    E_Float vol_rec = rvols[i];
    E_Float volcumul = accumulated_piece_vols_for_receptor[i];
    
    bool error_vol  = (::fabs(vol_rec-volcumul)/vol_rec > E_EPSILON);
    
    if (error_vol)
    {
      std::cout << "TRANSFER ERROR for receptor " << i << std::endl;
      std::cout << "volumes (rec vs acc) : " << vol_rec << "/" << volcumul << std::endl;
      return 1;
    }
  }
  
  for (size_t i=0; i < fields_don.cols(); ++i)
  {
    E_Float val_don   = fields_don(0, i);
    E_Float vol_don = dvols[i];
    E_Float mass_don = val_don*vol_don;
    
    E_Float masscumul = mass_don_by_pieces(0, i);
    E_Float volcumul    = accumulated_piece_vols_for_donnor[i];
        
    
    bool error_vol  = (::fabs(vol_don-volcumul)/vol_don > E_EPSILON);
    bool error_mass = (::fabs(mass_don-masscumul)/mass_don > E_EPSILON);
    
    if (error_vol || error_mass)
    {
      std::cout << "TRANSFER ERROR for donnor " << i << std::endl;
      std::cout << "volumes (don vs acc) : " << vol_don << "/" << volcumul << std::endl;
      std::cout << "relative volume error : " << ::fabs(vol_don-volcumul)/vol_don << std::endl;
      std::cout << "masses (don vs acc) : " << mass_don << "/" << masscumul << std::endl;
      std::cout << "relative mass error : " << ::fabs(mass_don-masscumul)/mass_don << std::endl << std::endl;
      return 1;
    }
  }
#endif
  
  return 0;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::total_mass(const K_FLD::FloatArray& crd, const K_FLD::IntArray& ng_cnt, const K_FLD::FloatArray& fields, const std::vector<E_Float>& vols, std::vector<E_Float>& masses)
{
  
  masses.clear();
  E_Int nb_fields(fields.rows()), nb_elts(fields.cols());
  masses.resize(nb_fields, 0.);
  
  for (size_t i=0; i < nb_elts; ++i)
  {
    for (size_t j=0; j < nb_fields; ++j)
      masses[j] += vols[i]*fields(j, i);
  }

  return 0;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::XcellN (std::vector<E_Float>& xcellN)
{
  // WARNING : second operand (_cNGON2) hides the first one => so we compute only cellN on the first one
  
  xcellN.clear();
  
  const ngon_type ng1 = _cNGON1;
  const ngon_type ng2 = _cNGON2;
  const K_FLD::FloatArray & crd1 = _aCoords1.array();
  //const K_FLD::FloatArray & crd2 = _aCoords2.array();

  //
  Coordinate_t coord;
  Connectivity_t connect;
  E_Int ret = Diff(coord, connect, SOLID_RIGHT, PRESERVE_RIGHT);
  
#ifdef DEBUG_BOOLEAN
  switch (ret)
  {
    case OK : std::cout << "diff returned OK" << std::endl; break;
    case EMPTY_X : std::cout << "diff returned EMPTY_X" << std::endl; break;
    case ERROR : std::cout << "diff returned ERROR" << std::endl; break;
    case IMMERSED_1 : std::cout << "diff returned IMMERSED_1" << std::endl; break;
  }

  if (connect.cols())
  {
    ngon_t<K_FLD::IntArray> ngo(connect);  
    K_FLD::IntArray cnto;
    ngo.export_to_array(cnto);
    MIO::write("diff.plt", coord, cnto, "NGON");
  }
#endif
  
  if ( (ret == EMPTY_X) || (ret == IMMERSED_1) )
  {
    E_Float VAL = (ret == EMPTY_X) ? VISIBLE : BLANKED; // in case of full immersion all the initial cells are blanked
    xcellN.resize(ng1.PHs.size(), VAL); //the prioritized mesh is outside the bgm
    return 0;
  }
  
  // initialize to all masked because
  // we use DIFF (rather than UNION for hpc)
  // so some element in 1 are missing, meaning that they are hidden.
  // they won't be in ngoper so we set everything to 0.
  xcellN.resize(ng1.PHs.size(), BLANKED); 
  
  
  if (!_ngoper)
  {
    std::cout << "Supermesh is missing" << std::endl;
    return 1;
  }

#ifdef DEBUG_BOOLEAN
//  {
//  E_Int anc0 = ng1.anc_PH(0,0);
//  std::cout << "ancestor de 0 : " <<  anc0 << std::endl;
//  NGON_DBG::draw_PH("focus0.plt", crd1, ng1, 0);
//  NGON_DBG::draw_PH("ancestorde0.plt", _aCoords1.array(), wNG1initial, anc0);
//  }
  K_FLD::FloatArray crdG=coord;
  E_Int shft=crdG.cols();
  crdG.pushBack(crd1);
#endif
  
  ng1.PGs.updateFacets();
  ng1.PHs.updateFacets();
  _ngoper->PGs.updateFacets();
  _ngoper->PHs.updateFacets();

  E_Int nb_phs = _ngoper->PHs.size();
  if (nb_phs == 0)  // Empty intersection.
    return 0;
  
  if (_ngoper->PHs._ancEs.cols() != nb_phs)
  {
    std::cout << "PH history is inconsistent" << std::endl;
    return 1;
  }
  
  DELAUNAY::Triangulator dt;

  std::vector<bool> keeinope(nb_phs, false);
  //
  for (E_Int PHi=0; PHi < nb_phs; ++PHi)
  {
    const E_Int& ancPH1i = _ngoper->PHs._ancEs(0, PHi);     
    assert (ancPH1i != E_IDX_NONE);
    
    bool unused = (_ngoper->PHs._type[PHi] == UNUSED);
    
    if (unused)
    {
      xcellN[ancPH1i] = VISIBLE;
      continue;
    }
    keeinope[PHi]=true;
    
    // the piece has 2 ancestors so it's a cutted cell
    // so compute the ratio
    E_Float vold, Gdum[3];
    K_MESH::Polyhedron<STAR_SHAPED>::metrics<DELAUNAY::Triangulator>(dt, crd1, ng1.PGs, ng1.PHs.get_facets_ptr(ancPH1i), ng1.PHs.stride(ancPH1i), vold, Gdum);
    E_Float vnew;
    K_MESH::Polyhedron<STAR_SHAPED>::metrics<DELAUNAY::Triangulator>(dt, coord, _ngoper->PGs, _ngoper->PHs.get_facets_ptr(PHi), _ngoper->PHs.stride(PHi), vnew, Gdum);
    
    xcellN[ancPH1i] = vnew/vold;
    if (xcellN[ancPH1i] > 1. - E_EPSILON) xcellN[ancPH1i] = VISIBLE; // cutoff to 1.

#ifdef DEBUG_BOOLEAN
    assert(xcellN[ancPH1i] > BLANKED);
    if (xcellN[ancPH1i] > BLANKED && xcellN[ancPH1i] < VISIBLE)
    {
      std::cout << "one partial value : " << vnew/vold << std::endl;
      ngon_type ng;
      ng.addPH(ng1, ancPH1i);
      ng.PGs.shift(shft);
      ng.addPH(*_ngoper, PHi);
      
      K_FLD::IntArray cnto;
      ng.export_to_array(cnto);
      
      //std::ostringstream o;
      //o << "fraction_" << PHi << ".plt";
      //MIO::write(o.str().c_str(), crdG, cnto, "NGON");
    }
#endif
  }
  
#ifdef DEBUG_BOOLEAN
// Fully hidden
  size_t sz = xcellN.size();
  std::vector<bool> keep(xcellN.size(), false);
  for (size_t i=0; i < xcellN.size(); ++i)
  {
    keep[i] = (0. == xcellN[i]);
  }
  assert(ng1.PHs.size() == xcellN.size());
  ngon_t<K_FLD::IntArray> bngo;
  Vector_t<E_Int> npgids;
  ng1.select_phs(ng1, keep, npgids, bngo);
  
  ngon_t<K_FLD::IntArray> ngo;
  ngo.PHs=bngo.PHs;
  ngo.PGs = bngo.PGs;
  
  K_FLD::IntArray cnto;
  ngo.export_to_array(cnto);
  MIO::write("fully_hidden.plt", crd1, cnto, "NGON");
#endif

  return 0;
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__init(eInterPolicy XPol, eMergePolicy MPol, eOperation Op)
{
  if ( (_XPol != XPol) || (_MPol != MPol)) _processed = false;
 
  _XPol=XPol;
  _MPol=MPol;
  _Op = Op;
  
  _build_hard=((_XPol == SOLID_RIGHT) && (_Op == UNION)) || (_Op == MODIFIED_SOLID);
  
  if (!_processed) // clear everything
  {
    _ngXs.clear();
    _ngXh.clear();
    _ng1.clear();
    _ng2.clear();
    _ngoper = NULL;
    _nb_cells2 = 0;
    
    _coord.clear();
    _zones.clear();

    _anc_PH_for_PHT3s.clear();
    _anc_PG.clear();

    _normals.clear();
    
    _palmares.clear();
    _tmp_set_int.clear();
    _tmp_vec.clear();
    _tmp_deq.clear();
    _tmp_IntArray2.clear();
    _tmp_FloatArray.clear();
    _tmp_set_oedge.clear();   
    _tmp_map_nn.clear();
  }
}

///
TEMPLATE_COORD_CONNECT
typename NGON_BOOLEAN_CLASS::eRetCode
NGON_BOOLEAN_CLASS::__get_working_PGs
(eInterPolicy XPol, eMergePolicy MPol, ngon_unit& wPGs, E_Int& nb_pgs1, ngon_unit& extrawPGs)
{
  ngon_type wNG1(_cNGON1), wNG2(_cNGON2);

#ifdef FLAG_STEP
  chrono c, c0;
  c.start();
  c0.start();
#endif

#ifdef DEBUG_BOOLEAN
  //fixme : should be done before the boolean call CODE_1
  E_Int nb_modifs = ngon_type::clean_connectivity(wNG1, _aCoords1.array());
  nb_modifs += ngon_type::clean_connectivity(wNG2, _crd2);
  
  if (nb_modifs)
    std::cout << "WARNING : INPUTS HAVE NOT BEEN CLEANED BEFORE CALLING THE BOOLEAN" << std::endl << std::endl;
#endif
  
  // Initialize history structures
  wNG1.PGs._ancEs.resize(2, 1);
  wNG1.PHs._ancEs.resize(2, 1);
  wNG2.PGs._ancEs.resize(2, 1);
  wNG2.PHs._ancEs.resize(2, 1);
  K_CONNECT::IdTool::init_inc(wNG1.PGs._ancEs, 0, wNG1.PGs.size());
  K_CONNECT::IdTool::init_inc(wNG1.PHs._ancEs, 0, wNG1.PHs.size());
  K_CONNECT::IdTool::init_inc(wNG2.PGs._ancEs, 1, wNG2.PGs.size());
  K_CONNECT::IdTool::init_inc(wNG2.PHs._ancEs, 1, wNG2.PHs.size());

#ifdef DEBUG_BOOLEAN
  assert (wNG1.is_consistent());
  assert (wNG2.is_consistent());

  std::cout << "GHOST creation : " << std::endl;
  nb_ghost(wNG2);
#endif
  
#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : init NGON data struct : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  // fixme : is designed for classic octrees i.e. with manifold PG contours.
  // it should be launched only in that case. closeness is guaranteed starting from basic elements meshes.
  // do we need a close for other topologies ? dealing with non manifold edges ?
  // it has been moved at the beginning (i.e. applied to the whole input) to prevent issues when attaching the modified bits with untouched ones (union)
  // issues are related to non conformity (some PGs supposed to be the same, sharing 2 PHs but one is refined).
  if (wNG1.PHs.size() != wNG1.PGs.size())//if it's not a surface : partial easy filter before the right one based on conform-ngon-or-not
    ngon_type::close_phs(wNG1, _aCoords1.array());
  if (wNG2.PHs.size() != wNG2.PGs.size())//if it's not a surface : partial easy filter before the right one based on conform-ngon-or-not
    ngon_type::close_phs(wNG2, _crd2);
  
#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : close_phs : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  //Externality : MUST BE PRECEEDED BY close_phs IN CASE OF OCTREES
  wNG1.flag_externals(INITIAL_SKIN);
  wNG2.flag_externals(INITIAL_SKIN);
  eRetCode err = __reorient_externals(XPol, wNG1, wNG2);
  if (err)
    return err;

  // Add ghost cells for prioritized overlapping PGs
  // fixme : HACK _pglist2 is used both for body walls (BdnWall) and here to extrude ghost cells on BndUser
  // MUST BE PRECEEDED BY REORIENT EXTERNALS (to create the ghost layer consistently)
  bool has_ghosts = false;
  if (XPol==SOLID_RIGHT && !_pglist2[1].empty())
  {
    //std::cout << "ghost list has been passed" << std::endl;
    // GHOST ARE THOSE IN wNG WITH IDS ABOVE _nb_cells2
    _nb_cells2 = wNG2.PHs.size(); // valued before appending ghost
    E_Int err = ngon_type::append_with_ghost_cells(_crd2, wNG2, _pglist2[1]);
    has_ghosts = true;
    wNG1.flag_externals(INITIAL_SKIN); // do it again to reflect ghost presence
    wNG2.flag_externals(INITIAL_SKIN);
    eRetCode er = __reorient_externals(XPol, wNG1, wNG2);
    if (er)
      return er;
  }
  
#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : externality : " << c.elapsed() << std::endl;
  c.start();
#endif

  if (_XPol != BOTH_SURFACE)
  {
    ngon_unit neighbors;
    wNG1.build_ph_neighborhood(neighbors);
    wNG1.build_F2E(neighbors, _F2E);
  }

  K_FLD::IntArray F2E2;
  if (_XPol == SOLID_NONE || _XPol == SOLID_RIGHT/*hack*/)
  {
    ngon_unit neighbors;
    wNG2.build_ph_neighborhood(neighbors);
    wNG2.build_F2E(neighbors, F2E2);
  }
  
#ifdef DEBUG_BOOLEAN
  std::cout << "GHOST apres F2E : " << std::endl;
  nb_ghost(wNG2);
#endif
  
#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : orienting (build_F2E) : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  if (XPol==SOLID_RIGHT) // discard inner holes. Is there an application requiring them ? fixme : something has to be done for BOTH_SURFACE and SURFACE_RIGHT modes.
  {
    //std::cout << "wall pglist state : " << _pglist2[0].size() << std::endl;
    if (!_pglist2[0].empty())
      __discard_prescribed_polygons(_crd2, wNG2, _pglist2[0]);
    else
      __discard_holes_by_box(_crd2, wNG2);
  }
  
#ifdef DEBUG_BOOLEAN
  std::cout << "GHOST apres discard holes : " << std::endl;
  nb_ghost(wNG2);
#endif
  
#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : __discard_holes : " << c.elapsed() << std::endl;
  c.start();
#endif

#ifdef FLAG_STEP
  std::cout << std::endl;
  std::cout << "__get_working_PGs : PREPARE INPUTS : " << c0.elapsed() << std::endl << std::endl;
#endif

  E_Int nb_ph10 = wNG1.PHs.size();
  E_Int nb_ph20 = wNG2.PHs.size();
  E_Int nb_pg10 = wNG1.PGs.size();
  E_Int nb_pg20 = wNG2.PGs.size();
    
#ifdef FLAG_STEP
  c.start();
#endif
  
#ifdef DEBUG_BOOLEAN
  std::cout << "avant focus" << std::endl;
  nb_ghost(wNG2);
#endif

if (_XPol != BOTH_SURFACE)
{
  eRetCode err = __focus_on_intersection_zone(XPol, wNG1, wNG2, _ng1, _ng2);
  if (err)
    return err;
  
#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : __focus_on_intersection_zone : " << c.elapsed() << std::endl;
  c.start();
#endif
  
#ifdef DEBUG_BOOLEAN
  std::cout << "apres focus wNG2" << std::endl;
  nb_ghost(wNG2);
  std::cout << "apres focus _ng2" << std::endl;
  nb_ghost(_ng2);
#endif
  
  
  // to accelerate usage of mapping : good to know if a cell has not been involved (e.g. cellN =1 straitforwardly)
  if (!has_ghosts) // fixme : better way to catch that we are in assembly mode
  {
    _ng1.PHs._type.clear();
    _ng1.PHs._type.resize(_ng1.PHs.size(), UNUSED);
    _ng2.PHs._type.clear();
    _ng2.PHs._type.resize(_ng2.PHs.size(), UNUSED);
  }

  // reduce F2E to the working sets
  if (nb_pg10 != wNG1.PGs.size())
  {
    E_Int nb_pgs = wNG1.PGs.size();
    K_FLD::IntArray f2e(2, nb_pgs);
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      E_Int ancPGi = wNG1.PGs._ancEs(0, i);
      f2e(0, i) = _F2E(0, ancPGi);
      f2e(1, i) = _F2E(1, ancPGi);
    }

    _F2E = f2e;
  }
  if (nb_pg20 != wNG2.PGs.size() && (F2E2.cols() != 0))
  {
    E_Int nb_pgs = wNG2.PGs.size();
    K_FLD::IntArray f2e(2, nb_pgs, E_IDX_NONE);
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      E_Int ancPGi = wNG2.PGs._ancEs(1, i);
      if (ancPGi == E_IDX_NONE) continue;
      f2e(0, i) = F2E2(0, ancPGi);
      f2e(1, i) = F2E2(1, ancPGi);
    }

    F2E2 = f2e;
  }
  
#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : update orienting after focus : " << c.elapsed() << std::endl;
#endif
}

//#ifdef DEBUG_BOOLEAN
//{
//  char* fname1 = "working1";
//  NGON_DBG::write(fname1, _aCoords1, wNG1);
//  char* fname2 = "left1";
//  NGON_DBG::write(fname2, _aCoords1, _ng1);
//  fname1 = "working2";
//  NGON_DBG::write(fname1, _aCoords2, wNG2);
//  fname2 = "left2";
//  NGON_DBG::write(fname2, _aCoords2, _ng2);
// 
//  NGON_DBG::draw_PGT3s(_aCoords1.array(), wNG1.PGs);
//  NGON_DBG::draw_PGT3s(_aCoords2.array(), wNG2.PGs);
//}
//#endif
  
#ifdef FLAG_STEP
  c.start();
#endif

  // Flag CONNEXION_SKIN type
  if ((XPol != BOTH_SURFACE) && (nb_ph10 != wNG1.PHs.size()))
  {
    // Set externality on the focused set
    ngon_type tmp(wNG1);
    tmp.flag_externals(INITIAL_SKIN);
  
    // Now compare this flag with the one on the complete config
    for (E_Int i = 0; i < wNG1.PGs.size(); ++i)
    {
      if ((tmp.PGs._type[i] == INITIAL_SKIN) && (wNG1.PGs._type[i] == INNER))
        wNG1.PGs._type[i] = CONNEXION_SKIN;
    }
    
    wNG1.flag_external_phs(CONNEXION_SKIN); 
  } 
  
  if ((XPol != SURFACE_RIGHT) && (XPol != BOTH_SURFACE) && (nb_ph20 != wNG2.PHs.size())) // i.e. if the second operand is not just a surface
  {
   ngon_type tmp(wNG2);
   tmp.flag_externals(INITIAL_SKIN);
   //PGs
   for (E_Int i=0; i < wNG2.PGs.size(); ++i)
   {
     if ((tmp.PGs._type[i] == INITIAL_SKIN) && (wNG2.PGs._type[i] == INNER))
       wNG2.PGs._type[i] = CONNEXION_SKIN;
   }
   
   wNG2.flag_external_phs(CONNEXION_SKIN); 
  }

#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : flag new skin type : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  // ONE COORDINATE CONTAINER (so reindex globally).
#ifdef FLDARR
  K_FLD::FloatArray crd1(coord1);
  K_FLD::FloatArray crd2(coord2);
#else
  const K_FLD::FloatArray& crd1 = _aCoords1.array();
  const K_FLD::FloatArray& crd2 = _crd2;
#endif
  
  _coord.clear();
  _coord.pushBack(crd1);
  E_Int nb_pts1= _coord.getSize();
  _coord.pushBack(crd2);
  
  //_crd2.release();not possible because currently if EMPTY_X, _crd2 is still required
    
  if (wNG2.PGs.size())
    wNG2.PGs.shift(nb_pts1);
  if (_ng2.PGs.size())
    _ng2.PGs.shift(nb_pts1);

#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : one container : " << c.elapsed() << std::endl;
#endif
  
#ifdef DEBUG_BOOLEAN
//{
//  K_FLD::ArrayAccessor<Coordinate_t> ac(_coord);
//  char* fname1 = "working1g";
//  NGON_DBG::write(fname1, ac, wNG1);
//  char* fname2 = "left1g";
//  NGON_DBG::write(fname2, ac, _ng1);
//  fname1 = "working2g";
//  NGON_DBG::write(fname1, ac, wNG2);
//  fname2 = "left2g";
//  NGON_DBG::write(fname2, ac, _ng2);
//
//  //NGON_DBG::draw_PGT3s(_coord, wNG1.PGs);
//  //NGON_DBG::draw_PGT3s(_coord, wNG2.PGs);
//}
#endif
  
  // Conformize PGs inside PHs
  
#ifdef DEBUG_BOOLEAN    
  //NGON_DBG::draw_PH("PHi.plt", _coord,*wng_first, 0);
#endif

#ifdef DEBUG_BOOLEAN
//  {
//  K_FLD::ArrayAccessor<Coordinate_t> ac(_coord);
//  char* fname1 = "working1h";
//  NGON_DBG::write(fname1, ac, wNG1);
//  char* fname2 = "left1h";
//  NGON_DBG::write(fname2, ac, _ng1);
//  fname1 = "working2h";
//  NGON_DBG::write(fname1, ac, wNG2);
//  fname2 = "left2h";
//  NGON_DBG::write(fname2, ac, _ng2);
//
//  //NGON_DBG::draw_PGT3s(_coord, wNG1.PGs);
//  //NGON_DBG::draw_PGT3s(_coord, wNG2.PGs);
//  }
#endif
  
#ifdef FLAG_STEP
  c.start();
#endif
    
  // WORKING PGs  
  ngon_unit *wpgS(0), tmp;
  Vector_t<E_Int> oPGids;

  K_FLD::IntArray *F2ES(0), F2Etmp;

  if (XPol == SOLID_RIGHT)
  {
    //separate skin and open layer PGs
    wNG2.PGs.extract_of_type(INITIAL_SKIN, tmp, oPGids);//get E32 skin faces
    //
    wpgS=&tmp;
    
    // skin
    F2Etmp.resize(2, wpgS->size());
    for (E_Int i = 0; i < F2Etmp.cols(); ++i)
    {
      F2Etmp(0,i) = F2E2(0, oPGids[i]);
      F2Etmp(1,i) = F2E2(1, oPGids[i]);
    }
    F2ES=&F2Etmp;
   
    //open layer
    std::sort(oPGids.begin(), oPGids.end());
    Vector_t<E_Int> pgnids, phnids;
    wNG2.remove_pgs(oPGids, pgnids, phnids); 
    
    //open layer
    extrawPGs = wNG2.PGs;   
    extrawPGs._type =  wNG2.PGs._type;
    extrawPGs._ancEs = wNG2.PGs._ancEs;

    _extraF2E.resize(2, extrawPGs.size());
    for (E_Int i = 0; i < F2E2.cols(); ++i)
    {
      if (pgnids[i] == E_IDX_NONE)
        continue;
      _extraF2E(0, pgnids[i])=F2E2(0, i);
      _extraF2E(1, pgnids[i])=F2E2(1, i);
    }
  }
  else if (XPol == SOLID_NONE)//get all E31 ans E32 faces
  {
    wpgS = &wNG2.PGs;
    F2ES = &F2E2;
  }
  
  _nb_pgs1 = nb_pgs1 = wNG1.PGs.size();
  
  wPGs = wNG1.PGs;
  wPGs.append(*wpgS);
  _anc_PG = wPGs._ancEs;
  _F2E.pushBack(*F2ES);
  assert(_F2E.cols() == wPGs.size());

#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__get_working_PGs : form wPGs set : " << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_BOOLEAN
{
  char* fname1 = "working";
  K_FLD::ArrayAccessor<Coordinate_t> ac(_coord);
  NGON_DBG::write(fname1, ac, wPGs);
  char* fname11 = "extra";
  NGON_DBG::write(fname11, ac, extrawPGs);
  char* fname2 = "left1";
  NGON_DBG::write(fname2, ac, _ng1);
  fname2 = "left2";
  NGON_DBG::write(fname2, ac, _ng2);
}
#endif
  
  return OK;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__conformize
(K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, Vector_t<E_Int>& colors, E_Int X0, Vector_t<E_Int>& priority)
{
    
#ifdef DEBUG_BOOLEAN
    //bool ok = NUGA::TRI_BooleanOperator::isClosed(coord, connectT3);
#endif

  K_FLD::FloatArray crd(coord);
  TRI_Conformizer<3> conformizer(true/* keep track of nodes history*/);
  E_Int err = conformizer.run(crd, connectT3, colors, &priority, _tolerance, X0, 1 /*one iter only*/);
  if (err)
    return err;
  
#ifdef DEBUG_BOOLEAN
    //ok = NUGA::TRI_BooleanOperator::isClosed(coord, connectT3);
#endif
  
  _nodes_history = conformizer.get_node_history();
  assert(_nodes_history.size() == coord.cols());
  
  //E_Int mid = K_CONNECT::IdTool::max(_nodes_history);
  //assert (mid < crd.cols());  
  
  E_Int shft = coord.getSize();
  
  //update node history
  for (size_t i = 0; i < _nodes_history.size(); ++i) 
     _nodes_history[i] = (_nodes_history[i] != E_IDX_NONE) ? _nodes_history[i] + shft : i;// Put back in history nodes not involved in conformizing
    
  connectT3.shift(shft);
  coord.pushBack(crd);
  
#ifdef DEBUG_BOOLEAN
  {
    //const char* fname = "conformized.mesh";
    //MIO::write(fname, coord, connectT3, "TRI");
  }
#endif
  
  return 0;
}

///
TEMPLATE_COORD_CONNECT
typename NGON_BOOLEAN_CLASS::eRetCode
NGON_BOOLEAN_CLASS::__compute()
{
  ngon_unit wPGs, extrawPGs;
  E_Int nb_pgs1;
  K_FLD::IntArray connectT3;
  Vector_t<E_Int> priority;
  
  //
  eRetCode er = __process_intersections(wPGs, nb_pgs1, extrawPGs, connectT3, _nT3_to_oPG, priority);
  if (er) return er;
  
  //
  if (_Op != MODIFIED_SOLID)
  {
    
#ifdef FLAG_STEP
    std::cout << "NGON Boolean : prepare is_skin"  << std::endl;
#endif
    Vector_t<E_Int> is_skin(connectT3.cols(), 0);
    for (E_Int i=0; i < connectT3.cols(); ++i)
    {
      const E_Int PGi = _nT3_to_oPG[i];
      is_skin[i]=wPGs._type[PGi];
      if (is_skin[i] == INITIAL_SKIN && PGi >= nb_pgs1)
        is_skin[i]=2;
    }

    for (size_t i = 0; i < priority.size(); ++i)
    {
      const E_Int Ti = priority[i];
      is_skin[Ti] = -is_skin[Ti];// required to distinguish duplicated skin when classifying and removing parasites
    }

#ifdef DEBUG_BOOLEAN
    {
      Vector_t<E_Int> colors=is_skin;
      if (!priority.empty())
      {
        for (E_Int i = 0; i < priority.size(); ++i)
          colors[priority[i]] = 4+::abs(is_skin[priority[i]]);
      }
      MIO::write("all_skins.mesh", _coord, connectT3, "TRI", 0, &colors);

      K_FLD::IntArray cT3;
      colors.clear();
      for (E_Int i = 0; i < is_skin.size(); ++i)
        if (is_skin[i] == 1 || is_skin[i] == 2)
        {
          cT3.pushBack(connectT3.col(i), connectT3.col(i) + 3);
          colors.push_back(is_skin[i]);
        }
      MIO::write("init_skins.mesh", _coord, cT3, "TRI", 0, &colors);
    }    
#endif

    // Duplicate the T3s to have both orientations
    K_FLD::IntArray connectT3o(connectT3);
    __duplicate_orient(connectT3o);

#ifdef FLAG_STEP
    std::cout << "NGON Boolean : duplicate orient done"  << std::endl;
#endif

    // Build the PHT3s by constructing the right neighborhood for each oriented T3.
    std::map<E_Int, Vector_t<E_Int> > PHT3s;
    E_Int err = __build_PHT3s(_coord, connectT3, connectT3o, is_skin, PHT3s);
    
#ifdef DEBUG_BOOLEAN
    {
      std::vector<E_Int> colors(PHT3s.size(), 0);
      NGON_DBG::draw_PHT3s("PHT3s.mesh", _coord, connectT3o, PHT3s, colors);
    }
#endif
    
    if (err)
    {
#ifdef FLAG_STEP
    std::cout << "NGON Boolean : ERROR in  : __build_PHT3s"  << std::endl;
#endif
      return ERROR;
    }
#ifdef FLAG_STEP
    std::cout << "NGON Boolean : nb PHT3s constructed : "  << PHT3s.size() << std::endl;
#endif

    E_Int shift = connectT3.cols();

    __classify_skin_PHT3s(PHT3s, is_skin, shift
#ifdef DEBUG_BOOLEAN
, connectT3o
#endif
            );// zoning for soft part only

    // remove duplicated T3s and associated normals.
    __remove_orientations(PHT3s, shift);

    // PHT3 -> PH (with aggregation).
    err= __build_PHs(PHT3s, connectT3, is_skin, _nT3_to_oPG, _ngXs);
    if (err)
      return ERROR;

#ifdef FLAG_STEP
    std::cout << "NGON Boolean : nb PHs constructed : "  << _ngXs.PHs.size() << std::endl;
#endif
  }

#ifdef DEBUG_BOOLEAN
  NGON_DBG::write("ngXs", ACoordinate_t(_coord), _ngXs);
#endif

  if (_build_hard)
  {
    K_FLD::IntArray connectHard;
    Vector_t<E_Int> is_skin;
    E_Int err = __build_connect_hard(_coord, extrawPGs, nb_pgs1, connectT3, connectHard, _nT3_to_oPG, is_skin);
    if (err)
      return ERROR;
    
#ifdef DEBUG_BOOLEAN
    //is_skin[0] = 1;
    MIO::write("connectHard.mesh", _coord, connectHard, "TRI", 0, &is_skin);
#endif
    
    // Duplicate the T3s to have both orientations
    K_FLD::IntArray connectT3o(connectHard);
    __duplicate_orient(connectT3o);
        
    // Build the PHT3s by constructing the right neighborhood for each oriented T3.
    std::map<E_Int, Vector_t<E_Int> > PHT3s;
    err = __build_PHT3s(_coord, connectHard, connectT3o, is_skin, PHT3s);
    if (err) return ERROR;

    E_Int shift = connectHard.cols();
    __set_PH_history(PHT3s, is_skin, shift, -1, _extraF2E, _anc_PH_for_PHT3s, false /* is not soft*/);

    // remove duplicated T3s and associated normals.
    __remove_orientations(PHT3s, connectHard.cols());
      
    // PHT3 -> PH (with aggregation).
    err = __build_PHs(PHT3s, connectHard, is_skin, _nT3_to_oPG, _ngXh);
    if (err)
      return ERROR;
    
#ifdef DEBUG_BOOLEAN
    std::cout << "_ngXh just built" << std::endl;
    nb_ghost(_ngXh);
#endif
    
#ifdef DEBUG_BOOLEAN
  NGON_DBG::write("ngXh", ACoordinate_t(_coord), _ngXh);
#endif
  }
  
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif
  
  // Aggregate the soft bits
#ifdef MANUAL
  {
    K_FLD::FloatArray coord=_coord;
    __compact_and_join(_ng1, coord);
    K_FLD::IntArray connect;
    _ng1.export_to_array(connect);
    MIO::write("ocext.plt", coord, connect, "NGON");
  }
#else
  _ngXs.append(_ng1);//tocheck : append to ngX to avoid to have to modify _zones. but might be more efficient to append the smaller to the bigger..
#endif
  _ng1.clear();

  if (_XPol == SOLID_NONE)
    _ngXs.append(_ng2);
  else if (_XPol == SOLID_RIGHT)
  {
#ifdef MANUAL
  {
    K_FLD::FloatArray coord=_coord;
    __compact_and_join(_ng2, coord);
    K_FLD::IntArray connect;
    _ng2.export_to_array(connect);
    MIO::write("sub.plt", coord, connect, "NGON");
  }
#else
//  {
//    K_FLD::FloatArray coord=_coord;
//    __compact_and_join(_ngXh, coord);
//    K_FLD::IntArray connect;
//    _ngXh.export_to_array(connect);
//    MIO::write("ngxhonelayer.plt", coord, connect, "NGON");
//  }
    _ngXh.append(_ng2);
#endif
  }

  _ng2.clear();

  //
  ngon_type::clean_connectivity(_ngXs, _coord, 3); //glue and clean the soft part :  required for classification
  ngon_type::clean_connectivity(_ngXh, _coord, 3); //glue and clean the soft part :  required for classification
  
#ifdef DEBUG_BOOLEAN
  NGON_DBG::write("ngXs1", ACoordinate_t(_coord), _ngXs);
   NGON_DBG::write("ngXh1", ACoordinate_t(_coord), _ngXh);
#endif

#ifdef FLAG_STEP
  std::cout << "NGON Boolean : aggregate the soft bits : " << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_BOOLEAN
  NGON_DBG::write("ngXsG", ACoordinate_t(_coord), _ngXs);
  NGON_DBG::write("ngXhG", ACoordinate_t(_coord), _ngXh);
  NGON_DBG::write_external_phs("ngXsG_ex.plt", ACoordinate_t(_coord), _ngXs);
  std::vector<bool> keep(_ngXs.PHs.size(), true);
  for (size_t i=0; i < _ngXs.PHs.size(); ++i) keep[i]=(_ngXs.PHs._type[i] != INITIAL_SKIN);
  NGON_DBG::write("ngXsGin", ACoordinate_t(_coord), _ngXs, &keep);
#endif
  
  return OK;
}

#define NODE_TO_NEIGH(i) ((i==0) ? 2 : (i==1) ? 0 : 1)

///
TEMPLATE_COORD_CONNECT
typename NGON_BOOLEAN_CLASS::eRetCode
NGON_BOOLEAN_CLASS::__process_intersections
(ngon_unit& wPGs, E_Int& nb_pgs1, ngon_unit& extrawPGs, K_FLD::IntArray& connectT3, Vector_t<E_Int>& nT3_to_PG, Vector_t<E_Int> & priority)
{
 #ifdef FLAG_STEP
  std::cout << "NGON Boolean : nb points : "  << _aCoords1.size() + _crd2.cols() << std::endl;
  std::cout << "NGON Boolean : nb PH(PG) M1 : "  << _cNGON1[2+_cNGON1[1]] << "(" << _cNGON1[0] << ")" << std::endl;
  std::cout << "NGON Boolean : nb PH(PG) M2 : "  << _cNGON2[2+_cNGON2[1]] << "(" << _cNGON2[0] << ")" << std::endl;
  
  chrono c, c0;
  c.start();
  c0.start();
#endif
  
  eRetCode er = __get_working_PGs(_XPol, _MPol, wPGs, nb_pgs1, extrawPGs);
  if (er)
    return er;
  
#ifdef FLAG_STEP
  std::cout << "NGON Boolean : __get_working_PGs : " << c.elapsed() << std::endl;
  c.start();
#endif

  Vector_t<E_Int> oT3_to_PG;
  E_Int err = ngon_type::template triangulate_pgs<DELAUNAY::Triangulator>(wPGs, _coord, connectT3, oT3_to_PG, _do_not_shuffle, _improve_quality_by_swapping);
  if (err)
    return ERROR;
  
  // Get the nb_tleft : triangle of priorized operand (left or right -if preserve right -)
  E_Int nb_tleft=0;
  for (E_Int i=0; i < connectT3.cols(); ++i)
  {
    if (oT3_to_PG[i]>= nb_pgs1)
    {nb_tleft = i; break;}
  }
  
   if (_MPol == PRESERVE_RIGHT)//swap PG1 and PG2 triangles
  {
    K_FLD::IntArray tmp;
    Vector_t<E_Int> o_t3_tpg_tmp;
    Vector_t<E_Int> selids;
    E_Int newnbtleft=0;    
    for (E_Int i=nb_tleft; i < connectT3.cols(); ++i)
    {
      tmp.pushBack(connectT3.col(i), connectT3.col(i)+3);
      o_t3_tpg_tmp.push_back(oT3_to_PG[i]);
    }
    newnbtleft = tmp.cols();
    for (E_Int i=0; i < nb_tleft; ++i)
    {
      tmp.pushBack(connectT3.col(i), connectT3.col(i)+3);
      o_t3_tpg_tmp.push_back(oT3_to_PG[i]);
    }
    nb_tleft=newnbtleft;
    connectT3=tmp;
    oT3_to_PG = o_t3_tpg_tmp;
  }
  
#ifdef DEBUG_BOOLEAN
  // Check what happened to a given PG
  /*E_Int PGi = 5080;
  // get all the oT3s for that PG
  std::vector<E_Int> oT3s;
  for (size_t i=0; i < connectT3.cols(); ++i)
    if (oT3_to_PG[i] == PGi)oT3s.push_back(i);
  NGON_DBG::draw_PG_to_T3(PGi, oT3_to_PG, _coord, connectT3);*/
#endif
  
#ifdef FLAG_STEP
  std::cout << "NGON Boolean : __triangulate " << wPGs.size() << " polygons onto " << connectT3.cols() << " triangles : " << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_BOOLEAN
//Check skin
{
  /*E_Int nT3nb = connectT3.cols(), PGi;
  Vector_t<bool> is_skin(nT3nb, false);
  for (size_t Ti = 0; Ti < nT3nb; ++Ti)
  {
    const E_Int&PGi = oT3_to_PG[Ti];
    is_skin[Ti] = wPGs._external[PGi];
  }
  // skin
  {
    K_FLD::IntArray cT3;
    Vector_t<E_Int> colors;
    for (size_t Ti = 0; Ti < nT3nb; ++Ti)
    {
      if (is_skin[Ti])
      {
        cT3.pushBack(connectT3.col(Ti), connectT3.col(Ti)+3);
        colors.push_back((oT3_to_PG[Ti] < nb_pgs1) ? 0 : 1);
      }
    }
    MIO::write("rawT3skin.mesh", _coord, cT3, "TRI", 0, &colors);
    NGON_DBG::write_wired("wired.mesh", _coord, cT3, true);
  }*/
  MIO::write("allT3.mesh", _coord, connectT3, "TRI", 0, &oT3_to_PG);
}
#endif
    
#ifdef FLAG_STEP
  c.start();
  if (chrono::verbose > 0) std::cout << "NGON Boolean : __conformize ..." << std::endl;
#endif

  // Set priorities for duplicates
  priority.resize(connectT3.cols(), 0);
  for (size_t i=0; i<priority.size(); ++i)
  {
    const E_Int& PGi = oT3_to_PG[i];
    priority[i] = (wPGs._type[PGi] != INITIAL_SKIN) ? 0 : (PGi < nb_pgs1) ? 1 : 2;
  }
  
  Vector_t<E_Int> nT3_to_oT3; // new triangles to old triangles.
  err = __conformize(_coord, connectT3, nT3_to_oT3, nb_tleft/*X0*/, priority);
  if (err) return ERROR;
  
#ifdef DEBUG_BOOLEAN

  E_Int imin;
  E_Float minq;
  K_CONNECT::GeomAlgo<K_MESH::Triangle>::min_quality(_coord, connectT3, minq, imin);
  
  if (minq < _tolerance)
  {
    std::cout << "ERROR : some swapping failed for Triangle (ancestor) : " << imin << "(" << nT3_to_oT3[imin] << ") coming from PG : " <<  oT3_to_PG[nT3_to_oT3[imin]] << std::endl;
    return ERROR;
  }

#ifdef DEBUG_CONFORMIZER
  E_Int Ti = TRI_debug::get_T3_index(connectT3, TRUPLE);
  E_Int a = (Ti != E_IDX_NONE) ? nT3_to_oT3[Ti] : E_IDX_NONE;
#endif
#endif
    
#ifdef FLAG_STEP
  std::cout << "NGON Boolean : __conformize : " << c.elapsed() << std::endl;
#endif
    
  E_Int nT3nb = connectT3.cols();
  
  // fixme : should split on connex sets of confomized triangles (using K_CONNECT::EltAlgo<ElementType>::getNeighbours)
  // for the general case for building poly cells from a general set of intersecting PGs.
  
  // T3 to original PG
  nT3_to_PG.resize(nT3nb);
  for (E_Int i = 0; i < nT3nb; ++i)
    nT3_to_PG[i]= oT3_to_PG[nT3_to_oT3[i]];
  
#ifdef DEBUG_BOOLEAN
  // Check what happened to a given PG
  //NGON_DBG::draw_PG_to_T3(PGi, nT3_to_PG, _coord, connectT3);
#endif
  
#ifdef FLAG_STEP
  c.start();
#endif
  
#ifdef DEBUG_BOOLEAN
  // process duplicates
  err = __process_duplicates(wPGs, connectT3, nT3_to_PG, priority);
  if (err) return ERROR;
#endif
  
  // Compute the normals and replace it by the PG ancestor for nearly degen triangles.
  K_FLD::ArrayAccessor<Coordinate_t> ac(_coord);
  K_FLD::ArrayAccessor<K_FLD::IntArray> acT3n(connectT3);
  __compute_normals(ac, acT3n, wPGs, nT3_to_PG, _normals);
  
#ifdef FLAG_STEP
  std::cout << "NGON Boolean : __process_duplicates : " << c.elapsed() << std::endl;
  std::cout << "NGON Boolean : nb of triangles in the conformal cloud : " << connectT3.cols() << std::endl;
#endif
    
#ifdef DEBUG_BOOLEAN
  // conformized
  {
    Vector_t<E_Int> partCol, PGCol;
    for (size_t Ti = 0; Ti < connectT3.cols(); ++Ti)
    {
      partCol.push_back((nT3_to_PG[Ti] < nb_pgs1) ? 0 : 1);
      PGCol.push_back(nT3_to_PG[Ti]);
    }
    MIO::write("conformized_partcol.mesh", _coord, connectT3, "TRI", 0, &partCol);
    MIO::write("conformized_PGcol.mesh", _coord, connectT3, "TRI", 0, &PGCol);
  }
#endif
  
#ifdef DEBUG_BOOLEAN
  {
  //NGON_DBG::draw_PGT3s(_coord, extrawPGs);
  //NGON_DBG::draw_PGT3s(_coord, wPGs);
  //NGON_DBG::draw_PG_to_T3(886, nT3_to_PG, _coord, connectT3);
  /*E_Int PGi = 5080;
  NGON_DBG::draw_PGT3(_coord, wPGs, PGi);
  std::vector<bool> flag(connectT3.cols(), false);
  for (size_t i=0; i < connectT3.cols(); ++i)flag[i]=(nT3_to_PG[i]==PGi);
  MIO::write("PGi.mesh", _coord, connectT3, "TRI", &flag);*/
  //TRI_debug::draw_connected_to_T3(_coord, connectT3, 12330, 12409, 7205);
  //E_Int Ti = TRI_DBG::get_T3_index(connectT3, TRUPLE);
  //E_Int a = (Ti != E_IDX_NONE) ? nT3_to_oT3[Ti] : E_IDX_NONE;
  }
#endif
  
  return OK;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__process_duplicates(const ngon_unit&wPGs, K_FLD::IntArray& connectT3, Vector_t<E_Int>& nT3_to_PG, Vector_t<E_Int>& priority)
{
  // Remove any duplicated to have a unique represent for each T3 before duplication for neighbouring (which can fail otherwise).
  std::vector<E_Int> dupIds;
  bool has_duplis = K_CONNECT::MeshTool::detectDuplicated(connectT3, dupIds, false);
  
  // THIS FUNCTION MUST BE USELESS AS DUPLICATES ARE PROCESSED AT THE CONFORMIZER STAGE
  assert (!has_duplis);
  

  if (!has_duplis)
    return 0;
  
#ifdef DEBUG_BOOLEAN
  Vector_t<bool> keep(connectT3.cols(), false);
  for (size_t i=0; i < dupIds.size(); ++i)
    keep[i]=(dupIds[i] != i);
  MIO::write("dups.mesh", _coord, connectT3, "TRI", &keep);
  TRI_debug::write_wired("dupnorms.mesh", _coord, connectT3, true, 0, &keep);
#endif

  std::vector<E_Int> priorityOut;
  bool old_is_skin, new_is_skin;
  for (size_t i=0; i < dupIds.size(); ++i)
  {
    E_Int & dupId = dupIds[i];
    if (dupId==E_Int(i)) continue;
    assert (dupId < i);
            
    //
    E_Int pgkeep = zABS(nT3_to_PG[dupId]);
    old_is_skin = (wPGs._type[nT3_to_PG[i]] == INITIAL_SKIN);
    new_is_skin = (wPGs._type[pgkeep] == INITIAL_SKIN);

    // if the one to remove is a skin and 
    // if the one to keep is not a skin or its PG ancestor has a smaller id, swap the one to remove and to keep
    if ( old_is_skin && (!new_is_skin || (new_is_skin && (pgkeep < nT3_to_PG[i])) ))
    {
      dupIds[dupId]=i;
      dupId=i;
    }
    // Flag this triangle as having a duplicate
    //nT3_to_PG[dupId] = -zABS(nT3_to_PG[dupId]);
    priorityOut.push_back(dupId);

  }
  
  //compress nT3_to_PG accordingly.
  K_CONNECT::unchanged pred(dupIds);
  K_CONNECT::IdTool::compress(connectT3, pred);
  std::vector<E_Int> nids;
  K_CONNECT::IdTool::compress(nT3_to_PG, pred, nids);
  
  for (size_t i = 0; i < priorityOut.size(); ++i)
    priorityOut[i] = nids[priorityOut[i]];
  
  priority.insert(priority.end(), priorityOut.begin(), priorityOut.end());
    
#ifdef DEBUG_BOOLEAN
  K_FLD::IntArray tmp;
  std::vector<E_Int> colors;
  for (E_Int i=0; i < priorityOut.size(); ++i)
  {
    const E_Int& Ki = priorityOut[i];
    tmp.pushBack(connectT3.col(Ki),connectT3.col(Ki)+3);
    colors.push_back(nT3_to_PG[Ki]);
  }
  MIO::write("ambiguousPGs.mesh", _coord, tmp, "TRI", 0, &colors);
  TRI_debug::write_wired("dupnorms.mesh", _coord, tmp, true);
#endif
  return 0;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__build_connect_hard
(const K_FLD::FloatArray& coord, ngon_unit& extrawPGs, E_Int nb_pgs1, const K_FLD::IntArray& connectT3,
 K_FLD::IntArray& connectHard, Vector_t<E_Int>& nT3_to_PG, Vector_t<E_Int>& is_skin)
{
  // Keep in connectT3 only solid skin T3s and flag hard_edges (to preserve soft impact AND original PG boundaries)
  Vector_t<E_Int> soft_colors;
  {
    E_Int nT3nb = connectT3.cols();
    Vector_t<bool> keep(nT3nb, false);
    
    for (E_Int Ti = 0; Ti < nT3nb; ++Ti)
    {
      E_Int PGi = nT3_to_PG[Ti];  
      keep[Ti]= (PGi >= nb_pgs1);
    }
    
    // flag hard edges (soft impact)
    std::set<K_MESH::NO_Edge> hedges;
    K_MESH::NO_Edge noE;
    K_FLD::IntArray::const_iterator pS;
    for (E_Int Ti = 0; Ti < nT3nb; ++Ti)
    {
      if (keep[Ti] == false) // Soft operand
      {
        pS = connectT3.col(Ti);
        for (size_t i=0;i < 3; ++i)
        {
          noE.setNodes(*(pS+i), *(pS+(i+1)%3));
          hedges.insert(noE);
        }
      }
    }
        
    Vector_t<E_Int> nids;
    connectHard = connectT3;
    K_FLD::IntArray::compact(connectHard, keep, nids);
    K_CONNECT::IdTool::compact(nT3_to_PG, nids); // PGs are referring to wNG/F2E
    K_CONNECT::IdTool::compact(_normals, nids);
    
    is_skin.clear();
    is_skin.resize(connectHard.cols(), 2);
    
#ifdef DEBUG_BOOLEAN
    MIO::write("solid_skin.mesh", _coord, connectHard, "TRI", 0, &nT3_to_PG);
    {
    K_FLD::IntArray conedges;
    Vector_t<E_Int> colors;
    E_Int nT3nb = connectHard.cols();
    for (size_t Ti = 0; Ti < nT3nb; ++Ti)
    {
      pS = connectHard.col(Ti);
      for (size_t i=0;i < 3; ++i)
      {
        noE.setNodes(*(pS+i), *(pS+(i+1)%3));
        conedges.pushBack(noE.begin(), noE.end());
        if (hedges.find(noE) != hedges.end())
          colors.push_back(1);
        else
          colors.push_back(0);
      }
    }
    MIO::write("wired_soft.mesh", _coord, conedges, "BAR", 0, &colors);
    }
#endif
    
    //Modify color to take into account soft impact
    K_FLD::IntArray neighbors;
    K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectHard, neighbors);
  
    // update neighbors by taking into account hedges
    nT3nb = connectHard.cols();
    for (E_Int Ti = 0; Ti < nT3nb; ++Ti)
    {
      pS = connectHard.col(Ti);
      for (size_t i=0;i < 3; ++i)
      {
        noE.setNodes(*(pS+i), *(pS+(i+1)%3));
        E_Int& Tneigh = neighbors(NODE_TO_NEIGH(i), Ti);
        if (Tneigh == E_IDX_NONE)
          continue;
        if ( (hedges.find(noE) != hedges.end()) || (nT3_to_PG[Tneigh] != nT3_to_PG[Ti]) ) // or if boundary between to old PGs
          Tneigh=E_IDX_NONE;
      }
    }
  
    //
    K_CONNECT::EltAlgo<K_MESH::Triangle>::coloring_pure(neighbors, soft_colors);
    
#ifdef DEBUG_BOOLEAN
    //MIO::write("solid_skin_with_soft.mesh", _coord, connectHard, "TRI", 0, &soft_colors);
#endif
  }
  
  // Move extrawPGs after conformizing (node _history).
  extrawPGs.change_indices(_nodes_history, true);
  extrawPGs.remove_consecutive_duplicated();//consecutive version to preserve PGs topology

#ifdef DEBUG_BOOLEAN
  //NGON_DBG::draw_PGT3s(_coord, extrawPGs);
  //E_Int Tj = TRI_debug::get_T3_index(cT3, 20389, 17508,17522);
  //NGON_DBG::draw_PG(_coord, extrawPGs, 8);
#endif
  
  // Refine the open polygons of the solid layer
  __refine_open_PGs(connectHard, nT3_to_PG, extrawPGs);

  // Triangulate extrawPGs and append to connectHard.
  Vector_t<E_Int> oT3_to_PG; //PG referring to extrawPGs/extraF2E
  K_FLD::IntArray cT3;
  E_Int err = ngon_type::template triangulate_pgs<DELAUNAY::Triangulator>(extrawPGs, _coord, cT3, oT3_to_PG, _do_not_shuffle, true /*improve quality*/); // improve qual of bulks triangulation improve robustness
  if (err) return err;

#ifdef DEBUG_BOOLEAN
  MIO::write("refined_open_layer.mesh", _coord, cT3, "TRI", 0, &oT3_to_PG);
#endif
    
  // Transfer skin info to the triangles
  {
    Vector_t<E_Int> skin_tmp;
    skin_tmp.reserve(cT3.cols()+is_skin.size());
    for (E_Int i=0; i < cT3.cols(); ++i)
    {
      const E_Int& flagi=extrawPGs._type[oT3_to_PG[i]];
      skin_tmp.push_back(flagi); //can be 0 or 3
    }
    skin_tmp.insert(skin_tmp.end(), is_skin.begin(), is_skin.end());
    is_skin=skin_tmp;
  }
  // compute the extra normals
  K_FLD::FloatArray extraNorms;
  K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd(_coord);
  K_FLD::ArrayAccessor<K_FLD::IntArray> acT3o(cT3);
  __compute_normals(acrd, acT3o, extrawPGs, oT3_to_PG, extraNorms);
    
  // append the hard skin to the open layer
  extraNorms.pushBack(_normals);
  _normals = extraNorms;
  cT3.pushBack(connectHard);
  connectHard=cT3;
  
  //shift soft_colors before appending
  E_Int shift = 1 + *std::max_element(oT3_to_PG.begin(), oT3_to_PG.end());//fixme : don't we already know this numer ? like nb_pgs...
  E_Int shift2 = 1 + *std::max_element(nT3_to_PG.begin(), nT3_to_PG.end());
  shift = std::max(shift, shift2);

  K_CONNECT::IdTool::shift(soft_colors, shift);
  oT3_to_PG.insert(oT3_to_PG.end(), soft_colors.begin(), soft_colors.end());

  E_Int max_id_new_pgs = 1 + *std::max_element(soft_colors.begin(), soft_colors.end());
  E_Int old_sz = _anc_PG.cols();
  assert(max_id_new_pgs >= old_sz);//ensure not reducing _anc_PG
  
  _anc_PG.resize(2, max_id_new_pgs, E_IDX_NONE);
  for (E_Int i = 0; i < nT3_to_PG.size(); ++i)
    _anc_PG(1, soft_colors[i]) = nT3_to_PG[i];
  
  assert (is_skin.size() == connectHard.cols());
    
  nT3_to_PG = oT3_to_PG;
    
#ifdef DEBUG_BOOLEAN
  MIO::write("tapped_layer.mesh", _coord, connectHard, "TRI", 0, &nT3_to_PG);
#endif
  return 0;
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__refine_open_PGs
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& nT3_to_PG, ngon_unit& extrawPGs)
{
  // Get the skin nodes : among them are the free nodes of extrawPGs.
  Vector_t<E_Int> skin_nodes;
  connectT3.uniqueVals(skin_nodes); //zero based
  // Flag the extrawPGs free nodes
  E_Int mN = 1+*std::max_element(skin_nodes.begin(), skin_nodes.end());
  
  Vector_t<bool> is_skin_node(mN, false), free_nodes(mN, false);
  for (size_t i = 0; i < skin_nodes.size(); ++i)
  {
    is_skin_node[skin_nodes[i]] = true;
  }
  
  Vector_t<E_Int> layer_nodes;
  extrawPGs.unique_indices(layer_nodes);
  
  for (size_t i = 0; i < layer_nodes.size(); ++i)
  {
    const E_Int& Ni = layer_nodes[i];
    free_nodes[Ni-1] = (Ni-1 < mN) ? is_skin_node[Ni-1] : false;
  }
    
  // Separate skin T3s by PG color.
  
#ifdef DEBUG_BOOLEAN
  E_Int colmax = 1 + *std::max_element(nT3_to_PG.begin(), nT3_to_PG.end());
  assert (colmax != E_IDX_NONE);
#endif
  
  std::map <E_Int, std::vector<E_Int> > skinPGT3s;
  for (E_Int i = 0; i < connectT3.cols(); ++i)
    skinPGT3s[nT3_to_PG[i]].push_back(i);
  
  K_CONT_DEF::int_pair_vector_type boundaries;
  K_FLD::IntArray connectPG, neighbors;
  std::map<E_Int, K_CONT_DEF::int_pair_type> node_to_nodes;
  Vector_t<E_Int> polyLine, sorted_nodes;
  std::map<K_MESH::NO_Edge, Vector_t<E_Int> > edge_to_refined_edge; //starts at 1
  K_MESH::NO_Edge E;

  std::map <E_Int, std::vector<E_Int> >::const_iterator itC;
  for (itC = skinPGT3s.begin(); itC != skinPGT3s.end(); ++itC)
  {  
    //E_Int col = itC->first;
    // Get the PG contour.
    connectPG.clear();
    connectPG.append_selection(connectT3, itC->second);
    
#ifdef DEBUG_BOOLEAN
    //MIO::write("connectPG.mesh", _coord, connectPG, "TRI");
#endif
    
    //Build the neighbourhood matrix.
    K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectPG, neighbors);
    // Get the contour.
    K_CONNECT::MeshTool::getBoundaryT3Mesh(connectPG, neighbors, boundaries);
    // Loop over the contour, create a single contour node list.
    node_to_nodes.clear();
    size_t nb_bounds(boundaries.size()), nb_free(0);
    bool free_consecutive = false;
    for (size_t i = 0; i < nb_bounds; ++i)
    {
      const E_Int& K = boundaries[i].first;
      const E_Int& n = boundaries[i].second;
      const E_Int& Ni = connectPG((n+1)%3,K);
      const E_Int& Nj = connectPG((n+2)%3, K);
      
      node_to_nodes[Ni].second=Nj;
      node_to_nodes[Nj].first=Ni;

      if (free_nodes[Ni])
        ++nb_free;
      free_consecutive |= (free_nodes[Ni] && free_nodes[Nj]);
    }

    //fast return : no refinement required if:
    //1) only 2 free nodes, that are consecutive in top of that
    //2) all are free nodes
    if (nb_free == 2 && free_consecutive)
      continue;
    if (nb_free == nb_bounds)
      continue;

    // Create a single sorted node list
    sorted_nodes.clear();
    std::map<E_Int, K_CONT_DEF::int_pair_type>::const_iterator it = node_to_nodes.begin();
    E_Int N=it->first, Nk;

    for (; it != node_to_nodes.end(); ++it)
    {
      Nk=node_to_nodes[N].second;
      sorted_nodes.push_back(Nk);
      N=Nk;
    }

    assert(sorted_nodes.size() == nb_bounds);
    
    // Go through the list, detect refine nodes and store them.
    bool started = false;
    polyLine.clear();
    
    //start from a free node
    size_t Ifree=0;
    for (size_t i = 0; i < nb_bounds; ++i)
    {
      const E_Int & Ni = sorted_nodes[i];
      if (free_nodes[Ni])
      {
        Ifree=i;
        break;
      }
    }
    
    const size_t& nb_nodes = nb_bounds;
    for (size_t i = 0; i < nb_nodes + 1; ++i)
    {
      const E_Int & Ni = sorted_nodes[(i + Ifree) % nb_nodes];
      if (!started && free_nodes[Ni]) //begin of polyLine
      {
        started=true;
        polyLine.clear();
        polyLine.push_back(Ni+1);
      }
      else if (started && free_nodes[Ni]) //end of polyLine
      {
        started = false;
        
        if (i < nb_nodes) --i; //to start next iteration from the end of this polyLine.
        
        if (polyLine.size() == 1) // no refinement required
          continue;
        
        polyLine.push_back(Ni+1);
        E_Int& Nstart=polyLine[0];
        E_Int& Nend=polyLine[polyLine.size()-1];
        
        if (Nend < Nstart) // need a swap
          K_CONNECT::IdTool::reverse_sorting(polyLine);
        
        E.setNodes(Nstart, Nend);//NO_Edge not necessary as we enure to have the first node smaller
        
        if (edge_to_refined_edge.find(E) != edge_to_refined_edge.end())
        {
#ifndef DEBUG_BOOLEAN
          continue;
#endif
          //keep the one with the smallest geom distance : WRONG. DOESNT GIVE CONSISTENT RESULTS SO USE DEVIATION INSTEAD
          /*E_Float d22 = 0.;
          size_t szd = edge_to_refined_edge[E].size();
          for (size_t k = 0; k < szd; ++k)
            d22 += K_FUNC::sqrDistance(_coord.col(edge_to_refined_edge[E][k]-1), _coord.col(edge_to_refined_edge[E][(k + 1) % szd]-1), 3);
          
          szd = polyLine.size();
          E_Float d21 = 0.;
          for (size_t k = 0; (k<szd) && (d21<d22); ++k)
            d21 += K_FUNC::sqrDistance(_coord.col(polyLine[k]-1), _coord.col(polyLine[(k + 1) % szd]-1), 3);*/

#ifdef DEBUG_BOOLEAN
          K_FLD::IntArray cnt1, cnt2;
#endif
          // keep the one with the smallest deviation
          E_Float d22 = 0., lambda;
          size_t szd1 = edge_to_refined_edge[E].size();
          for (size_t k = 0; k < szd1; ++k){
            d22 = std::max(K_MESH::Edge::edgePointMinDistance2<3>(_coord.col(E.node(0)-1), _coord.col(E.node(1)-1), _coord.col(edge_to_refined_edge[E][k]-1), lambda), d22);
#ifdef DEBUG_BOOLEAN            
            E_Int Ed[] = { edge_to_refined_edge[E][k]-1, edge_to_refined_edge[E][(k + 1) % szd1]-1 };
            cnt1.pushBack(Ed, Ed + 2);
#endif
          }
          size_t szd2 = polyLine.size();
          E_Float d21 = 0.;
          for (size_t k = 0; (k < szd2); ++k){
            d21 = std::max(K_MESH::Edge::edgePointMinDistance2<3>(_coord.col(E.node(0)-1), _coord.col(E.node(1)-1), _coord.col(polyLine[k]-1), lambda), d21);
#ifdef DEBUG_BOOLEAN
            E_Int Ed[] = { polyLine[k]-1, polyLine[(k + 1) % szd2]-1 };
            cnt2.pushBack(Ed, Ed + 2);
#endif
          }
         
#ifdef DEBUG_BOOLEAN
          bool a = (d22 < d21);
          bool b = (edge_to_refined_edge[E].size() < polyLine.size());
          if (a&&!b || b && !a)
          {
            MIO::write("pL1.mesh", _coord, cnt1, "TRI");
            MIO::write("pL2.mesh", _coord, cnt2, "TRI");
          }
#endif
          
          if (d22 < d21) //if the stored is better meaning having the smallest geom distance
            continue;
          else if (IS_ZERO(d22-d21) && szd1 <= szd2) //same distance => keep smallest nb of nodes
            continue;
        }
        
        edge_to_refined_edge[E]=polyLine;   
      }
      else if (started) //refine node
        polyLine.push_back(Ni+1);
    } 
  }
  
  // Refine edges according to skin.
  ngon_type::refine_pgs(edge_to_refined_edge, extrawPGs);
}

// skin and conenction skin
TEMPLATE_COORD_CONNECT
typename NGON_BOOLEAN_CLASS::eRetCode
NGON_BOOLEAN_CLASS::__reorient_externals(eInterPolicy XPol, ngon_type& wNG1, ngon_type& wNG2)
{
  DELAUNAY::Triangulator dt;
  if (XPol == BOTH_SURFACE)
  {
#ifdef FLAG_STEP
    chrono c;
    c.start();
#endif

    //size_t sz = wNG1.PHs.size();
    //size_t sz1 = wNG1.PGs.size();
    ngon_unit ngu = wNG1.PGs;
    wNG1 = ngon_type(ngu, true);// one ph for all to have a closed PH to find out the orientation

    bool has_been_rev;
    E_Int er=ngon_type::reorient_skins(dt, _aCoords1.array(), wNG1, has_been_rev);
    if (er)
      return ERROR;

    //now we build one ph for each pg to reduce the working set with __refine_working_area
    ngu = wNG1.PGs;
    wNG1 = ngon_type(ngu, false);//

#ifdef FLAG_STEP
    if (chrono::verbose >0) std::cout << "__get_working_PGs : SURFACE RIGHT MODE : reorient the right surface : " << c.elapsed() << std::endl;
#endif
  }
  else
  {
    bool has_been_rev;
    E_Int er = ngon_type::reorient_skins(dt, _aCoords1.array(), wNG1, has_been_rev);
    if (er)
      return ERROR;
  }

  if ((XPol == SURFACE_RIGHT) || (XPol == BOTH_SURFACE))
  {
#ifdef FLAG_STEP
    chrono c;
    c.start();
#endif

    //size_t sz = wNG2.PHs.size();
    //size_t sz1 = wNG2.PGs.size();
    ngon_unit ngu = wNG2.PGs;
    wNG2 = ngon_type(ngu, true);// one ph for all to have a closed PH to find out the orientation
    bool has_been_rev;
    E_Int er = ngon_type::reorient_skins(dt, _crd2, wNG2, has_been_rev);
    if (er)
      return ERROR;
    
    //now we build one ph for each pg to reduce the working set with __refine_working_area
    ngu = wNG2.PGs;
    wNG2 = ngon_type(ngu, false);//

#ifdef FLAG_STEP
    if (chrono::verbose >0) std::cout << "__get_working_PGs : SURFACE RIGHT MODE : reorient the right surface : " << c.elapsed() << std::endl;
#endif
  }
  else
  {
    bool has_been_rev;
    E_Int er = ngon_type::reorient_skins(dt, _crd2, wNG2, has_been_rev);
    if (er)
      return ERROR;
  }
  return OK;
}

///
TEMPLATE_COORD_CONNECT
typename NGON_BOOLEAN_CLASS::eRetCode
NGON_BOOLEAN_CLASS::__focus_on_intersection_zone
(eInterPolicy XPol, ngon_type& wNG1, ngon_type& wNG2, ngon_type& rNG1, ngon_type& rNG2)
{ 
  
  ACoordinate_t acrd2(_crd2);

#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif 
  
#ifdef DEBUG_BOOLEAN
  std::cout << "focus 1: " << std::endl;
  nb_ghost(wNG2);
#endif
  
  // E1 : RIGHT IS SHELL ?
  if (XPol == SOLID_RIGHT)
  {
#ifdef DEBUG_BOOLEAN
    std::cout << "size of type : " << wNG2.PHs._type.size()<< std::endl;
    std::cout << "nb of PHs : " << wNG2.PHs.size() << std::endl;
#endif

    Vector_t<bool> externals(wNG2.PHs._type.size(), false);
    for (size_t i=0; i < externals.size(); ++i)
      externals[i] = (wNG2.PHs._type[i] == INITIAL_SKIN);

    // append with element with nodes on the skin
    std::set<E_Int> skinNodes;
    for (E_Int i = 0; i < wNG2.PGs.size(); ++i)
    {
      if (wNG2.PGs._type[i] != INITIAL_SKIN)
        continue;
      const E_Int* nodes = wNG2.PGs.get_facets_ptr(i);
      E_Int stride = wNG2.PGs.stride(i);
      for (E_Int j = 0; j < stride; ++j)
        skinNodes.insert(*(nodes + j) - 1);
    }

    // now catch the attached elements
    for (E_Int i = 0; i < wNG2.PHs.size(); ++i)
    {
      if (externals[i])
        continue;

      const E_Int* pgs = wNG2.PHs.get_facets_ptr(i);
      E_Int nb_pgs = wNG2.PHs.stride(i);

      bool is_attached = false;

      for (E_Int p = 0; p < nb_pgs; ++p)
      {
        E_Int PGi = *(pgs + p) - 1;

        const E_Int* nodes = wNG2.PGs.get_facets_ptr(PGi);
        E_Int stride = wNG2.PGs.stride(PGi);
        for (E_Int j = 0; j < stride; ++j)
        {
          E_Int Ni = *(nodes + j) - 1;
          if (skinNodes.find(Ni) != skinNodes.end())
          {
            is_attached = true;
            break;
          }
        }

        if (is_attached) break;
      }

      externals[i] = is_attached;
    }

    __refine_working_area(wNG2, externals, _ng2);

#ifdef DEBUG_BOOLEAN
{
  NGON_DBG::write("working2_E1", acrd2, wNG2);
  NGON_DBG::write("left2_E1", acrd2, _ng2);
}
#endif
  }
  
#ifdef DEBUG_BOOLEAN
  std::cout << "focus 2: " << std::endl;
  nb_ghost(wNG2);
#endif

#ifdef FLAG_STEP
  if (chrono::verbose >0) std::cout << "__focus_on_intersection_zone : E1 : " << c.elapsed() << std::endl;
#endif
  

  // Init : box and tree creation
  Vector_t<K_SEARCH::BBox3D*> boxes1, boxes2;
  K_SEARCH::BBox3D *pool1, *pool2;
  K_SEARCH::BBox3D GBbox1, GBbox2, GBX;

  if (XPol == BOTH_SURFACE)
    __create_pg_boxes(wNG1, _aCoords1, boxes1, pool1, GBbox1);
  else
    __create_ph_boxes(wNG1, _aCoords1, boxes1, pool1, GBbox1);
  
  Vector_t<E_Int> pg_oids; //pg boxes are partial when (XPol == SOLID_RIGHT) so need to get back to original ids.
  if (XPol == SOLID_RIGHT)
    __create_pg_boxes(wNG2, acrd2, boxes2, pool2, GBbox2, &pg_oids, true/*only_skin*/); 
  else if ((XPol == SURFACE_RIGHT) || (XPol == BOTH_SURFACE))
    __create_pg_boxes(wNG2, acrd2, boxes2, pool2, GBbox2);
  else
    __create_ph_boxes(wNG2, acrd2, boxes2, pool2, GBbox2);
  
  bool intersect = K_SEARCH::BBox3D::intersection(GBbox1, GBbox2, GBX);
  if (!intersect){
    delete[] pool1; delete[] pool2;
    boxes1.clear(); boxes2.clear();
    return EMPTY_X;
  }
  #ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__focus_on_intersection_zone : create boxes : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  // E2 : Capture boxes that are inside the common part (in GBX)
  Vector_t<bool> is_in1, is_in2;
  size_t sz1(boxes1.size()), sz2(boxes2.size());

  //immersion case
  bool immersed1=false;
  if (GBX == GBbox1) // fully included
  {
    switch (_Op)
    {
      
      case DIFF :
      {
        switch (XPol)
        {
          case SOLID_RIGHT:
          case BOTH_SURFACE:
          case SURFACE_RIGHT: // empty answer : operand 1 is inside operand 2 
          {
            immersed1=true;
          }
          default: break;
        }
        break;
      }
      default: break; //todo : UNION, INTER
    }
  } 
  else // not fully included
  {
    K_CONNECT::DynZoneExtractor ze;
    ze.getInBox<3>(&GBX.minB[0], &GBX.maxB[0], E_EPSILON, boxes1, is_in1);
  }
  
  if (GBX == GBbox2){} // fully included
  else
  {
    K_CONNECT::DynZoneExtractor ze;
    ze.getInBox<3>(&GBX.minB[0], &GBX.maxB[0], E_EPSILON, boxes2, is_in2);
  }
  
#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__focus_on_intersection_zone : E2 : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  // E3 : 
  {
    bool tree_from_1 = (boxes1.size() < boxes2.size()); // create the tree from the smaller set
    Vector_t<K_SEARCH::BBox3D*>* tree_boxes = tree_from_1 ? &boxes1 : &boxes2;
    Vector_t<K_SEARCH::BBox3D*>* test_boxes = tree_from_1 ? &boxes2 : &boxes1;
    Vector_t<bool>* is_in = tree_from_1 ? &is_in2 : &is_in1;
    Vector_t<bool>* is_in_capt = tree_from_1 ? &is_in1 : &is_in2;
    size_t sz_in = tree_from_1 ? sz2 : sz1;
    size_t sz_capt = tree_from_1 ? sz1 : sz2;
    
    if (is_in->empty()) //i.e. E2 didn't do anything
      is_in->resize(sz_in, true); // true to test all the box
    
    Vector_t<E_Int> oids;
    // WARNING : if SOLID_RIGHT => sz_capt(=sz2) is smaller than the real nb of PGs in wNG2. So need to reinterpret to actual PGs
    if (is_in_capt->empty())
    { //i.e. E2 didn't do anything
      is_in_capt->resize(sz_capt, false); //false to mask those not appearing in captured_boxes
      
      oids.resize(sz_capt);
      for(size_t i=0; i< sz_capt; ++i)oids[i]=i;
    }
    else
    {
      Vector_t<E_Int> nids;  
      // Compact tree boxes to reflect E2 filtering  
      K_CONNECT::IdTool::build_indir(*is_in_capt, nids, oids);
      if (oids.empty()){ // Empty intersection
        delete[] pool1; delete[] pool2;
        boxes1.clear(); boxes2.clear();
        return EMPTY_X;
      }
      K_CONNECT::IdTool::compact(*tree_boxes, nids); // compact to the used boxes accordingly
    }
    

    // Now build the tree
    K_SEARCH::BbTree3D tree(*tree_boxes, E_EPSILON);
    Vector_t<E_Int> captured_boxes;
    
    size_t sz = test_boxes->size(), sz1(0);
    for (size_t i = 0; i < sz; ++i)
    {
      if (!(*is_in)[i]) continue;
      
      sz1 = captured_boxes.size();
      tree.getOverlappingBoxes((*test_boxes)[i]->minB, (*test_boxes)[i]->maxB, captured_boxes);
      (*is_in)[i]=(sz1 < captured_boxes.size());
    }
    
    delete[] pool1; delete[] pool2;
    boxes1.clear(); boxes2.clear();
    
    if (sz1 == 0) //nothing has been captured
    {
      return immersed1 ? IMMERSED_1 : EMPTY_X;
    }
    
    for (size_t i=0; i < captured_boxes.size(); ++i)
    {
      E_Int k = oids[captured_boxes[i]];
      (*is_in_capt)[k]=true;
    }
  }
  
  // Compute the right is_in2 (sized for wnG2.PHs)
  if (XPol == SOLID_RIGHT) // i.e. wNG2 is_in2 is sized for PGs
    __flag_PHs_sharing_nodes_with_selected_PGs(wNG2, is_in2, pg_oids);
  
#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__focus_on_intersection_zone : E3 : " << c.elapsed() << std::endl;
  c.start();
#endif
  
#ifdef DEBUG_BOOLEAN
  std::cout << "focus a: avant refine" << std::endl;
  nb_ghost(wNG2);
#endif
  
  __refine_working_area(wNG1, is_in1, _ng1); //fixme : améliorer cas remap pyramide : trop de travil faita lors qu'au final _ng1 et _ng2 ne contiennent rien...
  __refine_working_area(wNG2, is_in2, _ng2);
  
#ifdef DEBUG_BOOLEAN
  std::cout << "focus a: apres refine" << std::endl;
  nb_ghost(wNG2);
#endif

  wNG1.PHs.updateFacets();
  wNG1.PGs.updateFacets();
  wNG2.PHs.updateFacets();
  wNG2.PGs.updateFacets();
  
#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__focus_on_intersection_zone : __refine_working_area : " << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_BOOLEAN
  {
	NGON_DBG::write("working1_E3", _aCoords1, wNG1);
	NGON_DBG::write("left1_E3", _aCoords1, _ng1);
	NGON_DBG::write("working2_E3", acrd2, wNG2);
	NGON_DBG::write("left2_E3", acrd2, _ng2);
  }
#endif

  return OK;
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__flag_PHs_sharing_nodes_with_selected_PGs(const ngon_type& wNG, Vector_t<bool>& keep, const Vector_t<E_Int>& pg_oids)
{
  assert(keep.size() == pg_oids.size());

  Vector_t<bool> keepPG;
  keepPG.resize(wNG.PGs.size(), false);

  // 1. catch set of involved nodes
  std::set<E_Int> inodes;
  E_Int nb_reduce_pgs = pg_oids.size();
  for (E_Int i = 0; i < nb_reduce_pgs; ++i)
  {
    E_Int PGi = pg_oids[i];

    keepPG[PGi] = keep[i];

    if (!keep[i]) continue;

    E_Int nb_nodes = wNG.PGs.stride(PGi);
    const E_Int* nodes = wNG.PGs.get_facets_ptr(PGi);
    for (E_Int n = 0; n < nb_nodes; ++n)
      inodes.insert(*(nodes + n) - 1);
    
  }

  // 2. catch set of involved PGs
  E_Int nb_pgs = wNG.PGs.size();
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    if (keepPG[i]) continue;

    E_Int nb_nodes = wNG.PGs.stride(i);
    const E_Int* nodes = wNG.PGs.get_facets_ptr(i);

    bool keepit = false;

    for (E_Int n = 0; (n < nb_nodes) && !keepit; ++n)
      keepit = (inodes.find(*(nodes + n) - 1) != inodes.end());

    keepPG[i] = keepit;
  }
      
  wNG.flag_PHs_having_PGs(keepPG, keep); // keep is now referring to PHs
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__create_ph_boxes
(const ngon_type& wNG, const ACoordinate_t& coord, Vector_t<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BBox3D*& pool, K_SEARCH::BBox3D& GBbox)
{
  size_t nb_elts(wNG.PHs.size());
  Vector_t<E_Int> nodes;

  wNG.PGs.updateFacets();
  wNG.PHs.updateFacets();

  boxes.reserve(nb_elts);
  pool = new K_SEARCH::BBox3D[nb_elts];
  
  for (E_Int i = 0; i < 3; ++i){GBbox.minB[i] = K_CONST::E_MAX_FLOAT; GBbox.maxB[i] = -K_CONST::E_MAX_FLOAT;}

  for (size_t i = 0; i < nb_elts; ++i)
  {
    wNG.nodes_ph(i, nodes, true);
    //std::cout << nodes.size() << std::endl;
    K_SEARCH::BBox3D* box = &pool[i];
    box->compute(coord, nodes);
    boxes.push_back(box);
    
    for (E_Int j = 0; j < 3; ++j)
    {
      GBbox.minB[j] = (GBbox.minB[j] > box->minB[j]) ? box->minB[j] : GBbox.minB[j];
      GBbox.maxB[j] = (GBbox.maxB[j] < box->maxB[j]) ? box->maxB[j] : GBbox.maxB[j];
    }
  }
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__create_pg_boxes
(const ngon_type& wNG, const ACoordinate_t& coord, Vector_t<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BBox3D*& pool, K_SEARCH::BBox3D& GBbox, Vector_t<E_Int>*oids, bool only_skin)
{
  wNG.PGs.updateFacets();

  size_t sz, nb_faces(wNG.PGs.size());
  const E_Int* pNodes;
  Vector_t<E_Int> nodes, oidstmp;
  
  boxes.clear();
  
  if (only_skin)
  {
    //std::assert (oids);
    for (size_t i = 0; i < nb_faces; ++i)
    {
      if (wNG.PGs._type[i] != INITIAL_SKIN)
        continue;
      oidstmp.push_back(i);
    }
    *oids=oidstmp;
  }
  else
  {
    oidstmp.resize(nb_faces);
    for (size_t i=0; i < nb_faces; ++i)oidstmp[i]=i;
  }
  
  nb_faces = oidstmp.size(); //update it in case of "only skin"

  boxes.reserve(nb_faces);
  pool = new K_SEARCH::BBox3D[nb_faces];
  
  for (E_Int i = 0; i < 3; ++i){GBbox.minB[i] = K_CONST::E_MAX_FLOAT; GBbox.maxB[i] = -K_CONST::E_MAX_FLOAT;}
  
  size_t fcount(0);
  for (size_t i = 0; i < nb_faces; ++i)
  {
    nodes.clear();
    
    const E_Int& PGi = oidstmp[i];
    
    sz = wNG.PGs.stride(PGi);
    pNodes = wNG.PGs.get_facets_ptr(PGi);
    for (size_t j = 0; j < sz; ++j, ++pNodes)nodes.push_back((*pNodes)-1);//indices convention : start at 1 
    K_SEARCH::BBox3D* box = &pool[fcount];
    box->compute(coord, nodes);
    
    for (E_Int j = 0; j < 3; ++j) 
    {
      if (box->minB[j]==box->maxB[j]) //fixme : what is that ?
      {
        box->minB[j] -= E_EPSILON;
        box->maxB[j] += E_EPSILON;
      }
      GBbox.minB[j] = (GBbox.minB[j] > box->minB[j]) ? box->minB[j] : GBbox.minB[j];
      GBbox.maxB[j] = (GBbox.maxB[j] < box->maxB[j]) ? box->maxB[j] : GBbox.maxB[j];
    }
    
    boxes.push_back(box);
    ++fcount;
  }
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__destroy_tree
(Vector_t<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BbTree3D *tree, K_SEARCH::BBox3D* pool)
{
  delete[] pool;
  boxes.clear();
  delete tree;
  tree=0;
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__refine_working_area
(ngon_type& wNG, const Vector_t<E_Int>& indices, ngon_type& remainingNG)
{
  if (indices.empty())
    return;
  
  Vector_t<bool> flag(wNG.PHs.size(), false);
  for (size_t i = 0; i < indices.size(); ++i)flag[indices[i]]=true;
  
  __refine_working_area(wNG, flag, remainingNG);
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__refine_working_area
(ngon_type& wNG, const Vector_t<bool>& flag, ngon_type& remainingNG)
{
  ngon_type rNG;
  
#ifdef DEBUG_BOOLEAN
  std::cout << "before split:  " << std::endl;
  nb_ghost(wNG);
#endif
  
  wNG.split_phs(wNG, flag, wNG, rNG);
  
#ifdef DEBUG_BOOLEAN
  std::cout << "after split : " << std::endl;
  std::cout << "wNG" << std::endl;
  nb_ghost(wNG);
  std::cout << "rNG" << std::endl;
  nb_ghost(rNG);
#endif

  remainingNG.append(rNG);  
}


///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__discard_holes_by_box(const K_FLD::FloatArray& coord, ngon_type& wNG)
{
  assert (wNG.PHs.size() < wNG.PGs.size()); //i.e. not one ph per pg beacuse we need at least a closed PH

#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif
  
  // 1. Get the skin PGs
  Vector_t<E_Int> oids;
  ngon_unit pg_ext;
  wNG.PGs.extract_of_type (INITIAL_SKIN, pg_ext, oids);
  if (pg_ext.size() == 0)
    return 0; //should be only the case where One of the mesh contains completely the other

#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__discard_holes : extract_of_type : " << c.elapsed() << std::endl;
  c.start();
#endif
    
  // 2. Build the neighbourhood for skin PGs
  ngon_unit neighbors;
  K_MESH::Polygon::build_pg_neighborhood(pg_ext, neighbors);

#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__discard_holes : build_pg_neighborhood : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  // Color to get connex parts
  E_Int nb_connex = 1;
  Vector_t<E_Int> colors;
  if (wNG.PHs.size() > 1)
  {
    K_CONNECT::EltAlgo<K_MESH::Polygon>::coloring (neighbors, colors);
    nb_connex = 1+*std::max_element(colors.begin(), colors.end());
  }

#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__discard_holes : coloring : " << c.elapsed() << std::endl;
  c.start();
#endif
  
  if (nb_connex == 1) //no holes
    return 0;
  
  // Find out which is the external (the one with the biggest bbox)
  Vector_t<K_SEARCH::BBox3D> bbox_per_color(nb_connex);
  //init boxes
  for (E_Int i=0; i < nb_connex; ++i)
  {
    bbox_per_color[i].maxB[0]=bbox_per_color[i].maxB[1]=bbox_per_color[i].maxB[2]=-K_CONST::E_MAX_FLOAT;
    bbox_per_color[i].minB[0]=bbox_per_color[i].minB[1]=bbox_per_color[i].minB[2]=K_CONST::E_MAX_FLOAT;
  }
  
  Vector_t<E_Int> nodes;
  
  E_Int nb_pgex=colors.size();
  K_SEARCH::BBox3D box;
  for (E_Int i=0; i < nb_pgex; ++i)
  {
    const E_Int& s = pg_ext.stride(i);
    const E_Int* pN = pg_ext.get_facets_ptr(i);
    
    nodes.clear();
    for (E_Int j = 0; j < s; ++j, ++pN)
      nodes.push_back((*pN)-1);//indices convention : start at 1 

    box.compute(coord, nodes);
    
    K_SEARCH::BBox3D& b = bbox_per_color[colors[i]];
    
    for (size_t j = 0; j < 3; ++j) 
    {
      b.minB[j] = (b.minB[j] > box.minB[j]) ? box.minB[j] : b.minB[j];
      b.maxB[j] = (b.maxB[j] < box.maxB[j]) ? box.maxB[j] : b.maxB[j];
    }
  }
  // And the hole color are..
  std::set<E_Int> hole_colors;
  for (E_Int i = 0; i < nb_connex; ++i)
  {
    for (E_Int j = i + 1; j < nb_connex; ++j)
    {
      if (BbTree3D::box1IsIncludedinbox2(&bbox_per_color[i], &bbox_per_color[j], E_EPSILON))
        hole_colors.insert(i);
      else if (BbTree3D::box1IsIncludedinbox2(&bbox_per_color[j], &bbox_per_color[i], E_EPSILON))
        hole_colors.insert(j);
    }
  }

  // Reset flag for holes
  for (E_Int i=0; i < nb_pgex; ++i)
  {
    if (hole_colors.find(colors[i]) != hole_colors.end())
      wNG.PGs._type[oids[i]]=INNER;
  }
  wNG.flag_external_phs(INITIAL_SKIN);//update consistently PHs flags
  
  return 0;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__discard_prescribed_polygons(const K_FLD::FloatArray& coord, ngon_type& wNG, const Vector_t<E_Int>& PGlist)
{
  // Reset flag for holes
  E_Int nb_pgs = PGlist.size();
  if (nb_pgs == 0) return 0;

#if defined(DEBUG_BOOLEAN) || defined(DEBUG_W_PYTHON_LAYER)
  ngon_unit pgs;
  E_Int minid=100000000;
  E_Int maxid = -1;
  std::ofstream of("pglist2.txt");
#endif
  
  for (E_Int i=0; i < nb_pgs; ++i)
  {
    const E_Int& PGi = PGlist[i];
    wNG.PGs._type[PGi]=INNER;

#if defined(DEBUG_BOOLEAN) || defined(DEBUG_W_PYTHON_LAYER)
    pgs.add(wNG.PGs.stride(PGi), wNG.PGs.get_facets_ptr(PGi));
    minid = std::min(minid, PGi);
    maxid = std::max(maxid, PGi);
    of << PGi << std::endl;
#endif
  }
  wNG.flag_external_phs(INITIAL_SKIN);//update consistently PHs flags
  
#if defined(DEBUG_BOOLEAN)
  K_FLD::ArrayAccessor<Coordinate_t> ac(coord);
  NGON_DBG::write( "disabled.plt", ac, pgs);
  NGON_DBG::draw_PGT3s(coord, pgs);
#endif
  
#if defined(DEBUG_BOOLEAN) || defined(DEBUG_W_PYTHON_LAYER)
  std::cout << "range for pg walls : [" << minid << "," << maxid << "]" << std::endl;
  of.close();
#endif
  
  return 0;
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__duplicate_orient
(K_FLD::IntArray& connectT3)
{
  E_Int nbT30 = connectT3.cols();
  K_FLD::IntArray tmp = connectT3;
  __swap(tmp);
  connectT3.pushBack(tmp);
  
  K_FLD::FloatArray norms=_normals; //fixme : need currently a copy for back pushing itself
  _normals.pushBack(norms);
  for (E_Int i=nbT30; i < _normals.cols(); ++i)
  {
    _normals(0,i)=-_normals(0,i);
    _normals(1,i)=-_normals(1,i);
    _normals(2,i)=-_normals(2,i);
  }
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__compute_normals
(const K_FLD::ArrayAccessor<K_FLD::FloatArray>& acrd, const K_FLD::ArrayAccessor<K_FLD::IntArray>& acnt, 
 const ngon_unit& PGs, const Vector_t<E_Int> T3_to_PG, K_FLD::FloatArray& T3normals)
{
  size_t nb_t3s = acnt.size();
  
  T3normals.clear();
  T3normals.reserve(3, nb_t3s);
  
  std::map<E_Int, E_Int> pgid_to_PGnorm;
  std::map<E_Int, E_Int>::iterator it;
  K_FLD::FloatArray PGnormals;
  
  E_Int t3[3], PGi, Normi;
  E_Float p0[3], p1[3], p2[3], normal[3];
    
  //
  for (size_t i = 0; i < nb_t3s; ++i)
  {
    acnt.getEntry(i, t3);
    acrd.getEntry(t3[0], p0);
    acrd.getEntry(t3[1], p1);
    acrd.getEntry(t3[2], p2);
    
    K_MESH::Triangle::normal(p0, p1, p2, normal); 
    
    E_Float l2 = ::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    
    if (::fabs(l2 - 1.) < E_EPSILON) // NOT degen
      T3normals.pushBack(normal, normal+3);
    else
    {
      PGi = T3_to_PG[i];
      it = pgid_to_PGnorm.find(PGi);
      if (it == pgid_to_PGnorm.end()) // PG normal not computed yet
      {
        const E_Int* pN = PGs.get_facets_ptr(PGi);
        E_Int stride = PGs.stride(PGi);
    
        typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> CordAcc_t;
        K_MESH::Polygon::normal<CordAcc_t, 3>(acrd, pN, stride, 1, &normal[0]); 
        
        PGnormals.pushBack(normal, normal+3);
        pgid_to_PGnorm[PGi]=PGnormals.cols()-1;
      }
      
      Normi = pgid_to_PGnorm[PGi];
      const E_Float* pg_normal = PGnormals.col(Normi);
      
      T3normals.pushBack(pg_normal, pg_normal+3);
    }
  }  
  
//#ifdef DEBUG_BOOLEAN
//  {
//    E_Int PGi = 5;
//    ///NGON_DBG::draw_PGT3(acrd.array(), PGs, PGi, (PGnormals.cols() == 0) ? 0: &PGnormals);
//    Vector_t<E_Int> oids;
//    const K_FLD::IntArray& cT3 = acnt.array();
//    TRI_debug::get_same_ancestor_T3s(PGi, 1, _coord, cT3, T3_to_PG, oids);
//    Vector_t<bool> keep(cT3.cols(), false);
//    for (size_t i=0; i < oids.size(); ++i)keep[oids[i]]=true;
//    std::ostringstream o;
//    o << "normals_for_" << PGi << ".mesh";
//    TRI_debug::write_wired(o.str().c_str(), _coord, cT3, T3normals, 0, &keep, true);
//  }
//#endif
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__sort_T3_sharing_an_edge
(E_Int E0, E_Int E1, E_Int shift, 
 const K_FLD::FloatArray& normals, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3,
 Vector_t<E_Int>& T3indices)
{
  _palmares.clear();
  
  E_Int sz(T3indices.size()), i, N1;
  const E_Int *pSi;
  E_Float q;
  
  E_Int& K0 = T3indices[0];
  
#ifdef DEBUG_BOOLEAN
//#define zeEdge (E0 == 49935 && E1 == 47262) || (E0 == 47262 && E1 == 49935)
#endif
  
  for (E_Int j = 1; j < sz; ++j)
  {
#ifdef DEBUG_BOOLEAN
    /*if (zeEdge)
    {
      std::cout << "dummy for hooking a breakpoint" << std::endl;
    }*/
#endif
    E_Int& Ki = T3indices[j];
    //finding out Ap and use the right orientation for Ki.
    pSi = connectT3.col(Ki);
    i = K_MESH::Triangle::getLocalNodeId(pSi, E0);
    
    N1 = *(pSi+(i+1)%3);
    if (N1 == E1) // the side to consider is the opposite
      Ki += shift;

    
#ifdef DEBUG_BOOLEAN   
/*if (zeEdge)
  {
    std::cout << "Ki : " << connectT3(0,Ki) << "/" << connectT3(1,Ki) << "/" << connectT3(2,Ki) << std::endl;
    std::cout << "E0 : " << coord(0,E0) << "/" << coord(1,E0) << "/" << coord(2,E0) << std::endl;
    std::cout << "E1 : " << coord(0,E1) << "/" << coord(1,E1) << "/" << coord(2,E1) << std::endl;
    
    E_Float nk[3];
    K_FUNC::crossProduct<3>(normals.col(K0), normals.col(Ki), nk);
    E_Float s2 = K_FUNC::sqrNorm<3>(nk);
    E_Float c = K_FUNC::dot<3>(normals.col(K0), normals.col(Ki));
    E_Float alpha = ::atan2(::sqrt(s2), c); 
    //alpha = K_CONST::E_PI + alpha;
    std::cout << "alpha : " << alpha << std::endl;
    std::cout << std::endl;  
  }*/
#endif
    
    if (sz == 2) break; //no need to sort
    
    q = K_CONNECT::GeomAlgo<K_MESH::Triangle>::angle_measure(normals.col(K0), normals.col(Ki), coord.col(E0), coord.col(E1));
    _palmares.push_back(std::make_pair(q, Ki));
  }
  
  K0 +=shift;// to simplify __update_neigboring algo
  
  if (_palmares.size() > 1)
  {
    std::sort(_palmares.begin(), _palmares.end());
    
    T3indices.clear();
    T3indices.push_back(K0); 
    
    for (size_t i = 0; i < _palmares.size(); ++i)
    {
      E_Int& K = _palmares[i].second;
      T3indices.push_back(K);
#ifdef DEBUG_BOOLEAN
      if (_palmares[i].first == _palmares[(i + 1) % _palmares.size()].first)
      {
        std::cout << "ERROR at edge E0E1 : " << E0 << "/" << E1 << std::endl;
        std::cout << "The conformizer missed some intersections there." << std::endl;
      }
      assert(_palmares[i].first != _palmares[(i + 1) % _palmares.size()].first); //if two q are the same, overlap triangles exist so the conformizer missed some intersections
#endif
    }
  }
    
  
#ifdef DEBUG_BOOLEAN
  /*if ( zeEdge )
  {  
    K_FLD::IntArray sorted_cnt;
    Vector_t<E_Int> colors;
    Vector_t<bool> keep(connectT3.cols(), false);
    sorted_cnt.reserve(3, T3indices.size());
    colors.resize(T3indices.size(), 1);
    for (size_t i = 0; i < T3indices.size(); ++i)
    {
      sorted_cnt.pushBack(connectT3.col(T3indices[i]), connectT3.col(T3indices[i])+3);
      keep[T3indices[i]]=true;
      std::cout << "Triangle : " << T3indices[i] << std::endl;
      if (i > 0)colors[i] = 0;
    }
    std::ostringstream o;
    o << "sorted_on_edge_" << E0 << "_" << E1 << ".mesh";
    MIO::write(o.str().c_str(), coord, sorted_cnt, "TRI", 0, &colors);
    o.str("");
    o << "Wsorted_on_edge_" << E0 << "_" << E1 << ".mesh";
    TRI_debug::write_wired(o.str().c_str(), coord, connectT3, normals, 0, &keep,true);
  }*/
#endif
}

#define THIRD_NODE_POS(pS, E0, E1) ( (*pS != E0 && *pS != E1) ? 0 : (*(pS+1) != E0 && *(pS+1) != E1) ? 1 : 2) 
#define OPPOSITE(K, shift) ( (K < shift) ? K+shift : K-shift )

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__update_neigboring
(E_Int E0, E_Int E1, E_Int shift, const K_FLD::IntArray& connectT3,
 const Vector_t<E_Int>& T3s, K_FLD::IntArray& neighbors)
{
  E_Int Kip1, Kib, Lip1, Lib, sz(T3s.size());
  for (E_Int i = 0; i < sz; ++i)
  {
    const E_Int& Ki = T3s[i];
    Kib = OPPOSITE(Ki, shift);
    
#ifdef DEBUG_BOOLEAN
    //E_Int nib0 = connectT3(0,Kib);
    //E_Int nib1 = connectT3(1,Kib);
    //E_Int nib2 = connectT3(2,Kib);
#endif
    
    Lib = THIRD_NODE_POS(connectT3.col(Kib), E0, E1);
    Kip1 = T3s[(i+1)%sz];
 
#ifdef DEBUG_BOOLEAN
    //E_Int nip0 = connectT3(0,Kip1);
    //E_Int nip1 = connectT3(1,Kip1);
    //E_Int nip2 = connectT3(2,Kip1);
#endif  
    
    Lip1 = THIRD_NODE_POS(connectT3.col(Kip1), E0, E1);
    
    neighbors(Lib, Kib) = Kip1;
    neighbors(Lip1, Kip1) = Kib;
  }
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__build_PHT3s
(const K_FLD::FloatArray& coord,  const K_FLD::IntArray& connectT3, const K_FLD::IntArray& connectT3o, 
 const Vector_t<E_Int>& is_skin, std::map<E_Int, Vector_t<E_Int> >& PHT3s)
{
  
#ifdef DEBUG_BOOLEAN
 
  //E_Int Ti = TRI_DBG::get_T3_index(connectT3, TRUPLE);
  //TRI_debug::draw_connected_to_T3(coord, connectT3, Ti);
  
  //K_FLD::FloatArray crd;
  //K_FLD::IntArray cnt;
  //std::vector<E_Int> clr;
  //TRI_debug::connected_to_T3(coord, connectT3, TRUPLE, crd, cnt, clr);
  //NGON_DBG::remove_T3(crd, cnt, 3445, 5256, 3465);
  //NGON_DBG::remove_T3(crd, cnt, 3445, 5256, 3464);
  //NGON_DBG::remove_T3(crd, cnt, 3157, 5256, 3445);
  //NGON_DBG::remove_T3(crd, cnt, 3440, 3445, 5256);
  //NGON_DBG::remove_T3(crd, cnt, 5256, 3450, 6586);
  //MIO::write("c1.mesh", crd, cnt, "TRI");
  //TRI_DBG::get_T3_neighbors("neighx.mesh", Ti, coord, connectT3o, neighbors, true);
#endif
 
#ifdef FLAG_STEP
  chrono c,c1;
  c.start();c1.start();
#endif

  E_Int err(0);
  K_FLD::IntArray neighbors;
  __build_neighbors_table(coord, connectT3, connectT3o, neighbors);
  
#ifdef FLAG_STEP
  if (chrono::verbose > 0) std::cout << "NGON Boolean : __build_PHT3s : __build_neighbors_table : " << c1.elapsed() << std::endl;
#endif
    
#ifdef FLAG_STEP
  c1.start();
#endif

  err = __assemble_PHT3s(coord, connectT3o, neighbors, PHT3s);
  if (err)
  {
#ifdef FLAG_STEP
  if (chrono::verbose > 0) std::cout << "NGON Boolean : ERROR in __assemble_PHT3s : " << c1.elapsed() << std::endl;
#endif
    return err;
  }
  
#ifdef FLAG_STEP
  if (chrono::verbose > 0) std::cout << "NGON Boolean : __build_PHT3s : __assemble_PHT3s : " << c1.elapsed() << std::endl;
  c1.start();
#endif
  
  // Remove parasites PHs (those having only SKIN (initial and new) triangle : most external and those due to holes

  E_Int PHierr = __remove_parasite_PHT3s(is_skin, PHT3s, connectT3o);
  if (PHierr != -2) //error
  {
#ifdef FLAG_STEP
  if (chrono::verbose > 0) std::cout << "NGON Boolean : ERROR in __remove_parasite_PHT3s : " << c1.elapsed() << std::endl;
#endif
#ifdef DEBUG_BOOLEAN
    if (PHierr > -1)
      TRI_DBG::coloring_frames(coord, connectT3o, neighbors, PHT3s[PHierr][0]);
#endif
    return err;
  }
  
#ifdef FLAG_STEP
  if (chrono::verbose > 0) std::cout << "NGON Boolean : __build_PHT3s : __remove_parasite_PHT3s : " << c1.elapsed() << std::endl;
  c1.start();
#endif
   
#ifdef DEBUG_BOOLEAN
//  std::vector<E_Int> colors(connectT3o.cols(), 0);
//  std::vector<bool> flag(connectT3o.cols(), true);
//  //flag[984]=true;
//  E_Int i=0;
//  for (std::map<E_Int, Vector_t<E_Int> >::const_iterator it = PHT3s.begin(); it != PHT3s.end(); ++it,++i)
//     NGON_DBG::draw_PHT3(coord, connectT3o, PHT3s, i);
//  //TRI_debug::draw_connected_to_T3(coord, connectT3, 783);
//  //E_Int Ti = TRI_DBG::get_T3_index(connectT3o, 3445, 5256, 6586);
//  
//  //NGON_DBG::get_PHT3_neighbors(837, PHT3s, coord, connectT3o, neighbors, true/*both_orient*/);
#endif
  
#ifdef FLAG_STEP
  c1.start();
#endif
  
  // Check if all PHT3s are closed.
  bool check_also_manifoldness = false;
  PHierr = __check_PHT3s_closure(check_also_manifoldness, coord, connectT3o, PHT3s);
    
#ifdef FLAG_STEP
  if (chrono::verbose > 0) std::cout << "NGON Boolean : __build_PHT3s : __check_PHT3s_closure : " << c1.elapsed() << std::endl;
#endif

  if (PHierr != -1)
  {
    std::cout << " faulty (unclosed) PHT3 is : " << PHierr << "(" << PHT3s[PHierr].size() << " T3s)" << std::endl;

#ifdef DEBUG_BOOLEAN
    NGON_DBG::draw_PHT3(coord, connectT3o, PHT3s, PHierr);
    //TRI_DBG::coloring_frames(coord, connectT3o, neighbors, PHT3s[PHierr][0]);
    // its history :
    std::vector<E_Int> PHs1, PHs2; //from 1 and 2
    NGON_DBG::__get_historical_PHs(coord, connectT3o, PHT3s, PHierr, connectT3.cols()/*shift*/, _nb_pgs1, _F2E, _anc_PH_for_PHT3s, _nT3_to_oPG, PHs1, PHs2);
    
    std::cout << "list of ids for M1 : " << std::endl;
    for (size_t i=0; i < PHs1.size(); ++i)
      std::cout << PHs1[i] << " ";
    std::cout << std::endl;
    
    std::cout << "list of ids for M2 : " << std::endl;
    for (size_t i=0; i < PHs2.size(); ++i)
      std::cout << PHs2[i] << " ";
    std::cout << std::endl;
    
    //ngon_type ng1(_cNGON1), ng2(_cNGON2);
    //NGON_DBG::draw_PHs("faulty_Molec1.plt", coord, ng1, PHs1);//fixme : not working
    //NGON_DBG::draw_PHs("faulty_Molec2.plt", coord, ng2, PHs2);
#endif

    return 1;
  }
  
#ifdef FLAG_STEP
  std::cout << "NGON Boolean : __build_PHT3s : " << c.elapsed() << std::endl;
#endif

  return 0;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__build_PHs
(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const K_FLD::IntArray& connectT3,
 const Vector_t<E_Int>& is_skinT3, const Vector_t<E_Int>& nT3_to_PG, ngon_type& ngout)
{
  E_Int err(0);
  
  // Assemble PGT3s.
  typedef std::map<E_Int, Vector_t<E_Int> > pg_to_t3s_t;
  typedef std::map<E_Int, pg_to_t3s_t > ph_to_pgt3s_t;
  ph_to_pgt3s_t PH_to_PGT3s;
  
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif
  
  err = __assemble_PGT3s(PHT3s, nT3_to_PG, PH_to_PGT3s);
  if (err)
    return err;
  
#ifdef DEBUG_BOOLEAN
  //NGON_DBG::draw_PHT3(_coord, connectT3, PH_to_PGT3s, 357);
  //NGON_DBG::draw_PHT3(_coord, connectT3, PH_to_PGT3s, 599);
#endif
  
  // Pass to split PGT3s not exactly shared by one PHT3 on each side.
  //  ____
  // |    |___
  // |    |___
  // |____|
  //      ^
  //      |
  err = __split_multiply_shared_PGT3s(connectT3, PH_to_PGT3s);
  if (err)
    return err;
  
#ifdef DEBUG_BOOLEAN
  //NGON_DBG::draw_PHT3(_coord, connectT3, PH_to_PGT3s, 357);
  //NGON_DBG::draw_PHT3(_coord, connectT3, PH_to_PGT3s, 599);
#endif
  
  // Pass to split non connex PGT3s
  //  __  __ <--- non connex top
  //  | \/  |
  //  |_____|  
  err = __split_non_connex_PGT3s(connectT3, PH_to_PGT3s);
  if (err)
    return err;
 
#ifdef DEBUG_BOOLEAN
  //NGON_DBG::draw_PHT3(_coord, connectT3, PH_to_PGT3s, 357);
  //NGON_DBG::draw_PHT3(_coord, connectT3, PH_to_PGT3s, 599);
#endif
  
#ifdef DEBUG_BOOLEAN
/*{
  E_Int phi, pg_count(0), ph_count;
  Vector_t<E_Int> pg_colors, ph_colors;
  K_FLD::IntArray cT3;
  
  for (ph_to_pgt3s_t::const_iterator itPH = PH_to_PGT3s.begin(); itPH != PH_to_PGT3s.end(); ++itPH)
  {
    phi = itPH->first;
    if (phi != 1)
      continue;
    for (pg_to_t3s_t::const_iterator itPG = PH_to_PGT3s[phi].begin(); itPG != PH_to_PGT3s[phi].end(); ++itPG)
    {
      if (itPG->first != 291)
        continue;
      const Vector_t<E_Int>& T3s = itPG->second;
      for (size_t i = 0; i < T3s.size(); ++i)
        cT3.pushBack(connectT3.col(T3s[i]), connectT3.col(T3s[i])+3);
      pg_colors.resize(cT3.cols(), pg_count++);
    }
    ph_colors.resize(cT3.cols(), ph_count++);
    //std::ostringstream o;
    //o << "PHT3_" << phi << ".mesh";
    //MIO::write(o.str().c_str(), _coord, cT3, "TRI", 0, "&pg_colors);
  }
  MIO::write("PH_toPGT3s_pgcol.mesh", _coord, cT3, "TRI", 0, &pg_colors);
  MIO::write("PH_toPGT3s_phcol.mesh", _coord, cT3, "TRI", 0, &ph_colors);
}*/
#endif
  
#ifdef FLAG_STEP
  std::cout << "NGON_BOOLEAN_CLASS::__build_PHs: prepare agregates : " << c.elapsed() << std::endl;
#endif
    
  // Aggregate into PHs.
  switch (_AggPol)
  {
    case FULL  : err = __aggregate_PHs<FULL>(connectT3, is_skinT3, PH_to_PGT3s, ngout);break;
    case CONVEX: err = __aggregate_PHs<CONVEX>(connectT3, is_skinT3, PH_to_PGT3s, ngout);break;
    case NONE  : err = __aggregate_PHs<NONE>(connectT3, is_skinT3, PH_to_PGT3s, ngout);break;
    default    : err = 1;break;
  }
  
#ifdef DEBUG_BOOLEAN
  /*
  NGON_DBG::write_external_phs("external_PHs.plt", ACoordinate_t(_coord), ngout);
  ngon_type tmp = ngout;
  K_CONNECT::IdTool::negative(tmp.PHs._external);
  NGON_DBG::write_external_phs("internal_PHs.plt", ACoordinate_t(_coord), tmp);
  */
#endif
  
  return err;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__classify_skin_PHT3s
(std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift
#ifdef DEBUG_BOOLEAN
, const K_FLD::IntArray& connectT3o
#endif
)
{
  // flag elements as being Z_1 or Z_IN when connected to the skins.
  __set_skin_PHT3s_zones(PHT3s, is_skin, shift, _zones
#ifdef DEBUG_BOOLEAN
, connectT3o
#endif
);

  __set_PH_history(PHT3s, is_skin, shift, _nb_pgs1, _F2E, _anc_PH_for_PHT3s, true/*soft*/);

  // Remove any Z_2 appearing here (i.e solid skin is closed and has created irrelevant PHT3s.
  if (_XPol == SOLID_RIGHT)
  {
    Vector_t<bool> keep(PHT3s.size(), false);
    std::map<E_Int, Vector_t<E_Int> >::const_iterator it = PHT3s.begin();
    std::map<E_Int, Vector_t<E_Int> >::const_iterator ite = PHT3s.end();
    E_Int z=0;
    Vector_t<E_Int> thrashIds;
    for (; it != ite; ++it, ++z)
    {
      keep[z]=(_zones[z] != Z_2);
      if (_zones[z] == Z_2)
        thrashIds.push_back(it->first);
    }
    if (!thrashIds.empty())
    {
      for (size_t i=0; i < thrashIds.size(); ++i)
        PHT3s.erase(thrashIds[i]);
      
      K_CONNECT::keep pred(keep);
      K_CONNECT::IdTool::compress(_anc_PH_for_PHT3s, pred);
      K_CONNECT::IdTool::compress(_zones, pred);
    }
  }
  
#ifdef DEBUG_BOOLEAN
  Vector_t<E_Int> colz(_zones.size(), 0);
//  {
//    for (size_t i=0; i < _zones.size(); ++i)colz[i]=_zones[i];
//    NGON_DBG::draw_PHT3s("PHT3zones.mesh", _coord, connectT3o, PHT3s, colz);
//  }
  {
    Vector_t<bool> flag(_zones.size(), false);
    for (size_t i=0; i < _zones.size(); ++i)if (_zones[i] == Z_IN){flag[i]=true;}
    NGON_DBG::draw_PHT3s("PHT3zin.mesh", _coord, connectT3o, PHT3s, colz, &flag);
  }
  {
    Vector_t<bool> flag(_zones.size(), false);
    for (size_t i=0; i < _zones.size(); ++i)if (_zones[i] == Z_1){flag[i]=true;}
    NGON_DBG::draw_PHT3s("PHT3z1.mesh", _coord, connectT3o, PHT3s, colz, &flag);
  }
  {
    Vector_t<bool> flag(_zones.size(), false);
    for (size_t i=0; i < _zones.size(); ++i)if (_zones[i] == Z_2){flag[i]=true;}
    NGON_DBG::draw_PHT3s("PHT3z2.mesh", _coord, connectT3o, PHT3s, colz, &flag);
  }
#endif

  return 0;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__classify_soft()
{
  E_Int err(0);
  
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif
      
  // Build NGON neighborhood.
  ngon_unit neighbors;
  std::vector<bool> wall;
  _ngXs.PGs.flag_indices_of_type(INITIAL_SKIN, wall);

#ifdef DEBUG_BOOLEAN
  ngon_unit pg_ext;
  Vector_t<E_Int> oids;
  _ngXs.PGs.extract_of_type(INITIAL_SKIN, pg_ext, oids);
  NGON_DBG::write("walls", ACoordinate_t(_coord), pg_ext);
  NGON_DBG::draw_PGT3s(_coord, pg_ext);
#endif
  
  _ngXs.build_ph_neighborhood(neighbors, wall);
  
#ifdef DEBUG_BOOLEAN
  
  // check if neighbouring is ok
  /*
  size_t nb_phs = _ngXs.PHs.size();
  for (size_t i = 0; i < nb_phs; ++i)
  {
    for (size_t j=0; j < _ngXs.PHs.stride(i); ++j)
    {
      E_Int Neigh = neighbors.get_facet(i,j);
      E_Int Fid = _ngXs.PHs.get_facet(i,j) - 1;
      bool is_ex_pg = _ngXs.PGs._external[Fid];
      
      bool valid = (is_ex_pg && Neigh == E_IDX_NONE) || (!is_ex_pg && Neigh != E_IDX_NONE);
      
      assert (valid);
      
    }
  }
    
  NGON_DBG::draw_PGT3s(_coord, _ngXs.PGs);
  
  size_t sz = neighbors.size();
  std::vector<bool> external(sz, false);
  for (size_t i=0; i < sz; ++i)
  {
    E_Int nb_neigh = neighbors.stride(i);
    
    for (size_t j=0; j < nb_neigh; ++j)
    {
      E_Int Neigh = neighbors.get_facet(i,j);
      if (Neigh == E_IDX_NONE) // is supposed external
      {
        external[i]=true; break;
      }
    }
  }
  NGON_DBG::write("ext_from_neigh", ACoordinate_t(_coord), _ngXs, &external);
  K_CONNECT::IdTool::negative(external);
  NGON_DBG::write("in_from_neigh", ACoordinate_t(_coord), _ngXs, &external);
  
  // check is external PGs are OK
  ngon_unit pgex;
  std::vector<E_Int> oids;
  _ngXs.PGs.extract_externals(pgex, oids);
  NGON_DBG::write("pgex", ACoordinate_t(_coord), pgex);
  */
#endif
    
  // Color.
  Vector_t<E_Int> colors;
  //fixme : ElementType is not relevant for this algo, should be moved somewhere else
  K_CONNECT::EltAlgo<K_MESH::Polygon>::coloring (neighbors, colors);
  E_Int colmax = *std::max_element(colors.begin(), colors.end())+1;
 
#ifdef E_TIMER
  std::cout << "toal number of elements : " << colors.size() << std::endl;
#endif
  
#ifdef DEBUG_BOOLEAN
{
  Vector_t<bool> flag(colors.size());  
  for (size_t c=0; c < colmax; ++c)
  {
    for (size_t i=0; i < colors.size(); ++i)
      flag[i]=(colors[i]==c);
    
    std::ostringstream o;
    o << "zc_" << c ;
    
    NGON_DBG::write(o.str().c_str(), ACoordinate_t(_coord), _ngXs, &flag);
  }
}
#endif
  
  if (colmax == 1) // special case : when a single zone unset : correspond to 2 meshes of a same domain
  {
    bool zone_is_unset=true;
    size_t sz = _zones.size();
    //size_t sz1 = _ngXs.PHs.size();
    for (size_t i=0; (i < sz) && zone_is_unset; ++i)
      zone_is_unset = (_zones[i] == Z_NONE);
    
    if (zone_is_unset)
    {
      _zones.clear();
      _zones.resize(sz, Z_IN);
      return OK;
    }
  }
  
#ifdef DEBUG_BOOLEAN
#ifdef DEBUG_ALG
/*
{
  std::vector<bool> keep(_ngXs.PHs.size(), false);
  std::cout << K_CONNECT::EltAlgo<K_MESH::Triangle>::_pool.size() << std::endl;
  for (std::set<E_Int>::const_iterator i=K_CONNECT::EltAlgo<K_MESH::Triangle>::_pool.begin(); i!=K_CONNECT::EltAlgo<K_MESH::Triangle>::_pool.end(); ++i)
  {if (*i == 1213)keep[*i]=true;std::cout << *i << std::endl;}
  NGON_DBG::write("toto", ACoordinate_t(_coord), _ngXs, &keep);
}
*/
#endif
#endif

  
#ifdef DEBUG_BOOLEAN
  std::vector<bool> keep(_ngXs.PHs.size(), false);
#endif
  
  //Interpret colors as zones
  assert (colors.size() >= _zones.size());
  Vector_t<eZone> col_to_z(colmax, Z_NONE);
  for (size_t i=0; i < _zones.size(); ++i)
  {
    const eZone& zi = _zones[i];
    
    if (zi == Z_NONE)
      continue;
    
    if (_Op == DIFF && zi != Z_1) //fixme : hack filter to boost : no necessary to retriev all the zones when diffsurf. similar thing should be done for all the other cases.
      continue;
    
    const E_Int& ci = colors[i];
    //eZone& ctzi = col_to_z[ci];

    if (col_to_z[ci] != Z_NONE)
    {
#ifdef DEBUG_BOOLEAN
      if (col_to_z[ci] != zi)
      {
        //std::cout << "i/zi/ci/col_to_z : " << i << "/" << zi << "/" << ci << "/" <<  col_to_z[ci] << std::endl;
        keep[i]=true;
      }
#else
      assert (col_to_z[ci] == zi);
#endif
      continue;
    }

    col_to_z[ci]=zi;
  }
  
#ifdef DEBUG_BOOLEAN
  //NGON_DBG::write("wrong_col", ACoordinate_t(_coord), _ngXs, &keep);
#endif
  
  // Now update missing zones (Z_NONE)
  _zones.resize(_ngXs.PHs.size(), Z_NONE);
  for (size_t i=0; i < _zones.size(); ++i)
  {
    if (_zones[i] == Z_NONE)
      _zones[i]=col_to_z[colors[i]];
    
#ifdef DEBUG_BOOLEAN
    /*if (_XPol != SURFACE_RIGHT || _Op != DIFF) //fixme : due to previous filter
        assert(_zones[i] != Z_NONE);*/
#endif
  }
  
  // Do the split on ngX
  Vector_t<bool> flag(colors.size());
  ngon_type ngcpy(_ngXs), *ngi(0);
  _ngXs.clear(); // fixme !! ancestor are lost !!!!
  //WARNING : _ng1 and/or _ng2 must have been cleared properly before the following. 
  for (E_Int i=0; i < colmax; ++i)
  {
    eZone Zi = col_to_z[i];
    if (_XPol == SURFACE_RIGHT && _Op == DIFF && Zi != Z_1)
      continue;
    
    for (size_t j=0; j < colors.size(); ++j) flag[j]=(_zones[j]==Zi);
    
    switch (Zi)
    {
      case Z_1 : ngi=&_ng1; break;
      case Z_2 : ngi=&_ng2; break;
      case Z_IN: ngi=&_ngXs; break;
      default : ngi=NULL; break;
    }
    
    if (!ngi)
      continue;
    
    Vector_t<E_Int> npgids;
    ngcpy.select_phs(ngcpy, flag, npgids, *ngi);
 
#ifdef DEBUG_BOOLEAN
    if (Zi== Z_1)
      NGON_DBG::write("first", ACoordinate_t(_coord), *ngi);
    else if (Zi == Z_2)
      NGON_DBG::write("second", ACoordinate_t(_coord), *ngi);
    else if (Zi == Z_IN)
      NGON_DBG::write("in", ACoordinate_t(_coord), *ngi);

#endif
  }
  
#ifdef FLAG_STEP
  std::cout << "NGON Boolean : __classify_soft for " << _ngXs.PHs.size() << " elements : " << c.elapsed() << std::endl;
#endif

  return err;
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__remove_orientations(std::map<E_Int, Vector_t<E_Int> >& PHT3s, E_Int shift)
{
  _normals.resize(3, _normals.cols()/2); //remove duplicates
  
  std::map<E_Int, Vector_t<E_Int> >::iterator it(PHT3s.begin()), itEnd(PHT3s.end());
  for (; it!=itEnd; ++it)
  {
    Vector_t<E_Int>& pht3i = it->second;
    E_Int sz = pht3i.size();
    // Remove orientation info.
    for (E_Int i = 0; i < sz; ++i)
      pht3i[i] = (pht3i[i] < shift) ? pht3i[i] : pht3i[i]-shift;
  }
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__build_neighbors_table
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const K_FLD::IntArray& connectT3o, K_FLD::IntArray& neighbors)
{
  neighbors.clear();
  neighbors.resize(3, connectT3o.cols(), E_IDX_NONE);
  
  // GET THE MAP EDGE to T3s (non oriented)
  K_FLD::ArrayAccessor<K_FLD::IntArray> acT3(connectT3);
  typedef K_CONNECT::EltAlgo<K_MESH::Triangle> algoT3;
  algoT3::BoundToEltType noE_to_oTs; // non oriented edge to oriented triangles (1 to n).
  
  //E_Int maxN = 
  algoT3::getBoundToElements(acT3, noE_to_oTs);
  //std::cout << maxN << std::endl; 
  
  E_Int K0, E0, E1, N1, sz, Kb, i, j, shift(connectT3.cols());
  //E_Float q, qb; 
  algoT3::BoundToEltType::iterator it, itEnd(noE_to_oTs.end());
  K_FLD::IntArray::const_iterator pS;
  for (it = noE_to_oTs.begin(); it != itEnd; ++it)
  {
    const K_MESH::Triangle::boundary_type& E = it->first;
    Vector_t<E_Int>& T3s = it->second;
        
    K0 = T3s[0];
    E0 = E.node(0);
    E1 = E.node(1); 
    
#ifdef DEBUG_BOOLEAN
    /*if (zeEdge)
    {
      std::cout << "dummy for hooking a breakpoint" << std::endl;
    }*/
#endif
    
    pS = connectT3o.col(K0);
    i = K_MESH::Triangle::getLocalNodeId(pS, E0);
    N1 = *(pS+(i+1)%3);
    // Finding out E0, E1
    if (N1 != E1)
      std::swap(E0, E1); // real orientation of E in K0 is the opposite of the non oriented shared edge.
    
    sz = T3s.size();
    if (sz < 2) //baffle
    {
      assert (sz == 1);
      Kb=K0+shift;
      //pSb = connectT3o.col(Kb);
      i = K_MESH::Triangle::getLocalNodeId(pS, E0);
      j = (i==0) ? 0 : (i==1) ? 2 : 1;//K_MESH::Triangle::getLocalNodeId(pSb, E0);//optim : since they are opposed..
      
      neighbors(i, K0) = Kb;
      neighbors(j, Kb) = K0;
    }
    else
    {
      // Loop over the other Ki and sort : upon exit contains {K0b, k1, ..., kn} with kj = Kj or Kjb
      __sort_T3_sharing_an_edge(E0, E1, shift, _normals, coord, connectT3o, T3s);
      
      // Now update the neighbors table
      __update_neigboring(E0, E1, shift, connectT3o, T3s, neighbors);
    }
  }

#ifdef DEBUG_BOOLEAN
 /*std::cout << neighbors << std::endl;
  Vector_t<E_Int> colors(neighbors.cols(), 0);
  bool found = false;
  for (size_t i = 0; i < neighbors.rows(); ++i)
    for (size_t j = 0; j < neighbors.cols(); ++j)
      if (neighbors(i,j) == E_IDX_NONE)
      {
        std::cout << "(i,j) : " << "(" << i << "," << j << ")" << std::endl;
        found = true;
        colors[j]=1;
      }
  //std::cout << "found ? : " << found << std::endl;
  MIO::write("neihbors.mesh", coord, connectT3o, "TRI", 0, &colors);
  
  K_FLD::IntArray ctmp;
  E_Int N0 = 3549;
  ctmp.pushBack(connectT3o.col(N0), connectT3o.col(N0)+3);
  for (size_t i = 0; i < neighbors.rows(); ++i)
    ctmp.pushBack(connectT3o.col(neighbors(i, N0)), connectT3o.col(neighbors(i, N0))+3);
  MIO::write("Neigi.mesh", coord, ctmp, "TRI");*/
#endif
  
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__assemble_PHT3s
(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3o, const K_FLD::IntArray& neighbors, std::map<E_Int, Vector_t<E_Int> >& PHT3s)
{
  Vector_t<E_Int> T3_to_PHT3;
  K_CONNECT::EltAlgo<K_MESH::Triangle>::coloring(neighbors, T3_to_PHT3);
  
#ifdef DEBUG_BOOLEAN
  //MIO::write("PHT3.mesh", coord, connectT3o, "TRI", 0, &T3_to_PHT3);
#endif
  
  PHT3s.clear();
  E_Int sz(connectT3o.getSize());
  for (E_Int i = 0; i < sz; ++i)
    PHT3s[T3_to_PHT3[i]].push_back(i);
  
#ifdef DEBUG_BOOLEAN
  //NGON_DBG::draw_PHT3(coord, connectT3o, PHT3s, T3_to_PHT3[1454]);
  //NGON_DBG::draw_PHT3(coord, connectT3o, PHT3s, T3_to_PHT3[1454+4712]);
#endif
  
  return 0;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__remove_parasite_PHT3s
(const Vector_t<E_Int>& is_skin, std::map<E_Int, Vector_t<E_Int> >& PHT3s, const K_FLD::IntArray& connectT3o
)
{
  // WARNING : this algorithm assume the connexity of PH set.
  E_Int shift(connectT3o.cols()/2), oTi;
  std::map<E_Int, Vector_t<E_Int> >::iterator it(PHT3s.begin()), itEnd(PHT3s.end());
  
#ifdef DEBUG_BOOLEAN
  std::vector<E_Int> faulty_PHT3, faulty_colors;
  //get the biggest one (supposed to be the envelop parasite)
  E_Int PHi = -1, maxT = -1, nb_T3s;
  for (std::map<E_Int, Vector_t<E_Int> >::const_iterator it1 = PHT3s.begin(); it1 != PHT3s.end(); ++it1)
  {
    nb_T3s = it1->second.size();
    if (maxT < nb_T3s)
    {
      maxT = it1->second.size();
      PHi = it1->first;
    }
  }
  NGON_DBG::draw_PHT3(_coord, connectT3o, PHT3s, PHi);
#endif
    
  _tmp_vec.clear();
  Vector_t<bool> treat;
  
  bool pure_skin;
  for (; it!=itEnd; ++it)
  {
    Vector_t<E_Int>& pht3i = it->second;
    E_Int sz = pht3i.size();
    //remove baffles (they disturb the logic)
    _tmp_set_int.clear();
    treat.clear();
    treat.resize(sz);
    for (E_Int i = 0; i < sz; ++i)
    {
      const E_Int& Ti = pht3i[i];
      oTi = (Ti < shift) ? Ti : Ti-shift;
      if (!_tmp_set_int.insert(oTi).second)
        _tmp_set_int.erase(oTi);
    }
    for (E_Int i = 0; i < sz; ++i)
    {
      const E_Int& Ti = pht3i[i];
      oTi = (Ti < shift) ? Ti : Ti-shift;
      treat[i] = (_tmp_set_int.find(oTi) != _tmp_set_int.end());
    }
    
#ifdef DEBUG_BOOLEAN
    faulty_PHT3.clear();
#endif
    
    pure_skin = true;
    for (E_Int i = 0; (i < sz)
            
#ifndef DEBUG_BOOLEAN
            && pure_skin
#endif
            ; ++i)
    {
      if (!treat[i])
        continue;

      const E_Int& Ti = pht3i[i];
      //const E_Int& isk = is_skin[(Ti < shift) ? Ti : Ti - shift];
      pure_skin &= (Ti < shift) ? (is_skin[Ti] != 0) : (is_skin[Ti-shift] == 3) || (is_skin[Ti-shift] < 0); // only skin ans if not original must be 3 (for inital skin, must be oriented toward exterior)
      
#ifdef DEBUG_BOOLEAN
      bool b = (Ti < shift) ? (is_skin[Ti] != 0) : (is_skin[Ti-shift] == 3) || (is_skin[Ti-shift] < 0);
      if (PHi == it->first)
      { E_Int Ci;
        if (!b)
        {
          Ci = is_skin[Ti < shift ? Ti : Ti - shift];
          std::cout << "Ti : " << Ti << " skin val : " << Ci << std::endl;
        }
        else //regular
          Ci=10;
        
        faulty_PHT3.push_back(Ti);
        faulty_colors.push_back(Ci);
      }
#endif
    }
    
    if (pure_skin) //parasite
      _tmp_vec.push_back(it->first);
#ifdef DEBUG_BOOLEAN
    else
    {
      //
      bool is_closed, is_manifold;
      //separate analysis by zone
      bool healthy = TRI_DBG::analyze_T3_set(PHi, _coord, connectT3o, faulty_PHT3, is_manifold, is_closed, &faulty_colors);
      if (!healthy)
        std::cout << "PHi " << PHi << " is not healthy." << std::endl;
      // analyze the PHT3 (except no treat parts) as a whole
      //healthy = TRI_DBG::analyze_T3_set(PHi, _coord, connectT3o, pht3i, is_manifold, is_closed);*/
    }
#endif
  }
 
  if (_tmp_vec.empty())
    return -1;
  
  for (size_t i=0; i < _tmp_vec.size(); ++i)
  {
#ifdef DEBUG_BOOLEAN
    //NGON_DBG::draw_PHT3(_coord, connectT3o, PHT3s, _tmp_vec[i]);
#endif
    PHT3s.erase(_tmp_vec[i]);
  }
  
  return -2;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__check_PHT3s_closure
(bool and_manifoldness, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3o, std::map<E_Int, Vector_t<E_Int> >& PHT3s)
{
  typedef std::map<E_Int, Vector_t<E_Int> > PHT3s_t;
  typedef PHT3s_t::iterator iterator;
  
  Vector_t<E_Int> unclosed_PHT3s;
  const E_Int* pS;
  std::map<K_MESH::NO_Edge, E_Int> emap;
  std::map<K_MESH::NO_Edge, E_Int>::iterator itm;
  E_Int sz;
  K_MESH::NO_Edge noE;
  
  for (iterator i = PHT3s.begin(); (i != PHT3s.end()) ; ++i)
  {
    emap.clear();
    const E_Int& PHi = i->first;
    Vector_t<E_Int>& T3s = i->second;
    sz = T3s.size();
    if (sz < 4) //cannot be a closed PHT3s shape, so 
    {
#ifdef DEBUG_BOOLEAN
      unclosed_PHT3s.push_back(PHi); // store it in debug
      continue;
#else
      return PHi; // return an error in release
#endif
    }
    
    for (E_Int j = 0; j < sz; ++j)
    {
      const E_Int& Tj=T3s[j];
      pS = connectT3o.col(Tj);
            
#ifdef DEBUG_BOOLEAN
      /*const E_Int & n0=connectT3o(0,Tj);
      const E_Int & n1=connectT3o(1,Tj);
      const E_Int & n2=connectT3o(2,Tj);*/
#endif
      E_Float h;
      for (size_t n=0; n < 3; ++n)
      {
        noE.setNodes(*(pS+n), *(pS+(n+1)%3));
        //h = noE.hash_func();
        itm = emap.find(noE);
        if (itm == emap.end())
          emap[noE]=1;
        else
          ++(itm->second);
      }
    }

    bool ok=true;
    for (itm = emap.begin(); itm != emap.end(); ++itm)
    {
      //const K_MESH::NO_Edge& e = itm->first; 
      if ( (itm->second == 1) || (and_manifoldness && itm->second != 2) ){
        ok=false; break;
      }
    }
    if (!ok) // Error 
    {
#ifdef DEBUG_BOOLEAN
    
      Vector_t<E_Int> colors;
      colors.resize(T3s.size(), 0);
      
      for (E_Int j = 0; j < sz; ++j)
      {
        const E_Int& Tj=T3s[j];
        pS = connectT3o.col(Tj);
            
        /*const E_Int & n0=connectT3o(0,Tj);
        const E_Int & n1=connectT3o(1,Tj);
        const E_Int & n2=connectT3o(2,Tj);*/
        
        E_Float h;
        for (size_t n=0; n < 3; ++n)
        {
          noE.setNodes(*(pS+n), *(pS+(n+1)%3));
          itm = emap.find(noE);
          if ( (itm->second == 1) || (and_manifoldness && itm->second != 2) )
            colors[j]=1;
        }
      }

      K_FLD::IntArray tmp;
      tmp.append_selection(connectT3o, T3s);
      MIO::write("unclosed_PH.mesh", coord, tmp, "TRI", 0, &colors);
      
      unclosed_PHT3s.push_back(PHi);
#else
      return PHi;
#endif
    }
  }

  sz = 0;
#ifdef DEBUG_BOOLEAN  
  sz = unclosed_PHT3s.size();
  for (size_t ii = 0; ii < sz; ++ii) 
  {
    NGON_DBG::draw_PHT3(coord, connectT3o, PHT3s, unclosed_PHT3s[ii]);
  }
#endif

  return (sz == 0) ? -1/*OK*/ : unclosed_PHT3s[0];
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__assemble_PGT3s
(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& T3_to_PG, std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s)
{
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif  

  E_Int err(0);
  for (std::map<E_Int, Vector_t<E_Int> >::const_iterator i = PHT3s.begin(); i != PHT3s.end(); ++i)
  {
    const E_Int& PHi = i->first;
    const Vector_t<E_Int>& T3i = i->second;
    for (size_t j = 0; j < T3i.size(); ++j)
    {
      PH_to_PGT3s[PHi][T3_to_PG[T3i[j]]].push_back(T3i[j]);
    }
  }
  
#ifdef FLAG_STEP
  if (chrono::verbose>0) std::cout << "NGON Boolean : __build_PHs : __assemble_PGT3s : " << c.elapsed() << std::endl;
#endif
   
  return err;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__split_non_connex_PGT3s
(const K_FLD::IntArray& connectT3, std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s)
{
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif  
  E_Int err(0), colmax, maxPGId;
  std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >::const_iterator it;
  std::map<E_Int, Vector_t<E_Int> >::const_iterator itPG, ritPG;
  K_FLD::IntArray cT3, neighbors;
  K_FLD::IntArray::const_iterator pS;
  Vector_t<E_Int> colors, T3scopy;
  size_t szT3;
  
  for (it = PH_to_PGT3s.begin(); it != PH_to_PGT3s.end(); ++it)
  {
    const E_Int& PHi = it->first;
#ifdef DEBUG_BOOLEAN
    //std::cout << "PHi : " << PHi << std::endl;
#endif

    const std::map<E_Int, Vector_t<E_Int> >& PGT3s = it->second;
    maxPGId = PGT3s.rbegin()->first;

    for (itPG = PGT3s.begin(); itPG != PGT3s.end(); ++itPG)
    {
      const E_Int& PGi = itPG->first;
      const Vector_t<E_Int> & T3s = itPG->second;
      
      // Put the triangles in an IntArray
      cT3.clear();
      szT3 = T3s.size();
      
      if (T3s.size() == 1)
        continue;
      
      for (size_t i=0; i < szT3; ++i)
      {
        pS = connectT3.col(T3s[i]);
        cT3.pushBack(pS, pS+3);
      }
      
      // find the neighboring matrix
      K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(cT3, neighbors);
      // do the coloring to see non connexity
      K_CONNECT::EltAlgo<K_MESH::Triangle>::coloring_pure (neighbors, colors);
      //
      colmax = *std::max_element(colors.begin(), colors.end())+1;
      if (colmax==1) //connex
        continue;
      
#ifdef DEBUG_BOOLEAN
        //MIO::write("nonconnex.mesh", _coord, cT3, "TRI");
#endif
      
      //Split
      T3scopy = T3s;
      PH_to_PGT3s[PHi][PGi].clear();
      
      for (size_t i=0; i < szT3; ++i)
        PH_to_PGT3s[PHi][PGi + colors[i] * maxPGId].push_back(T3scopy[i]); //use the initial entry for color 0, then append with a unique id
    }
  }
  
#ifdef FLAG_STEP
  if (chrono::verbose>0) std::cout << "NGON Boolean : __build_PHs : __split_non_connex_PGT3s : " << c.elapsed() << std::endl;
#endif
  
  return err;
}

#define OPP_PH(neighbors, Ti, PHi) (neighbors(0,Ti) == PHi ? neighbors(1,Ti) : neighbors(0,Ti))

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__split_multiply_shared_PGT3s
(const K_FLD::IntArray& connectT3, std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s)
{
  E_Int err(0);
  
#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif  
  
  // Build info T3->PH
  K_FLD::IntArray neighbors;
  __build_T3_to_PH(connectT3, PH_to_PGT3s, neighbors);
  
  //
  std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >::const_iterator it;
  std::map<E_Int, Vector_t<E_Int> >::const_iterator itPG;
  size_t szT3;
  E_Int PHcur(E_IDX_NONE), PHopp, maxPGId;
  bool same_PH;
  Vector_t<E_Int> T3scopy;
  std::map<E_Int, E_Int> for_unique_id;
  std::map<E_Int, E_Int>::iterator itUID;
          
  for (it = PH_to_PGT3s.begin(); it != PH_to_PGT3s.end(); ++it)
  {
    const E_Int& PHi = it->first;
    const std::map<E_Int, Vector_t<E_Int> >& PGT3s = it->second;
    maxPGId = PGT3s.rbegin()->first;
    for_unique_id.clear();
    E_Int uid = 0;
    
    for (itPG = PGT3s.begin(); itPG != PGT3s.end(); ++itPG)
    {
      const E_Int& PGi = itPG->first;
      const Vector_t<E_Int> & T3s = itPG->second;
      
      assert (!T3s.empty());
      szT3 = T3s.size();
      
      if (szT3 == 1)
        continue;
      
      same_PH=true;
      PHcur=OPP_PH(neighbors, T3s[0], PHi);
      for (size_t i=1; (i < szT3) && same_PH; ++i)
        same_PH=(PHcur == OPP_PH(neighbors, T3s[i], PHi));
      
      if (same_PH)
        continue;

      uid=0;
      for_unique_id.clear();
      for_unique_id[PHcur] = uid++; //0 for PHcur : use the initial entry for this color
      
      // Need a split
      T3scopy = T3s;
      PH_to_PGT3s[PHi][PGi].clear();
      
      for (size_t i=0; i < szT3; ++i)
      {
        PHopp = OPP_PH(neighbors, T3scopy[i], PHi);

        itUID = for_unique_id.find(PHopp);
        if (itUID == for_unique_id.end())
            for_unique_id[PHopp] = uid++;

        E_Int ukey = PGi + for_unique_id[PHopp] * maxPGId; //use the initial entry for color 0, then append with a unique id
        PH_to_PGT3s[PHi][ukey].push_back(T3scopy[i]);
      }
    }
  }
  
#ifdef FLAG_STEP
  if (chrono::verbose>0) std::cout << "NGON Boolean : __build_PHs : __split_multiply_shared_PGT3s : " << c.elapsed() << std::endl;
#endif
  
  return err;
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__build_T3_to_PH
(const K_FLD::IntArray& connectT3, std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >& PH_to_PGT3s, K_FLD::IntArray& neighbors)
{
  neighbors.clear();
  neighbors.resize(2, connectT3.cols(), E_IDX_NONE);
  
  std::map<E_Int, std::map<E_Int, Vector_t<E_Int> > >::const_iterator it;
  std::map<E_Int, Vector_t<E_Int> >::const_iterator itPG;
  size_t szT3;
          
  for (it = PH_to_PGT3s.begin(); it != PH_to_PGT3s.end(); ++it)
  {
    const E_Int& PHi = it->first;
    const std::map<E_Int, Vector_t<E_Int> >& PGT3s = it->second;
    
    for (itPG = PGT3s.begin(); itPG != PGT3s.end(); ++itPG)
    {
      //const E_Int& PGi = itPG->first;
      const Vector_t<E_Int> & T3s = itPG->second;
      
      szT3 = T3s.size();
      for (size_t i=0; i < szT3; ++i)
      {
        const E_Int& Ti = T3s[i];
        if (neighbors(0, Ti) == E_IDX_NONE)
          neighbors(0, Ti) = PHi;
        else
          neighbors(1, Ti) = PHi;
      }
    }
  }
}

///
template <>
template <> inline
E_Int NGON_BOOLEAN_FLD_SPECIALIZED::__aggregate<NGON_BOOLEAN_FLD_SPECIALIZED::NONE>
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs)
{
  // Just pass the T3s to the NGON format.
  __aggregate_none(connectT3, indices, agg_pgs);
  return 0;
}

///
template <>
template <> inline 
E_Int NGON_BOOLEAN_DYN_SPECIALIZED::__aggregate<NGON_BOOLEAN_DYN_SPECIALIZED::NONE>
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs)
{
  // Just pass the T3s to the NGON format.
  __aggregate_none(connectT3, indices, agg_pgs);
  return 0;
}

TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__aggregate_none
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs)
{
  agg_pgs.clear();
  
  Vector_t<E_Int>& molecule = _tmp_vec;

  molecule.resize(4);
  molecule[0] = 3;
  
  const E_Int* pt;
  for (size_t i = 0; i < indices.size(); ++i)
  {
    pt = connectT3.col(indices[i]);
    
    molecule[1] =1+ *(pt++); //ngon starts at 1.
    molecule[2]=1+*(pt++);   //ngon starts at 1.
    molecule[3]=(*pt)+1;     //ngon starts at 1.
    
    agg_pgs.add(molecule);
  }
  
  return 0;
}

///
template <>
template <> inline 
E_Int NGON_BOOLEAN_DYN_SPECIALIZED::__aggregate<NGON_BOOLEAN_DYN_SPECIALIZED::FULL>
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs)
{
  E_Int err = 0;
  if (indices.size() > 1)
  {
    err = __aggregate_full(connectT3, indices, agg_pgs);
    if (err) // if only one or it could not do it, pass the T3s
      err = __aggregate<NONE>(connectT3, indices, agg_pgs);
  }
  else
    err = __aggregate<NONE>(connectT3, indices, agg_pgs);
  return err;
}

TEMPLATE_COORD_CONNECT
E_Int  NGON_BOOLEAN_CLASS::__aggregate_full
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs)
{
  //WARNING : ASSUME HERE TO HAVE A UNIQUE CONTOUR, THIS RELIES ON __split_non_connex_PGT3s
  agg_pgs.clear();
  
  std::deque<E_Int>& PGi = _tmp_deq;
  E_Int err = K_CONNECT::MeshTool::get_polygonal_boundary(connectT3, indices, PGi, _tmp_set_oedge, _tmp_map_nn);
  if (err)
    return 1;

  K_CONNECT::IdTool::shift(PGi, 1);
  
  size_t sz = PGi.size();
  PGi.push_front(sz);
  agg_pgs.add(PGi);
  return 0;
}

///
template <>
template <> inline
E_Int NGON_BOOLEAN_FLD_SPECIALIZED::__aggregate<NGON_BOOLEAN_FLD_SPECIALIZED::CONVEX>
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs)
{
  E_Int err = __aggregate_convex(connectT3, indices, agg_pgs);
  if (err) // if could not do it pass the T3s
    err = __aggregate<NONE>(connectT3, indices, agg_pgs);
  return err;
}

///
template <>
template <> inline 
E_Int NGON_BOOLEAN_DYN_SPECIALIZED::__aggregate<NGON_BOOLEAN_DYN_SPECIALIZED::CONVEX>
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs)
{
  E_Int err = __aggregate_convex(connectT3, indices, agg_pgs);
   if (err) // if could not do it pass the T3s
    err = __aggregate<NONE>(connectT3, indices, agg_pgs);
  return err;
}

TEMPLATE_COORD_CONNECT
E_Int  NGON_BOOLEAN_CLASS::__aggregate_convex
(const K_FLD::IntArray& connectT3, const Vector_t<E_Int>& indices, ngon_unit& agg_pgs)
{
  agg_pgs.clear();
  
  // WARNING : ASSUME THAT THE INPUT MESH IS CONNEX
  K_FLD::IntArray& neighbors = _tmp_IntArray;
  std::vector<K_FLD::IntArray*> pool, trash;
  
  // Init
  E_Int sz(indices.size());
  K_FLD::IntArray* connectM = new K_FLD::IntArray;
  connectM->reserve(3, sz);
  for (E_Int i = 0; i < sz; ++i)
    connectM->pushBack(connectT3.col(indices[i]), connectT3.col(indices[i])+3);
  
#ifdef DEBUG_BOOLEAN
    //enabling_write("cM.mesh", _coord, *connectM, "TRI");
#endif
  
  pool.push_back(connectM);
  trash.push_back(connectM);//to delete at the end
 
  std::map< E_Int, std::pair<E_Int, E_Int> > node_to_nodes_m;
  K_CONT_DEF::int_pair_vector_type boundaries;
  K_FLD::IntArray connectB;

#ifdef DEBUG_BOOLEAN
  E_Int count=0;
  if (_enabled)
    count=0;
#endif
  E_Int err=0;
  while (!pool.empty() && !err)
  {
    K_FLD::IntArray* ci = pool.back(); pool.pop_back();
    K_FLD::IntArray& connectMi = *ci;
    
    if (connectMi.cols() == 0)
      continue;
    
#ifdef DEBUG_BOOLEAN
    std::ostringstream o;
    o << "cM" << count++ << ".mesh";
    if (_enabled)
      MIO::write(o.str().c_str(), _coord, connectMi, "TRI");
#endif
        
    //Build the neighbourhood matrix.
    K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectMi, neighbors);

#ifdef DEBUG_BOOLEAN  
    if (_enabled)
      std::cout << neighbors << std::endl;
#endif
    
    // Get the contour.
    K_CONNECT::MeshTool::getBoundaryT3Mesh(connectMi, neighbors, boundaries);
    connectB.clear();
    E_Int sz = boundaries.size();
    
    if (sz == 0)
      return 1;
    
    K_MESH::Edge E;
    for (E_Int i = 0; i < sz; ++i)
    {
      const E_Int& K = boundaries[i].first;
      const E_Int& n = boundaries[i].second; 
      E.setNodes(connectMi((n+1)%3,K), connectMi((n+2)%3, K));
      connectB.pushBack(E.begin(), E.begin()+2);
    }
    
#ifdef DEBUG_BOOLEAN  
    if (_enabled)
      std::cout << connectB << std::endl;
#endif
            
    // Sort the nodes
    err = BARSplitter::get_node_to_nodes(connectB, node_to_nodes_m);
    if (err)
      break;
    
    // get worst concavity
    E_Int K0, n0, Eip1; //returned K0 is local
    std::deque<E_Int>& PGi = _tmp_deq;
    bool convex = __is_convex(boundaries, node_to_nodes_m, _coord, connectMi, _normals, indices, PGi, K0, n0, Eip1);
    
    if (convex)
    {
      PGi.push_front(PGi.size());
      agg_pgs.add(PGi);
    }
    else // do the split starting at Nworst
    {
      K_FLD::IntArray *c1(new K_FLD::IntArray), *c2(new K_FLD::IntArray);
      
      err = __cut_mesh_convex_line(connectMi, _normals, indices, connectB, K0, n0, Eip1, neighbors, *c1, *c2);
      if (!err)
      {
        pool.push_back(c1); pool.push_back(c2);
      }
      trash.push_back(c1);
      trash.push_back(c2);
    }
  }
  
  //cleaning
  for (size_t i = 0; i < trash.size(); ++i)
    if (trash[i])delete trash[i];
  
  return err;
}

///
TEMPLATE_COORD_CONNECT
E_Int NGON_BOOLEAN_CLASS::__cut_mesh_convex_line
(const K_FLD::IntArray& connectT3,
 const K_FLD::FloatArray& normals, const Vector_t<E_Int>& globalId,
 const K_FLD::IntArray& connectB,
 E_Int K0, E_Int n0, E_Int B, K_FLD::IntArray& neighbors, 
 K_FLD::IntArray& connectT31, K_FLD::IntArray& connectT32)
{
  assert (K0 != E_IDX_NONE);
    
  std::set<E_Int>& bnodes = _tmp_set_int;
  bnodes.clear();
  for (E_Int i = 0; i < connectB.cols(); ++i) bnodes.insert(connectB(0,i));
  
  E_Float Eip1Aip1[3], Ek[3], Dir[3], ps;
  E_Int Ai=0, Aip1, Ki(K0), Kip1, Kim1(K0), ni((n0+1)%3), nim1=0, Ei(connectT3((n0+1)%3, K0)), Eip1(connectT3((n0+2)%3, K0));
  K_FLD::IntArray::const_iterator pKi;
  
  //K_FUNC::diff<3>(_coord.col(Eip1), _coord.col(Ei), &Dir1[0]);
  //K_FUNC::diff<3>(_coord.col(Eip1), _coord.col(B), &Dir2[0]);
  //K_FUNC::sum<3>(0.5, &Dir1[0], 0.5, &Dir2[0], &Dir[0]);//ax+by
  //K_FUNC::normalize<3>(Dir);
  
  //
  while (1)
  {
    K_FUNC::diff<3>(_coord.col(Eip1), _coord.col(Ei), &Dir[0]);
    K_FUNC::normalize<3>(Dir);
    
    while (1)
    {
      pKi = connectT3.col(Ki);
      Aip1 = *(pKi + (ni+2)%3);
      Kip1 = neighbors(ni, Ki);
      
      if (Kip1 == E_IDX_NONE) break;
      
      K_FUNC::diff<3>(_coord.col(Aip1), _coord.col(Eip1), &Eip1Aip1[0]);
      K_FUNC::normalize<3>(Eip1Aip1);
      K_FUNC::crossProduct<3>(Dir, Eip1Aip1, Ek);
      ps = K_FUNC::dot<3>(Ek, normals.col(globalId[Ki]));
      
      if (ps < -E_EPSILON) break;
       
      // current best
      Kim1 = Ki;
      Ai = Aip1;
      nim1 = ni;
      
      //nex step
      Ki = Kip1;
      pKi = connectT3.col(Ki);
      ni = K_MESH::Triangle::getLocalNodeId(pKi, Aip1);
    };
    
    // cut edge : (Eip1, Ai)
    
    //do the cut for this edge by updating neighbors.
    neighbors(nim1, Kim1) = E_IDX_NONE;
    neighbors((ni+2)%3, Ki) = E_IDX_NONE;
    /*if (Kip1 != E_IDX_NONE)
    {
      pKi = connectT3.col(Kip1);
      ni = K_MESH::Triangle::getLocalNodeId(pKi, Ai);
      neighbors((ni+2)%3, Kip1) = E_IDX_NONE;
    }*/
    
    if (bnodes.find(Ai) != bnodes.end()) break;
    
    Ei = Eip1;
    Eip1 = Ai;
    Ki = Kim1;
    ni = (nim1+1)%3; //new Ei
  };
  
  Vector_t<E_Int> colors;  
  K_CONNECT::EltAlgo<K_MESH::Triangle>::coloring(neighbors, colors);//hpc
  
  E_Int colmax = *std::max_element(colors.begin(), colors.end());
  if (colmax != 1)
  {
#ifdef DEBUG_BOOLEAN
    //MIO::write("part.plt", _coord, connectT3, "TRI");
#endif
    return 1;
  }
  
  for (size_t i = 0; i < colors.size(); ++i)
  {
    if (colors[i] == 0)
      connectT31.pushBack(connectT3.col(i), connectT3.col(i)+3);
    else
    {
      assert (colors[i]==1); //we must have 2 parts only.
      connectT32.pushBack(connectT3.col(i), connectT3.col(i)+3);
    }
  }

  return 0;
}

///
TEMPLATE_COORD_CONNECT
bool NGON_BOOLEAN_CLASS::__is_convex
(const K_CONT_DEF::int_pair_vector_type &boundaries,
 const std::map< E_Int, std::pair<E_Int, E_Int> >& node_to_nodes, 
 const K_FLD::FloatArray& coord, 
 const K_FLD::IntArray& connectT3, const K_FLD::FloatArray& normals, const Vector_t<E_Int>& globalId, 
       std::deque<E_Int>& sorted_contour, E_Int& K0, E_Int& n0, E_Int& Eip1)
{
  E_Int eim1, eip1, ei, K, n;
  E_Float worst_ps = 1., Ei[3], Ej[3], /*Ek[3],*/ ps;
  bool convex = true;
  Eip1=K0=n0=E_IDX_NONE;
  std::map< E_Int, std::pair<E_Int, E_Int> >::const_iterator itN;
  
  E_Float angle_max = K_CONST::E_PI* _convexity_tol; // a fraction betwen 0 and Pi
  E_Float cos_min = ::cos(angle_max); // cos is decreasing on [0; Pi]
  
  E_Float Z[3];
  for (size_t i = 0; i < boundaries.size(); ++i)
  {
    K = boundaries[i].first;
    n = boundaries[i].second;
    ei = connectT3((n+2)%3, K);
    eim1=connectT3((n+1)%3, K);
    itN = node_to_nodes.find(ei);
    assert (itN != node_to_nodes.end());
    eip1=itN->second.second;/*cT30((n+2)%3, K);//*/
    
    K_FUNC::diff<3>(coord.col(ei), coord.col(eim1), &Ei[0]);
    K_FUNC::diff<3>(coord.col(eip1), coord.col(ei), &Ej[0]);
    K_FUNC::sum<3>(normals.col(globalId[K]), coord.col(ei), Z);
    
    E_Float det = K_FUNC::zzdet4(coord.col(eim1), coord.col(ei), coord.col(eip1), Z);
    
    if (det >= 0.) continue; // convex
    
    K_FUNC::normalize<3>(Ei);
    K_FUNC::normalize<3>(Ej);

    E_Float c = K_FUNC::dot<3>(Ei, Ej);
    
    if (c < cos_min) // angle > anle max
    {
      convex = false;
  
      if ((K0 == E_IDX_NONE) || (c < worst_ps))
      {
        K0 = K; //returned K is local (for __cut_mesh_convex_line)
        n0 = n;
        worst_ps = c;
        Eip1 = eip1;
      }
    }
  }
  
  if (convex)
  {
    //sort the contour and format to ngon (starting at 1)
    sorted_contour.clear();
    K = boundaries[0].first;
    n = boundaries[0].second;
    ei = connectT3((n+2)%3, K);
    size_t sz = boundaries.size();
    for (size_t i = 0; i < sz; ++i)
    {
      sorted_contour.push_back(ei+1);
      itN = node_to_nodes.find(ei);
      assert (itN != node_to_nodes.end());
      ei = itN->second.second;//node_to_nodes[ei].second;
    }
  }
  
  return convex;
}


///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__swap
(K_FLD::IntArray& connectT3)
{
  E_Int sz = connectT3.getSize();
  for (E_Int i = 0; i < sz; ++i)
    std::swap(connectT3(1,i), connectT3(2,i));
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS:: __compact_and_join(ngon_type& ngio, K_FLD::FloatArray& coord)
{
  if (ngio.PGs.size() == 0)
  {
    //coord.clear();//fixme : to enable when Generator.py will have a proper handling of exceptions
    return;
  }
  K_FLD::ArrayAccessor<K_FLD::FloatArray> ca(coord);
  
  Vector_t<E_Int> nids;
  E_Int nb_merges = ::merge(ca, E_EPSILON, nids);
  
  if (nb_merges)
    ngio.PGs.change_indices(nids, true/*zero based*/); 
  
  ngon_type::compact_to_used_nodes(ngio.PGs, coord);
}

/// T3 -> PHT3
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__set_skin_PHT3s_zones
(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift, Vector_t<eZone>& zones
#ifdef DEBUG_BOOLEAN
, const K_FLD::IntArray& connectT3o
#endif
)
{   
  //

  zones.clear();
  zones.resize(PHT3s.size(),Z_NONE);//fixme hpc : only need one element per skin'connex parts...=> zones vector would reduce massively
  
  Vector_t<eZone> T3z(shift, Z_NONE);
  
  __flag_border_OUT_zones(PHT3s, is_skin, shift, Z_2, 1, zones, T3z);
  
#ifdef DEBUG_BOOLEAN
  typedef std::map<E_Int, Vector_t<E_Int> > cont_t;
  typedef cont_t::const_iterator const_it_t;
  const_it_t it;
{
  Vector_t<bool> keep(2*shift, false);
  E_Int i = 0;
  for (it = PHT3s.begin(); it != PHT3s.end(); ++it,  ++i)
  {
    if (zones[i] != Z_2) continue;
    for (size_t j=0; j < it->second.size(); ++j)
      keep[it->second[j]]=true;
  }
  
  MIO::write("Z_2.mesh", _coord, connectT3o, "TRI", &keep);
}
#endif
  
  __flag_border_OUT_zones(PHT3s, is_skin, shift, Z_1, 2, zones, T3z);
  
#ifdef DEBUG_BOOLEAN
{
  Vector_t<bool> keep(2*shift, false);
  E_Int i = 0;
  for (it = PHT3s.begin(); it != PHT3s.end(); ++it, ++i)
  {
    if (zones[i] != Z_1) continue;
    for (size_t j=0; j < it->second.size(); ++j)
    { E_Int id = it->second[j];
      keep[it->second[j]]=true;
    }
  }
  MIO::write("Z_1.mesh", _coord, connectT3o, "TRI", &keep);
}
{
  K_FLD::IntArray cT3;
  for (E_Int i = 0; i < shift; ++i)
    if (T3z[i] != Z_NONE)
    {
      cT3.pushBack(connectT3o.col(i), connectT3o.col(i) + 3);
    }
  MIO::write("T3z.mesh", _coord, cT3, "TRI");
}
#endif
  
  __flag_border_IN_zones(PHT3s, is_skin, shift, T3z, zones);
  
#ifdef DEBUG_BOOLEAN
{
  Vector_t<bool> keep(2*shift, false);
  E_Int i = 0;
  for (it = PHT3s.begin(); it != PHT3s.end(); ++it, ++i)
  {
    if (zones[i] != Z_IN) continue;
    for (size_t j=0; j < it->second.size(); ++j)
      keep[it->second[j]]=true;
  }
  MIO::write("Z_IN.mesh", _coord, connectT3o, "TRI", &keep);
}
#endif
  
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__flag_border_OUT_zones
(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift, eZone Z, E_Int col, Vector_t<eZone>& zones, Vector_t<eZone>& T3z)
{   
  //
  typedef std::map<E_Int, Vector_t<E_Int> > map_t;
  typedef Vector_t<E_Int> vec_t;
  map_t::const_iterator it;
  
  size_t i=0;
  for (it = PHT3s.begin(); it != PHT3s.end(); ++it, ++i)
  {
    if (zones[i] != Z_NONE)
      continue;
    
    const vec_t& T3s = it->second;
    //
    size_t sz= T3s.size();
    for (size_t j=0; (j< sz); ++j)
    {
      const E_Int& t = T3s[j];
      //we only consider elements having original skin triangles (before duplication) because their orientation is
      // toward the exterior (because of reorienting step) and therefore it means that the element is an outside one (Z_1 or Z_2).
      if ((t < shift) && (zABS(is_skin[t]) == col))
        T3z[t]=zones[i]=Z;
    }
  }
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__flag_border_IN_zones
(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift, const Vector_t<eZone>& T3z, Vector_t<eZone>& zones)
{   
  // Flag inside elemnts connected to Z_1 or Z_2 elements.
  //
  typedef std::map<E_Int, Vector_t<E_Int> > map_t;
  typedef Vector_t<E_Int> vec_t;
  map_t::const_iterator it;
  
  size_t i=0;
  // Z_IN are neighbors of Z having a skin face in common.
  for (it = PHT3s.begin(); it != PHT3s.end(); ++it, ++i)
  {
    if (zones[i] != Z_NONE)
      continue;
    
    const vec_t& T3s = it->second;
    size_t sz = T3s.size();

    for (size_t j=0; (j< sz); ++j)
    {
      const E_Int& t = T3s[j];
      if (t < shift)
        continue;
      if (// the opposite is an original flagged as "out" and it's attached to the corresponding skin
              ((T3z[t-shift] == Z_1) && (zABS(is_skin[t-shift]) ==2))
              ||
              ((T3z[t-shift] == Z_2) && (zABS(is_skin[t-shift]) ==1))
         )
              
      {
        zones[i] = Z_IN;
        break;
      }
    }
  }
}

///
TEMPLATE_COORD_CONNECT
void NGON_BOOLEAN_CLASS::__set_PH_history
(const std::map<E_Int, Vector_t<E_Int> >& PHT3s, const Vector_t<E_Int>& is_skin, E_Int shift, E_Int nb_pgs1, const K_FLD::IntArray& F2E, K_FLD::IntArray& anc_PH, bool soft)
{
  // Where does each PHT3 comes from ?
  //

  anc_PH.resize(2, PHT3s.size(), E_IDX_NONE);
  
  typedef std::map<E_Int, Vector_t<E_Int> > map_t;
  typedef Vector_t<E_Int> vec_t;
  map_t::const_iterator it;

  size_t i = 0;
  // Z_IN are neighbors of Z having a skin face in common.
  for (it = PHT3s.begin(); it != PHT3s.end(); ++it, ++i)
  {
    const vec_t& T3s = it->second;
    size_t sz = T3s.size();

#ifdef DEBUG_BOOLEAN
    E_Int hPH01(-1), hPH02(-1);
#endif

    for (size_t j = 0; (j< sz); ++j)
    {
      const E_Int& t = T3s[j];

      // set the historical PH
      E_Int T = (t<shift) ? t : t - shift;
      E_Int I = (t<shift) ? 1 : 0;
      E_Int wPG = _nT3_to_oPG[T];

      if ((_XPol == SOLID_RIGHT) && (wPG >= nb_pgs1)  && soft) // a soft element cut by hard polygon : we always keep the ouside part of it so no second ancestor in this case
        continue;

#ifdef DEBUG_BOOLEAN
      assert(!soft || (soft && wPG < F2E.cols()) ); // either hard (and then open layer bulkhead are not in F2E) or soft (and then the polygon must be in)
#endif
      
      if (!soft && wPG >= F2E.cols()) continue; //ignore open layer bulkheads 
      
      E_Int hPH = F2E(I, wPG);

      if (hPH == E_IDX_NONE)
        continue;

      anc_PH((wPG < nb_pgs1) ? 0 : 1, i) = hPH; //if nb_pgs1 < 0 it means that set_hitory is called for hard part which is always right operand

#ifdef DEBUG_BOOLEAN
      if ((wPG < nb_pgs1) && (hPH01 == -1))hPH01 = hPH;
      else if ((wPG >= nb_pgs1) && (hPH02 == -1))hPH02 = hPH;

      if (wPG < nb_pgs1)
      {
        assert(hPH01 == hPH);
        /*if (hPH01 != hPH)
        {
          std::cout << "stored (left operand): " << hPH01 << " and current : " << hPH << std::endl;
          std::cout << "other side is : " << F2E((I + 1) % 2, wPG) << std::endl;
        }*/
      }
      else
      {
        assert(hPH02 == hPH);
        /*if (hPH02 != hPH)
        {
          std::cout << "stored (right operand): " << hPH02 << " and current : " << hPH << std::endl;
          std::cout << "other side is : " << F2E((I+1)%2, wPG) << std::endl;
        }*/
      }
#else
      if ((anc_PH(0, i) != E_IDX_NONE) && (anc_PH(1, i) != E_IDX_NONE))
        break;
#endif
    }
  }
}

}

#endif	/* NGON_BOOLEANOPERATOR_H */

