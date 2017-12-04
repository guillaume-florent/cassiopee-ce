#ifndef NUGA_AGGLOMERATOR_H
#define NUGA_AGGLOMERATOR_H

#include "Fld/DynArray.h"
#include "Fld/ngon_t.hxx"
#ifdef DEBUG_AGGLOMERATOR
#include "NGON_debug.h"
#endif
/*#include <vector>
#include <set>
#include <map>
#include <deque>
#include "Def/DefTypes.h"*/


//#include  "Nuga/Boolean/NGON_debug.h"

#define NONEVAL -2
#define UNCHANGED -1

namespace NUGA
{

  class Agglomerator
  {

  public:
    typedef ngon_t<K_FLD::IntArray>                 ngon_type;
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    typedef std::set<K_MESH::NO_Edge>               eset_t;
    typedef std::set<E_Int>                         iset_t;
    typedef std::vector<E_Int>                      ivec_t;
    typedef std::map<E_Int, iset_t >                id_to_ids_t;
    typedef std::map<E_Int, E_Int>                  id_to_id_t;
    typedef std::vector<std::deque<E_Int> >         idlists_t;
    typedef std::map<K_MESH::NO_Edge, E_Float>      edge_angle_t;
    
  public: //WRAPPERS
   
    template<typename TriangulatorType>
    inline static void agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs);
    
//    template<typename TriangulatorType>
//    inline static void agglomerate_uncomputable_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo);
    
    template<typename TriangulatorType>
    inline static void collapse_uncomputable_pgs(K_FLD::FloatArray& crd, ngon_type& ngio);
    
    inline static E_Int collapse_pgs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids);
    inline static E_Int collapse_pgs2(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids);

  public:
    /// agglomerate superfluous polygons (multiply-shared by the same polyhedra. within the flatness tolerance only for area-computable polygons)
    inline static void simplify_phs
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
       ngon_type& ngo, ngon_unit& oriento, ngon_unit& phneighborso, const Vector_t<E_Int>* PHlist = 0);

    ///
//    template<typename TriangulatorType>
//    inline static void agglomerate_phs(const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& neighborsi, const Vector_t<E_Int>& PHlist, ngon_type& ngo);

    ///
    template<typename TriangulatorType>
    inline static void agglomerate_phs (const K_FLD::FloatArray& crd, 
                                        const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
                                        ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs);
  
  
  private:
      
    inline static void __simplify_phs
      (const K_FLD::FloatArray& crd, const ngon_type& ngi, E_Int PHi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
      ngon_unit& gagg_pgs, std::vector<E_Int>& nids, ngon_unit& wagg_pgs, std::map<E_Int, std::vector<E_Int> >& wneigh_to_faces);

  };

  ///
  void NUGA::Agglomerator::simplify_phs
  (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
   ngon_type& ngo, ngon_unit& oriento, ngon_unit& phneighborso, const Vector_t<E_Int>* PHlist)
  {
    typedef std::map<E_Int, std::vector<E_Int> > map_t;
    map_t neigh_to_faces;
    map_t::const_iterator it;

    std::vector<E_Int> nids;
    ngon_unit lagg_pgs, gagg_pgs;

    nids.clear();
    nids.resize(ngi.PGs.size(), NONEVAL);
    
    std::cout << "simplify_phs : initial nb of pgs : " << ngi.PGs.size() << std::endl;

    if (PHlist)
    {
      size_t nb_phs = PHlist->size();
      for (size_t i = 0; i < nb_phs; ++i)
      {
        E_Int PHi = (*PHlist)[i];
        NUGA::Agglomerator::__simplify_phs
          (crd, ngi, PHi, orienti, phneighborsi, angular_max, process_externals, gagg_pgs, nids, lagg_pgs, neigh_to_faces);
      }
    }
    else
    {
      E_Int nb_phs = ngi.PHs.size();
      for (E_Int PHi = 0; PHi < nb_phs; ++PHi)
      {
        NUGA::Agglomerator::__simplify_phs
          (crd, ngi, PHi, orienti, phneighborsi, angular_max, process_externals, gagg_pgs, nids, lagg_pgs, neigh_to_faces);
      }
    }
    
#ifdef DEBUG_AGGLOMERATOR
  {
    E_Int idmax=-1;
    for (size_t i=0; i < nids.size(); ++i) if (nids[i] != E_IDX_NONE) idmax = std::max(nids[i], idmax);
    assert (idmax == gagg_pgs.size()-1);
    K_FLD::IntArray cnto;
    ngon_type ng(gagg_pgs);
    ng.export_to_array(cnto);
    MIO::write("agg.plt", crd, cnto, "NGON");
  }
#endif
    
    ngo = ngi;
    
    //if (gagg_pgs.size() == 0)
      //return;

    // append new pgs
    E_Int nb_pgs_init = ngo.PGs.size();
    ngo.PGs.append(gagg_pgs);
    
#ifdef DEBUG_AGGLOMERATOR
  {
    K_FLD::IntArray cnto;
    ngon_type ng(ngo.PGs);
    ng.export_to_array(cnto);
    MIO::write("aggandinit.plt", crd, cnto, "NGON");
  }
#endif
   
    // sync nids
    size_t nsz = nids.size();
    for (size_t i=0; i < nsz;++i)
    {
      E_Int& ni = nids[i];
      ni = (ni == E_IDX_NONE) ? E_IDX_NONE : (ni <0/*UNCHANGED OR NONEVAL*/) ? i : ni + nb_pgs_init;
    }
    
    //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH("ph220.plt", crd, ngi, 220);

    // sync PHs
    E_Int nb_phs = ngo.PHs.size();
    std::vector<E_Int> dummy;
    ngo.PHs.remove_facets(nids, dummy);
    nb_phs = ngo.PHs.size();
    
    Vector_t<E_Int> pgnids, phnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);

#ifdef DEBUG_AGGLOMERATOR
    std::cout << "simplify_phs : final nb of pgs : " << ngo.PGs.size() << std::endl;
#endif

    // clean superfluous nodes
    ngon_type::simplify_pgs(ngo, crd);
  }

//  template<typename TriangulatorType>
//  void NUGA::Agglomerator::agglomerate_phs
//  (const K_FLD::FloatArray& crd, const ngon_type& ngi, const ngon_unit& neighborsi, const Vector_t<E_Int>& PHlist, ngon_type& ngo)
//  {
//    ngo.clear();
//    
//    // Fast return
//    if (PHlist.size() * ngi.PHs.size() * ngi.PGs.size() == 0) 
//    {
//      return;
//    }
//    
//    //std::cout << "ngi : initial nb of phs : " << ngi.PHs.size() << std::endl;
//    
//    std::vector<bool> frozen(ngi.PHs.size(), false); //process one agglo at a time
//    std::vector<E_Int> shared_pgs;
//
//    ngon_unit all_agg_phs, agg_phs;
//
//    E_Int nb_phs = PHlist.size();
//    for (E_Int ii = 0; ii < nb_phs; ++ii)
//    {
//      E_Int i = PHlist[ii];
//
//      if (frozen[i])
//        continue;
//
//      E_Int bestn = E_IDX_NONE;
//      size_t maxf = 0;
//
//      E_Int nb_neighs = ngi.PHs.stride(i);
//      const E_Int* neighs = neighborsi.get_facets_ptr(i);
//
//      // the best is the one sharing the most number of faces
//      for (E_Int n = 0; (n < nb_neighs); ++n)
//      {
//        E_Int j = *(neighs + n);
//        if (j == E_IDX_NONE)
//          continue;
//
//        if (frozen[j])
//          continue;
//
//        ngon_type::shared_faces_list(ngi, i, j, shared_pgs);
//        size_t nbf = shared_pgs.size();
//
//        if (nbf < maxf)
//          continue;
//
//        maxf = nbf;
//        bestn = n;
//      }
//
//      if (bestn == E_IDX_NONE) continue;
//
//      E_Int j = *(neighs + bestn);
//      frozen[i] = frozen[j] = true;
//
//      K_MESH::Polyhedron<UNKNOWN>::merge_two_phs(crd, ngi.PGs, 
//                                                 ngi.PHs.get_facets_ptr(i), ngi.PHs.stride(i),
//                                                 ngi.PHs.get_facets_ptr(j), ngi.PHs.stride(j), 
//                                                 agg_phs);
//
//      all_agg_phs.append(agg_phs);
//
//    }
//
//    ngo.PGs = ngi.PGs;
//    ngo.PHs = all_agg_phs;
//    
//    //std::cout << ngo.PHs.size()/2 << " uncomputable cells have been agglomerated." << std::endl;
//
//    // now add untouched ones
//    for (size_t i = 0; i < frozen.size(); ++i)
//    {
//      if (!frozen[i])
//      {
//        ngo.PHs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
//        if (!ngi.PHs._type.empty()) ngo.PHs._type.push_back(ngi.PHs._type[i]);
//      }
//    }
//    
//    //std::cout << "effective nb of phs in agglomerated ngo : " << ngo.PHs.size() << std::endl;
//
//    std::vector<E_Int> pgnids, phnids;
//    ngo.remove_unreferenced_pgs(pgnids, phnids);
//
//  }
  
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_phs
  (const K_FLD::FloatArray& crd, 
   const ngon_type& ngi, const ngon_unit& neighborsi, const ngon_unit& orienti, const Vector_t<E_Int>& PHlist,
   ngon_type& ngo, ngon_unit& oriento, E_Int& nb_aggs)
  {
    ngo.clear();
    oriento.clear();
    nb_aggs = 0;
    
    //std::cout << "ngi : initial nb of phs : " << ngi.PHs.size() << std::endl;
  
    // Fast return
    if (PHlist.size() * ngi.PHs.size() * ngi.PGs.size() == 0) return;
  
    TriangulatorType dt;
  
    std::vector<bool> frozen(ngi.PHs.size(), false); //process one agglo at a time
    std::vector<E_Int> shared_pgs;
    std::map<K_MESH::NO_Edge, E_Float> reflex_edges;
    std::set<K_MESH::NO_Edge> convex_edges;

    ngon_unit all_agg, cur_agg, all_ori, cur_ori, best_agg, best_ori;

    E_Int nb_phs = PHlist.size();
    for (E_Int ii = 0; ii < nb_phs; ++ii)
    {
      E_Int i = PHlist[ii];

      if (frozen[i])
        continue;

      E_Int bestn = E_IDX_NONE;
      size_t nbfmax=0;
      E_Float smax=0.;
      E_Float worst_reflex_angle=0.;
      E_Float qmax=0.;
      
      
      E_Int nb_neighs = ngi.PHs.stride(i);
      const E_Int* neighs = neighborsi.get_facets_ptr(i);
      const E_Int* pgsi = ngi.PHs.get_facets_ptr(i);

      // the best is the one sharing the most number of faces
      for (E_Int n = 0; (n < nb_neighs); ++n)
      {
        E_Int j = *(neighs + n);
        if (j == E_IDX_NONE)
          continue;

        if (frozen[j]) 
        {
          // continue; fixme : hack to have the targeted logic (but bad perfo) : 
          // a cell has to be aggregated with its best mate, i.e the one with most surface in common, not the best available one.
          bestn = E_IDX_NONE;
          break;
        }

        const E_Int* pgsj = ngi.PHs.get_facets_ptr(j);
        E_Int nb_pgsj = ngi.PHs.stride(j);

        cur_agg.clear();
        cur_ori.clear();

        K_MESH::Polyhedron<UNKNOWN>::merge_two_phs(crd, ngi.PGs, pgsi, nb_neighs, pgsj, nb_pgsj, orienti.get_facets_ptr(i), orienti.get_facets_ptr(j), cur_agg, cur_ori);

        bool ok = (K_MESH::Polyhedron<UNKNOWN>::is_centroid_star_shaped(dt, crd, ngi.PGs, cur_agg.get_facets_ptr(0), cur_agg.stride(0), cur_ori.get_facets_ptr(0)) == 0); 
        if (!ok)
          continue;
        
        reflex_edges.clear();
        convex_edges.clear();
        bool concave =false;
        E_Int err = K_MESH::Polyhedron<UNKNOWN>::is_concave(crd, ngi.PGs, cur_agg.get_facets_ptr(0), cur_agg.stride(0), false, cur_ori.get_facets_ptr(0), concave, 
                                                                      reflex_edges, convex_edges, 1./3 /*concave_threshold*/, 0.05/*convex_threshold*/);
        // is it better ?
        ngon_type::shared_faces_list(ngi, i, j, shared_pgs);
        size_t nbf = shared_pgs.size();
        // compute shared surface
        E_Float s=0.;
        for (size_t f = 0; f < nbf; ++f)
        {
          E_Int PGi = shared_pgs[f]-1;
          s += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }
        
        //fixme proto : compute the total surface of PHi
        E_Float stot=0.;
        for (size_t f = 0; f < nb_neighs; ++f)
        {
          E_Int PGi = *(pgsi+f)-1;
          stot += K_MESH::Polygon::surface<K_FLD::FloatArray, 3>(crd, ngi.PGs.get_facets_ptr(PGi), ngi.PGs.stride(PGi), 1);
        }
        
        E_Float face_ratio = E_Float(nbf) / E_Float (nb_neighs);
        E_Float surface_ratio = s /stot;
        
        E_Float worst_reflex_a=2.*K_CONST::E_PI;
                
        if (!reflex_edges.empty())
        {
          for (std::map<K_MESH::NO_Edge, E_Float>::const_iterator it = reflex_edges.begin(); it != reflex_edges.end(); ++it) 
            worst_reflex_a = std::min(worst_reflex_a, it->second);
        
          //std::cout << "worst angle : " << worst_reflex_a << std::endl;
        }
        
        E_Float q = face_ratio * surface_ratio * worst_reflex_a;
      
        //maximum surface is best. if equality, count the number of shared faces.
        //bool is_better = (s > smax) || ((::fabs(s-smax) < E_EPSILON) && (nbf > nbfmax));
        bool is_better = (q > qmax);
      
        if (is_better)
        {
          bestn = n;

          smax = s;
          nbfmax = nbf;
          best_agg = cur_agg;
          best_ori = cur_ori;
          
          qmax = q;
          
#ifdef DEBUG_AGGLOMERATOR 
         if (worst_reflex_angle > worst_reflex_a)
            std::cout << "better but generating more concavity !!" << std::endl;
#endif          
          worst_reflex_angle = worst_reflex_a;
          
#ifdef DEBUG_AGGLOMERATOR
          //if (ii==73)
          /*{
            K_FLD::IntArray cnto;
            ngon_type ng(ngi.PGs, best_agg);
            ng.export_to_array(cnto);
            std::ostringstream o;
            o << "agg_" << n << ".tp";
            MIO::write(o.str().c_str(), crd, cnto, "NGON"); 
          }*/
#endif
        } 
      }

      if (bestn == E_IDX_NONE) continue;
      
      E_Int j = *(neighs+bestn);
      frozen[i] = frozen[j] = true;
      
#ifdef DEBUG_AGGLOMERATOR
      //std::cout << "AGGLOMERATION : " << i << " with " << j << std::endl;
      //if (ii==73)
      {
        std::ostringstream o;
        o << "best_agg_" << ii << ".tp";
        K_FLD::IntArray cnto;
        ngon_type ng(ngi.PGs, best_agg);
        ng.export_to_array(cnto);
        MIO::write(o.str().c_str(), crd, cnto, "NGON");
      }
#endif

      all_agg.append(best_agg);
      all_ori.append(best_ori);
    }

    ngo.PGs = ngi.PGs;
    ngo.PHs = all_agg;
    
    //std::cout << ngo.PHs.size() << " small cells have been agglomerated." << std::endl;

    // now add untouched ones
    for (size_t i = 0; i < frozen.size(); ++i)
    {
      if (!frozen[i])
      {
        ngo.PHs.add(ngi.PHs.stride(i), ngi.PHs.get_facets_ptr(i));
        if (!ngi.PHs._type.empty()) ngo.PHs._type.push_back(ngi.PHs._type[i]);
      }
    }
    
    ngo.PGs.updateFacets();
    ngo.PHs.updateFacets();
    
    //std::cout << "effective nb of phs in agglomerated ngo : " << ngo.PHs.size() << std::endl;

    std::vector<E_Int> pgnids, phnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);

  }
  
  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::agglomerate_small_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, E_Float vmin, E_Float vratio, ngon_type& ngo, E_Int& nb_aggs)
  {
    ngon_unit neighborsi;
    ngi.build_ph_neighborhood(neighborsi);
    
    ngon_unit orienti;
    ngon_type::build_orientation_ngu<TriangulatorType>(crd, ngi, orienti);
    
    std::vector<E_Int> PHlist;
    ngon_type::detect_bad_volumes<TriangulatorType>(crd, ngi, neighborsi, vmin, vratio, PHlist);
    
    if (PHlist.empty()) 
    {
      ngo = ngi;
      return;
    }

    ngon_unit oriento;
    NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, orienti, PHlist, ngo, oriento, nb_aggs);
  }
  
//  ///
//  template<typename TriangulatorType>
//  void NUGA::Agglomerator::agglomerate_uncomputable_phs(const K_FLD::FloatArray& crd, ngon_type& ngi, ngon_type& ngo)
//  {
//    std::vector<E_Int> PHlist;
//    ngon_type::detect_uncomputable_phs<TriangulatorType>(crd, ngi, PHlist);
//    
//    if (PHlist.empty()) 
//    {
//      std::cout << "OK : There are no uncomputable phs" << std::endl;
//      ngo = ngi;
//      return;
//    }
//    
//    std::cout << PHlist.size() << " uncomputable phs were detected" << std::endl;
//    
//    ngon_unit neighborsi;
//    ngi.build_ph_neighborhood(neighborsi);
//
//    NUGA::Agglomerator::agglomerate_phs<TriangulatorType>(crd, ngi, neighborsi, PHlist, ngo);
//  }
 
  ///
  template<typename TriangulatorType>
  void NUGA::Agglomerator::collapse_uncomputable_pgs(K_FLD::FloatArray& crd, ngon_type& ngio)
  {
    std::vector<E_Int> PGlist;
    std::vector<ngon_type::ePathoPG> flagPG;
    ngon_type::detect_uncomputable_pgs<TriangulatorType>(crd, ngio.PGs, flagPG);

#ifdef DEBUG_AGGLOMERATOR
    std::vector<E_Int> badPGs;
    for (size_t i=0; i < flagPG.size(); ++i) if (flagPG[i] != ngon_type::NONE/*== ngon_type::SPIKE*/) badPGs.push_back(i);
    
    NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs(crd, ngio.PGs, badPGs);
    
    ngon_type tmp = ngio;
    tmp.keep_PHs_having_PGs(badPGs);
    K_FLD::IntArray ctmp;
    tmp.export_to_array(ctmp);
    MIO::write("bad_phs.plt", crd, ctmp, "NGON");
#endif
    
    ///  ONLY DEAL WITH T3 SPIKES AND MESHER FAILURES : same repsonse : 
    ///  COLLAPSING THE SMALLEST EDGE => NUGA::Agglomerator::collapse_pgs2
    E_Int nb_pathos(0);
    for (size_t i=0; i < flagPG.size(); ++i)
    {
      if (flagPG[i] != ngon_type::PATHO_PG_NONE) ++nb_pathos;
      
      if (flagPG[i] == ngon_type::SPIKE)
        PGlist.push_back(i);
      else if (flagPG[i] == ngon_type::DELAUNAY_FAILURE)
        PGlist.push_back(i);
    }
    
    if (nb_pathos == 0) 
      std::cout << "OK : COMPUTABLE MESH" << std::endl;
    else
      std::cout << PGlist.size() << "uncomputable pgs will be fixed over " <<  nb_pathos << " pathologies." << std::endl;
    
    //do the healing
    NUGA::Agglomerator::collapse_pgs2(crd, ngio, PGlist);
    
  }
  
  /// collapse to the barycenter
  E_Int NUGA::Agglomerator::collapse_pgs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids)
  {
    if (pgids.empty())
      return 0;
    
    E_Float G[3];
    std::vector<E_Int> nid(crd.cols());
    for (size_t i=0; i < nid.size(); ++i) nid[i]=i;
    E_Int ng_pgs = pgids.size();
  
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray > acrd_t;
    acrd_t acrd(crd);
  
    for (size_t i=0; i<pgids.size(); ++i)
    {
      E_Int PGi = pgids[i];
        
      K_MESH::Polygon::iso_barycenter<acrd_t, 3>(acrd, ng.PGs.get_facets_ptr(PGi), ng.PGs.stride(PGi), 1, G);
    
     E_Int Ni = ng.PGs.get_facet(PGi, 0)-1;
      for (size_t j=0; j < ng.PGs.stride(PGi); ++j)
        nid[ng.PGs.get_facet(PGi, j)-1]=Ni;
    
      crd(0,Ni)=G[0]; crd(1,Ni)=G[1]; crd(2,Ni)=G[2];
    }
  
    ng.PGs.change_indices(nid, true);
  
    ngon_type::clean_connectivity(ng, crd);
  
  return 0;
}
  /// collaspe the smallest edge keeping the oldest node.
  E_Int NUGA::Agglomerator::collapse_pgs2(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& pgids)
  {
    if (pgids.empty())
      return 0;
    
    E_Float G[3];
    std::vector<E_Int> nids(crd.cols());
    for (size_t i=0; i < nids.size(); ++i) nids[i]=i;
    E_Int ng_pgs = pgids.size();
  
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray > acrd_t;
    acrd_t acrd(crd);
  
    for (size_t i=0; i<pgids.size(); ++i)
    {
      E_Int PGi = pgids[i];
      const E_Int* nodes = ng.PGs.get_facets_ptr(PGi);
      E_Int nb_nodes = ng.PGs.stride(PGi);
      E_Float d2min = K_CONST::E_MAX_FLOAT;
      E_Int Nsource(E_IDX_NONE), Ntarget(E_IDX_NONE);
      
      // for the smallest edge of the current polygon, assign the non-moved and biggest id to the smallest one
      for (E_Int n=0; n<nb_nodes; ++n)
      {
        E_Int Ni = *(nodes + n) - 1;
        E_Int Nj = *(nodes + (n+1) % nb_nodes) - 1;
        E_Float d2 = K_FUNC::sqrDistance(crd.col(Ni), crd.col(Nj), 3);
        
        if (d2 < d2min)
        {
          Ntarget = (nids[Ni] == Ni && nids[Nj] == Nj) ? MIN(Ni, Nj) : (nids[Ni] == Ni) ? Ni : Nj;
          d2min=d2;
          Nsource = (Ntarget == Ni) ? Nj : Ni;//the other one
        }
      }
      
      if (nids[Nsource] != Nsource || nids[Ntarget] != Ntarget) continue; //both have moved so this polygon is discarded
      //degenerate an edge of this polygon the pg by setting all the nodes to 
      nids[Nsource]= Ntarget;
    }
  
    ng.PGs.change_indices(nids, true);
  
    ngon_type::clean_connectivity(ng, crd);
  
  return 0;
}

/// PRIVATE ///

  ///
  void NUGA::Agglomerator::__simplify_phs
  (const K_FLD::FloatArray& crd, const ngon_type& ngi, E_Int PHi, const ngon_unit& orienti, const ngon_unit& phneighborsi, E_Float angular_max, bool process_externals,
   ngon_unit& gagg_pgs, std::vector<E_Int>& nids, ngon_unit& wagg_pgs, std::map<E_Int, std::vector<E_Int> >& wneigh_to_faces)
  {
    // nids[i] == -1         =>  UNCHANGED
    // nids[i] == E_IDX_NONE =>  DELETED
    // nids[i] == j          =>  j-th PG in agg_pgs
    
    typedef std::map<E_Int, std::vector<E_Int> > map_t;
    map_t::const_iterator it;

    std::vector<E_Int> lnids, ori, pgsi;

    wneigh_to_faces.clear();
 
    E_Int nb_neighs = phneighborsi.stride(PHi);
    const E_Int* neighs = phneighborsi.get_facets_ptr(PHi);
    const E_Int* pgs = ngi.PHs.get_facets_ptr(PHi);
    const E_Int* orient = orienti.get_facets_ptr(PHi);

    for (E_Int n = 0; n < nb_neighs; ++n)
    {
      E_Int PHn = *(neighs + n);
      if (!process_externals && PHn == E_IDX_NONE)
        continue;
      //E_Int PGi = *(pgs + n) - 1;
      wneigh_to_faces[PHn].push_back(n);
    }

    if ((E_Int)wneigh_to_faces.size() == nb_neighs) //no possible agglo
      return;

    for (it = wneigh_to_faces.begin(); it != wneigh_to_faces.end(); ++it)
    {
      const E_Int& PHn = it->first;
      const std::vector<E_Int>& common_pg_pos = it->second;
      size_t sz = common_pg_pos.size();
      if (sz == 1)
        continue;

      ori.clear();  ori.resize(sz);
      pgsi.clear(); pgsi.resize(sz);

      size_t already_done = 0;
      for (size_t k = 0; k < sz; ++k)
      {
        E_Int ni = common_pg_pos[k];
        ori[k] = orient[ni];
        pgsi[k] = pgs[ni];
        if (nids[pgsi[k] - 1] != NONEVAL) ++already_done;
        
       // std::cout << "Face  : " << pgsi[k] << " / id : " << nids[pgsi[k] - 1] << " . done already ? : " << (nids[pgsi[k] - 1] != -1) << std::endl;
      }
      
#ifdef DEBUG_AGGLOMERATOR
      assert(already_done == 0 || already_done == sz);
#endif
      
      if (already_done)
        continue;

      K_MESH::Polygon::full_agglomerate(crd, ngi.PGs, &pgsi[0], E_Int(sz), angular_max, &ori[0], wagg_pgs, lnids /*,normals*/);
      
#ifdef DEBUG_AGGLOMERATOR
      assert (lnids.size() == sz);
      if (PHn == 31923)
      {
        std::vector<E_Int> tmppg = pgsi;
        K_CONNECT::IdTool::shift(tmppg, -1);
        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs(crd, ngi.PGs, tmppg);
        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH("PHi.tp", crd, ngi, PHi);
        NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH("PHn.tp", crd, ngi, PHn);
      }
#endif

      //local to global
      E_Int shft = gagg_pgs.size();
      assert (lnids.size() == sz);
      for (E_Int k=0; k < sz; ++k)
      {
        E_Int PGi = pgsi[k] - 1;
        E_Int nid = lnids[k];

        nids[PGi] = (nid == E_IDX_NONE) ? E_IDX_NONE : (nid >= sz) ? nid -sz +shft : UNCHANGED; // deleted ? agglomerated ? unchanged
        
#ifdef DEBUG_AGGLOMERATOR  
        std::string tmp = (nids[PGi] == E_IDX_NONE) ? "DELETED" : (nids[PGi] == UNCHANGED) ? "UNCHANGED" : "NEW";
        if (tmp != "NEW") std::cout << PGi << " status : " <<  tmp << std::endl;
        else std::cout << PGi << " status : " << nid -sz +shft << std::endl;
#endif
        
      }

      if (wagg_pgs.size())
      {
        gagg_pgs.append(wagg_pgs);
        gagg_pgs._type.resize(gagg_pgs.size(), INNER);
        gagg_pgs._ancEs.resize(2, gagg_pgs.size(), E_IDX_NONE);
      }
    }
    
#ifdef DEBUG_AGGLOMERATOR
    assert (nids.size() == ngi.PGs.size());
#endif
  }
}

#endif