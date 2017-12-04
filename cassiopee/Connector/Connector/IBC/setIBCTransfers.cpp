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
# include "connector.h"
# include <math.h>
using namespace std;
using namespace K_FLD;

# include "IBC/commonLaws.h"

//=============================================================================
//Retourne -2: incoherence entre meshtype et le type d interpolation
//         -1: type invalide
//          1: ok
// Entree/Sortie: des variables conservatives ( + ronutildeSA ) 
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar1(
  E_Int bctype,
  E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
  E_Float* xPC, E_Float* yPC, E_Float* zPC,
  E_Float* xPW, E_Float* yPW, E_Float* zPW,
  E_Float* xPI, E_Float* yPI, E_Float* zPI, E_Float* densPtr, E_Float* pressPtr, E_Float* utauPtr, E_Float* yplusPtr,
  E_Float* tmp, E_Int& size,
  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts,
  vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields)
{
  E_Float gam1 = gamma-1.;
  /* lois de paroi */
  E_Float roext, uext, pext, text, muext, yext, yplus, yibc;
  E_Float uscaln, un, vn, wn, ut, vt, wt, utau, utauv, utau0, umod;
  E_Float expy, denoml10, ax, l1, l2, l3;
  E_Float aa, bb, dd, fp, tp, f1v, f1p;
  E_Float ucible0, ucible, vcible, wcible, nutcible, nutilde, signibc;
  E_Int npass;
  //Lois de paroi : criteres d arret pour estimer le frottement par Newton
  E_Float newtoneps = 1.e-7; // critere d arret pour u+
  E_Float newtonepsprime = 1.e-12;// critere d arret pour la derivee  
  E_Float cvgaminv = 1./(cv*gam1);
  E_Float coefSuth = muS * (1.+Cs/Ts);
  E_Float Tsinv = 1./Ts;
  E_Float kappa = 0.4; // Constante de Von Karman
  E_Float kappainv = 1./kappa;
  E_Float cc = 5.2;//pour la loi log
  /* fin parametres loi de parois */

  E_Int nvars = vectOfDnrFields.size();
  //E_Int nbRcvPts = rcvPtsI.getSize();

  E_Float a0,a1,a2,b0,b1,b2,n0,n1,n2;
  E_Float normb, ro, u,v,w, p, roE;
  E_Float vnc, alpha, beta, alphasbeta;
  E_Float* roOut  = vectOfRcvFields[0];// ro
  E_Float* rouOut = vectOfRcvFields[1];// rou
  E_Float* rovOut = vectOfRcvFields[2];// rov
  E_Float* rowOut = vectOfRcvFields[3];// row
  E_Float* roEOut = vectOfRcvFields[4];// roE
  E_Float* varSAOut = NULL;
  if (nvars == 6) varSAOut = vectOfRcvFields[5]; // ronutildeSA

  //E_Int* rcvPts = rcvPtsI.begin();

  if (bctype == 0) // wallslip
  { 
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
           E_Int indR = rcvPts[noind+ideb];

            roext  = roOut[indR];     //ro
            roE    = roEOut[indR];    //roE
            E_Float roinv = 1./roext;
            // vitesse
            u = rouOut[indR]*roinv; v = rovOut[indR]*roinv; w = rowOut[indR]*roinv;

# include "IBC/commonBCType0.h"
      
           rouOut[indR] = ucible*roext; rovOut[indR] = vcible*roext; rowOut[indR] = wcible*roext;      
           // pression
           pext            = gam1*(roE-0.5*roext*(u*u+v*v+w*w));//PI
           roEOut[indR]    = pext/gam1+0.5*roext*(ucible*ucible+vcible*vcible+wcible*wcible);
           pressPtr[noind +ideb] = pext; 
           densPtr[ noind +ideb] = roext;

           if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;
        }
  }
  else if (bctype == 1) // adherence
  {
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
           E_Int indR = rcvPts[noind+ideb];

            roext  = roOut[indR];      //ro
            roE    = roEOut[indR];     //roE
            E_Float roinv = 1./roext;
            // vitesse
            u = rouOut[indR]*roinv; v = rovOut[indR]*roinv; w = rowOut[indR]*roinv;

# include "IBC/commonBCType1.h"
      
           rouOut[indR] = ucible*roext; rovOut[indR] = vcible*roext; rowOut[indR] = wcible*roext;      
           // pression
           pext         = gam1*(roE-0.5*roext*(u*u+v*v+w*w));//PI
           roEOut[indR] = pext/gam1+0.5*roext*(ucible*ucible+vcible*vcible+wcible*wcible);

           pressPtr[noind +ideb] = pext; 
           densPtr[ noind +ideb] = roext;
           //utau pas calcule ni y+
           //
           if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;
        }
  }
  else if (bctype == 2) // loi de paroi en log
  {
#     include "IBC/pointer.h" 

      E_Int err  = 0;
      E_Int skip = 0; 
      //initilisation parametre geometrique et utau
#ifdef _OPENMP4
      #pragma omp simd
#endif 
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
       {
          //E_Int indR = rcvPts[noind];
          E_Int indR = rcvPts[noind+ideb];
 
          roext = roOut[indR]; // densite du point interpole
          roE   = roEOut[indR];    //roE

          E_Float roinv = 1./roext;
          // vitesse du pt ext
          u = rouOut[indR]*roinv;
          v = rovOut[indR]*roinv;
          w = rowOut[indR]*roinv;

          pext = gam1*(roE-0.5*roext*(u*u+v*v+w*w));//PI
          text = pext / roext * cvgaminv;

#         include "IBC/commonLogLaw_init.h" 
          // out= utau  et err 
       }

      // Newton pour utau
#     include "IBC/commonLogLaw_Newton.h" 

      //initialisation Newton SA  + vitesse cible
#     include "IBC/commonLogLaw_cible.h" 
      if (nvars == 6)
      {
          // Newton pour mut
#         include "IBC/nutildeSA_Newton.h" 

          // mise a jour des variable
#ifdef _OPENMP4
         #pragma omp simd
#endif 
          for (E_Int noind = 0; noind < ifin-ideb; noind++)
          {
           E_Int indR = rcvPts[noind+ideb];

           roE   = roEOut[indR];//roE

           E_Float roinv = 1./ro_vec[noind];
           // vitesse du pt ext
           u = rouOut[indR]*roinv;
           v = rovOut[indR]*roinv;
           w = rowOut[indR]*roinv;

           rouOut[indR] = ucible_vec[noind]*ro_vec[noind]; 
           rovOut[indR] = vcible_vec[noind]*ro_vec[noind];
           rowOut[indR] = wcible_vec[noind]*ro_vec[noind];      

           roEOut[indR] = roE + 0.5*ro_vec[noind]*( ucible_vec[noind]*ucible_vec[noind] - u*u
                                                   +vcible_vec[noind]*vcible_vec[noind] - v*v
                                                   +wcible_vec[noind]*wcible_vec[noind] - w*w );

           varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind]*ro_vec[noind];                                    //nutilde*signibc    
          }
      }
      else //5eq 
      {
        // mise a jour des variable
#ifdef _OPENMP4
         #pragma omp simd
#endif 
          for (E_Int noind = 0; noind < ifin-ideb; noind++)
          {
           E_Int indR = rcvPts[noind+ideb];

           E_Float roinv = 1./ro_vec[noind];
           // vitesse du pt ext
           u = rouOut[indR]*roinv;
           v = rovOut[indR]*roinv;
           w = rowOut[indR]*roinv;

           rouOut[indR] = ucible_vec[noind]*ro_vec[noind]; rovOut[indR] = vcible_vec[noind]*ro_vec[noind]; rowOut[indR] = wcible_vec[noind]*ro_vec[noind];      
           roE          = roEOut[indR];//roE
           roEOut[indR] = roE + 0.5*ro_vec[noind]*( ucible_vec[noind]*ucible_vec[noind] - u*u
                                                   +vcible_vec[noind]*vcible_vec[noind] - v*v
                                                   +wcible_vec[noind]*wcible_vec[noind] - w*w );
          }
      }
  }
  else if (bctype == 3)// loi de paroi  Musker
  {
#   include "IBC/pointer.h" 

    E_Int err  = 0;
    E_Int skip = 0; 
    //initilisation parametre geometrique et utau
#ifdef _OPENMP4
    #pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
        //E_Int indR = rcvPts[noind];
        E_Int indR = rcvPts[noind+ideb];
 
        roext = roOut[indR]; // densite du point interpole
        roE   = roEOut[indR];    //roE

        E_Float roinv = 1./roext;
        // vitesse du pt ext
        u = rouOut[indR]*roinv;
        v = rovOut[indR]*roinv;
        w = rowOut[indR]*roinv;

        
        pext = gam1*(roE-0.5*roext*(u*u+v*v+w*w));//PI
        text = pext / roext * cvgaminv;

#       include "IBC/commonMuskerLaw_init.h" 
        // out= utau  et err 
    }

    // Newton pour utau
#   include "IBC/commonMuskerLaw_Newton.h" 

    //initialisation Newton SA  + vitesse cible
#   include "IBC/commonMuskerLaw_cible.h" 
    if (nvars == 6)
    { 
      err  = 0; skip = 0;
// Newton pour mut
#     include "IBC/nutildeSA_Newton.h" 

      // mise a jour des variable
#ifdef _OPENMP4
      #pragma omp simd
#endif 
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
         E_Int indR = rcvPts[noind+ideb];

         roE   = roEOut[indR];//roE

         E_Float roinv = 1./ro_vec[noind];
         // vitesse du pt ext
         u = rouOut[indR]*roinv;
         v = rovOut[indR]*roinv;
         w = rowOut[indR]*roinv;

         rouOut[indR] = ucible_vec[noind]*ro_vec[noind]; 
         rovOut[indR] = vcible_vec[noind]*ro_vec[noind];
         rowOut[indR] = wcible_vec[noind]*ro_vec[noind];      

         roEOut[indR] = roE + 0.5*ro_vec[noind]*( ucible_vec[noind]*ucible_vec[noind] - u*u
                                                 +vcible_vec[noind]*vcible_vec[noind] - v*v
                                                 +wcible_vec[noind]*wcible_vec[noind] - w*w );

         varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind]*ro_vec[noind];                                    //nutilde*signibc    
      }
    }
    else //5eq 
    { // mise a jour des variable
#ifdef _OPENMP4
      #pragma omp simd
#endif 
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
         E_Int indR = rcvPts[noind+ideb];

         roE   = roEOut[indR];//roE

         E_Float roinv = 1./ro_vec[noind];
         // vitesse du pt ext
         u = rouOut[indR]*roinv; 
         v = rovOut[indR]*roinv;
         w = rowOut[indR]*roinv;

         rouOut[indR] = ucible_vec[noind]*ro_vec[noind]; rovOut[indR] = vcible_vec[noind]*ro_vec[noind]; rowOut[indR] = wcible_vec[noind]*ro_vec[noind];      

         roEOut[indR] = roE + 0.5*ro_vec[noind]*( ucible_vec[noind]*ucible_vec[noind] - u*u
                                                 +vcible_vec[noind]*vcible_vec[noind] - v*v
                                                 +wcible_vec[noind]*wcible_vec[noind] - w*w );
      }
    }
  }
  return 1;
}
//=============================================================================
//Retourne -2: incoherence entre meshtype et le type d'interpolation
//         -1: type invalide
//          1: ok
// Entree/Sortie:  (ro,u,v,w,t) ( + nutildeSA ) 
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar2(
  E_Int bctype,
  E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
  E_Float* xPC, E_Float* yPC, E_Float* zPC,
  E_Float* xPW, E_Float* yPW, E_Float* zPW,
  E_Float* xPI, E_Float* yPI, E_Float* zPI, E_Float* densPtr, E_Float* pressPtr, E_Float* utauPtr, E_Float* yplusPtr,
  E_Float* tmp, E_Int& size,
  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts,
  vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields)
{
  /* lois de paroi */
  E_Float roext, uext, pext, text, muext, yext, yplus, yibc;
  E_Float uscaln, un, vn, wn, ut, vt, wt, utau, utauv, utau0, umod;
  E_Float aa, bb, dd, fp, tp, f1v, f1p;
  E_Float expy, denoml10,ax,l1,l2, l3;
  E_Float ucible0, ucible, vcible, wcible, nutcible, nutilde, signibc;
  E_Int npass;
  // Lois de paroi: criteres d'arret pour estimer le frottement par Newton
  E_Float newtoneps = 1.e-7; // critere d'arret pour u+
  E_Float newtonepsprime = 1.e-12;// critere d'arret pour la derivee  
  E_Float cvgam = cv*(gamma-1.);
  E_Float coefSuth = muS * (1.+Cs/Ts);
  E_Float Tsinv = 1./Ts;
  E_Float kappa = 0.4; // Constante de Von Karman
  E_Float kappainv = 1./kappa;
  E_Float cc = 5.2; //pour la loi log
  /* fin parametres loi de parois*/


  E_Int nvars = vectOfDnrFields.size();

  E_Float a0,a1,a2,b0,b1,b2,n0,n1,n2;
  E_Float normb, u, v, w;
  E_Float vnc, alpha, beta, alphasbeta;
  E_Float* roOut = vectOfRcvFields[0];// ro
  E_Float* uOut  = vectOfRcvFields[1];// u
  E_Float* vOut  = vectOfRcvFields[2];// v
  E_Float* wOut  = vectOfRcvFields[3];// w
  E_Float* tOut  = vectOfRcvFields[4];// temperature
  E_Float* varSAOut = NULL;
  if (nvars == 6) varSAOut = vectOfRcvFields[5]; // nutildeSA

  if (bctype == 0) // wallslip
  { 
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         // vitesse
         u = uOut[indR];
         v = vOut[indR];
         w = wOut[indR];
# include "IBC/commonBCType0.h"
         uOut[indR] = ucible; 
         vOut[indR] = vcible; 
         wOut[indR] = wcible;      
         if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

         pressPtr[noind +ideb] = roOut[indR]* tOut[indR]*cvgam; 
         densPtr[ noind +ideb] = roOut[indR];
        }
  }
  else if (bctype == 1) // adherence (lineaire)
  {
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         // vitesse
         u = uOut[indR]; 
         v = vOut[indR];
         w = wOut[indR];
# include "IBC/commonBCType1.h"
         uOut[indR] = ucible; 
         vOut[indR] = vcible; 
         wOut[indR] = wcible;      
         if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

         pressPtr[noind +ideb] = roOut[indR]* tOut[indR]*cvgam; 
         densPtr[ noind +ideb] = roOut[indR];
        }
  }
  else if (bctype == 2) // loi de paroi log
  {
#   include "IBC/pointer.h" 

    E_Int err  = 0;
    E_Int skip = 0; 
    //initilisation parametre geometrique et utau
#ifdef _OPENMP4
    #pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
     {
        //E_Int indR = rcvPts[noind];
        E_Int indR = rcvPts[noind+ideb];
 
        roext = roOut[indR]; // densite du point interpole
        text  = tOut[indR];  // pression du point interpole
        pext  = text*roext*cvgam;

        // vitesse du pt ext
        u = uOut[indR];
        v = vOut[indR];
        w = wOut[indR];

#       include "IBC/commonLogLaw_init.h" 
        // out= utau  et err 
     }

     // Newton pour utau
#    include "IBC/commonLogLaw_Newton.h" 

     //initialisation Newton SA  + vitesse cible
#    include "IBC/commonLogLaw_cible.h" 
    if (nvars == 6)
      {
        // Newton pour mut
#       include "IBC/nutildeSA_Newton.h" 

        // mise a jour des variable
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         uOut[indR]     = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR]  = wcible_vec[noind];
         varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];                                         //nutilde*signibc    
        }
      }
    else //5eq 
      {
        // mise a jour des variable
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         uOut[indR] = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];
        }
      }
  }
  else if (bctype == 3) // loi de paroi Musker
  {
#   include "IBC/pointer.h" 

    E_Int err  = 0;
    E_Int skip = 0; 
    //initialisation parametre geometrique et utau
#ifdef _OPENMP4
    #pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
     {
        //E_Int indR = rcvPts[noind];
        E_Int indR = rcvPts[noind+ideb];
 
        roext = roOut[indR]; // densite du point interpole
        text  = tOut[indR];  // pression du point interpole
        pext  = text*roext*cvgam;

        // vitesse du pt ext
        u = uOut[indR];
        v = vOut[indR]; 
        w = wOut[indR];
        //printf("IN WALL LAW: %f %f %f %f %f \n",roext, text, u,v,w);
#       include "IBC/commonMuskerLaw_init.h" 
        // out= utau  et err 
      
     }  

     // Newton pour utau
#    include "IBC/commonMuskerLaw_Newton.h" 

     //initialisation Newton SA  + vitesse cible
#    include "IBC/commonMuskerLaw_cible.h" 

    if (nvars == 6)
      {
        // Newton pour mut
#       include "IBC/nutildeSA_Newton.h" 

        // mise a jour des variable
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
          E_Int indR = rcvPts[noind+ideb];

          uOut[indR]     = ucible_vec[noind];
          vOut[indR]     = vcible_vec[noind];
          wOut[indR]     = wcible_vec[noind];
          varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];  //nutilde*signibc
          // printf("OUT WALL LAW: %f %f %f %f\n",uOut[indR],vOut[indR],wOut[indR],varSAOut[indR]);
        }
      }
    else //5eq 
      {
        // mise a jour des variable
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         uOut[indR]     = ucible_vec[noind];
         vOut[indR]     = vcible_vec[noind];
         wOut[indR]     = wcible_vec[noind];
        }
      }

  }//bctype

  return 1;
}
//=============================================================================
//Retourne -2 : incoherence entre meshtype et le type d interpolation
//         -1 : type invalide
//          1 : ok
// Entree/Sortie :  (ro,u,v,w,p) ( + ronutildeSA ) 
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar3(
  E_Int bctype,
  E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin,  E_Int& ithread,
  E_Float* xPC, E_Float* yPC, E_Float* zPC,
  E_Float* xPW, E_Float* yPW, E_Float* zPW,
  E_Float* xPI, E_Float* yPI, E_Float* zPI, E_Float* densPtr, E_Float* pressPtr, E_Float* utauPtr, E_Float* yplusPtr,
  E_Float* tmp, E_Int& size,
  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts,
  vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields)
{
  /* lois de paroi */
  E_Float roext, uext, pext, text, muext, yext, yplus, yibc;
  E_Float uscaln, un, vn, wn, ut, vt, wt, utau, utauv, utau0, umod;
  E_Float aa, bb, dd, fp, tp, f1v, f1p;
  E_Float expy, denoml10,ax,l1,l2, l3;
  E_Float ucible0, ucible, vcible, wcible, nutcible, nutilde, signibc;
  E_Int npass;
  //Lois de paroi : criteres d arret pour estimer le frottement par Newton
  E_Float newtoneps = 1.e-7; // critere d arret pour u+
  E_Float newtonepsprime = 1.e-12;// critere d arret pour la derivee  
  E_Float cvgaminv = 1./(cv*(gamma-1.));
  E_Float coefSuth = muS * (1.+Cs/Ts);
  E_Float Tsinv = 1./Ts;
  E_Float kappa = 0.4; // Constante de Von Karman
  E_Float kappainv = 1./kappa;
  E_Float cc = 5.2;//pour la loi log
  /* fin parametres loi de parois*/

  E_Int nvars = vectOfDnrFields.size();

  //E_Float a[3], b[3], n[3];
  E_Float a0,a1,a2,b0,b1,b2,n0,n1,n2;
  E_Float normb, u,v,w;
  E_Float vnc, alpha, beta, alphasbeta;
  E_Float* roOut = vectOfRcvFields[0];// density
  E_Float* uOut = vectOfRcvFields[1];// ux
  E_Float* vOut = vectOfRcvFields[2];// uy
  E_Float* wOut = vectOfRcvFields[3];// uz
  E_Float* pOut = vectOfRcvFields[4];// pressure
  E_Float* varSAOut = NULL;
  if ( nvars == 6 ) varSAOut = vectOfRcvFields[5];//ronutildeSA

  if (bctype == 0) // wallslip
  {
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         // vitesse
         u = uOut[indR];
         v = vOut[indR];
         w = wOut[indR];
# include "IBC/commonBCType0.h"
         uOut[indR] = ucible; 
         vOut[indR] = vcible; 
         wOut[indR] = wcible;      
         if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

         pressPtr[noind +ideb] = pOut[indR]; 
         densPtr[ noind +ideb] = roOut[indR];
        }
  }
  else if (bctype == 1) // adherence
  {
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         // vitesse
         u = uOut[indR];
         v = vOut[indR];
         w = wOut[indR];
# include "IBC/commonBCType1.h"
         uOut[indR] = ucible; 
         vOut[indR] = vcible; 
         wOut[indR] = wcible;      
         if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

         pressPtr[noind +ideb] = pOut[indR]; 
         densPtr[ noind +ideb] = roOut[indR];

        }
  }
  else if (bctype == 2) // loi de paroi en log
  {
#   include "IBC/pointer.h" 

    E_Int err  = 0;
    E_Int skip = 0; 
    //initilisation parametre geometrique et utau
#ifdef _OPENMP4
    #pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
        E_Int indR = rcvPts[noind+ideb];
 
        roext= roOut[indR]; // densite du point interpole
        pext = pOut[indR];  // pression du point interpole
        text = pext / roext * cvgaminv;

        // vitesse du pt ext
        u = uOut[indR];
        v = vOut[indR];
        w = wOut[indR];

#       include "IBC/commonLogLaw_init.h" 
        // out= utau  et err 
    }

    // Newton pour utau
#   include "IBC/commonLogLaw_Newton.h" 

    //initialisation Newton SA  + vitesse cible
#   include "IBC/commonLogLaw_cible.h" 
    if (nvars == 6)
    {
        // Newton pour mut
#       include "IBC/nutildeSA_Newton.h" 

        // mise a jour des variable
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         uOut[indR]     = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR]  = wcible_vec[noind];
         varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];                                         //nutilde*signibc    
        }
    }
    else //5eq 
    {
        // mise a jour des variable
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         uOut[indR] = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];
        }
    }
  }
  else if (bctype == 3)// loi de paroi de Musker
  {
#   include "IBC/pointer.h" 

    E_Int err  = 0;
    E_Int skip = 0; 
    //initilisation parametre geometrique et utau
#ifdef _OPENMP4
    #pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
        E_Int indR = rcvPts[noind+ideb];
 
        roext= roOut[indR]; // densite du point interpole
        pext = pOut[indR];  // pression du point interpole
        text = pext / roext * cvgaminv;

        // vitesse du pt ext
        u = uOut[indR];
        v = vOut[indR];
        w = wOut[indR];

#       include "IBC/commonMuskerLaw_init.h" 
        // out= utau  et err 
    }

    // Newton pour utau
#   include "IBC/commonMuskerLaw_Newton.h" 

    //initialisation Newton SA  + vitesse cible
#   include "IBC/commonMuskerLaw_cible.h" 
    if (nvars == 6)
    {
        // Newton pour mut
#       include "IBC/nutildeSA_Newton.h" 
       // mise a jour des variable
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         uOut[indR]     = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR]  = wcible_vec[noind];
         varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];                                         //nutilde*signibc    
        }
    }
    else //5eq 
    {
        // mise a jour des variable
#ifdef _OPENMP4
       #pragma omp simd
#endif 
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         uOut[indR] = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];
        }
    }
  }
  return 1;
}
//=============================================================================
/* Effectue les transfers IBC */
//=============================================================================
//Stephanie: fonction absente du pytree?
PyObject* K_CONNECTOR::setIBCTransfers(PyObject* self, PyObject* args)
{
  PyObject *arrayR, *arrayD;
  PyObject *pyVariables; 
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayXPC, *pyArrayXPI, *pyArrayXPW;
  PyObject *pyArrayYPC, *pyArrayYPI, *pyArrayYPW;
  PyObject *pyArrayZPC, *pyArrayZPI, *pyArrayZPW;
  PyObject *pyArrayPressure, *pyArrayUtau, *pyArrayYplus, *pyArrayDens;
  E_Int bctype;
  E_Int vartype;
  E_Float gamma, cv, muS, Cs, Ts;

  if (!PYPARSETUPLE(args,
                    "OOOOOOOOOOOOOOOOOOOOllddddd", "OOOOOOOOOOOOOOOOOOOOiiddddd",
                    "OOOOOOOOOOOOOOOOOOOOllfffff", "OOOOOOOOOOOOOOOOOOOOiifffff",
                    &arrayR, &arrayD,  &pyVariables,
                    &pyIndRcv, &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, 
                    &pyArrayXPC, &pyArrayYPC, &pyArrayZPC,
                    &pyArrayXPW, &pyArrayYPW, &pyArrayZPW,
                    &pyArrayXPI, &pyArrayYPI, &pyArrayZPI,
                    &pyArrayPressure, &pyArrayUtau, &pyArrayYplus, &pyArrayDens, 
                    &bctype, &vartype, &gamma, &cv, &muS, &Cs, &Ts))
  {
      return NULL;
  }
    
  E_Int bcType = E_Int(bctype); // 0 : wallslip; 1: noslip; 2: log law of wall; 3: Musker law of wall
  /* varType : 
     1  : conservatives, 
     11 : conservatives + ronutildeSA 
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA 
     3  : (ro,u,v,w,p)     
     31 : (ro,u,v,w,p) + ronutildeSA */
  E_Int varType = E_Int(vartype); 

  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringR; char* eltTypeR;
  E_Int resr = K_ARRAY::getFromArray(arrayR, varStringR, fr, 
                                     imr, jmr, kmr, cnr, eltTypeR, true); 
  if (resr != 2 && resr != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "setIBCTransfers: 1st arg is not a valid array.");
    return NULL; 
  }
  /*---------------------------------------------*/
  /* Extraction des infos sur le domaine donneur */
  /*---------------------------------------------*/
  E_Int imd, jmd, kmd, imdjmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringD; char* eltTypeD;
  E_Int resd = K_ARRAY::getFromArray(arrayD, varStringD, fd, 
                                     imd, jmd, kmd, cnd, eltTypeD, true); 
  E_Int* ptrcnd = NULL;
  E_Int cndSize = 0; E_Int cnNfldD = 0;
  if (resd != 2 && resd != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "setIBCTransfers: 2nd arg is not a valid array.");
    RELEASESHAREDB(resr, arrayR, fr, cnr); 
    return NULL; 
  }
  E_Int meshtype = resd;// 1 : structure, 2 non structure
  if (resd == 2)
  {
    if (strcmp(eltTypeD,"TETRA") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "setIBCTransfers: donor array must be a TETRA if unstructured.");
      RELEASESHAREDB(resr, arrayR, fr, cnr); 
      RELEASESHAREDB(resd, arrayD, fd, cnd); 
      return NULL; 
    }
    ptrcnd  = cnd->begin();
    cndSize = cnd->getSize();
    cnNfldD = cnd->getNfld();
  }

  E_Int nvars;
  if ( varType == 1 || varType == 2 || varType == 3 ) 
    nvars = 5;
  else if ( varType == 11 || varType == 21 || varType == 31 )
    nvars = 6;
  else 
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr); 
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfers: varType value is not valid.");
    return NULL;
  }

# include "extract_interpD.h"
  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  E_Int res_rcv = K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI, true);
  nbRcvPts      = rcvPtsI->getSize();
  E_Int* rcvPts = rcvPtsI->begin();

  if (res_donor*res_type*res_coef*res_rcv ==0) 
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr); 
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    if (res_donor != 0) { RELEASESHAREDN(pyIndDonor  , donorPtsI  );}
    if (res_type  != 0) { RELEASESHAREDN(pyArrayTypes, typesI     );}
    if (res_coef  != 0) { RELEASESHAREDN(pyArrayCoefs, donorCoefsF);}
    if (res_rcv   != 0) { RELEASESHAREDN(pyIndRcv    , rcvPtsI    );}
    PyErr_SetString(PyExc_TypeError,"setInterpTransfersD: 4th to 6th arg must be a numpy of integers. 7th arg a numpy floats ");
    return NULL;
  }

# include "IBC/extract_IBC.h"

  if (okc1*okc2*okc3 == 0 )
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    RELEASESHAREDN(pyIndDonor  , donorPtsI  );
    RELEASESHAREDN(pyArrayTypes, typesI     );
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    if (okc1 != 0) { RELEASESHAREDN(pyArrayXPC  , coordxPC  );}
    if (okc2 != 0) { RELEASESHAREDN(pyArrayYPC  , coordyPC  );}
    if (okc3 != 0) { RELEASESHAREDN(pyArrayZPC  , coordzPC  );}
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfersD: coordinates of corrected points are invalid.");
    return NULL;
  }
  if (okw1*okw2*okw3 == 0 )
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    RELEASESHAREDN(pyIndDonor  , donorPtsI  );
    RELEASESHAREDN(pyArrayTypes, typesI     );
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    RELEASESHAREDN(pyArrayXPC  , coordxPC  );
    RELEASESHAREDN(pyArrayYPC  , coordyPC  );
    RELEASESHAREDN(pyArrayZPC  , coordzPC  );
    if (okw1 != 0) { RELEASESHAREDN(pyArrayXPW  , coordxPW  );}
    if (okw2 != 0) { RELEASESHAREDN(pyArrayYPW  , coordyPW  );}
    if (okw3 != 0) { RELEASESHAREDN(pyArrayZPW  , coordzPW  );}
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfersD: coordinates of wall points are invalid.");
    return NULL;
  }
  if (oki1*oki2*oki3 == 0 )
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    RELEASESHAREDN(pyIndDonor  , donorPtsI  );
    RELEASESHAREDN(pyArrayTypes, typesI     );
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    RELEASESHAREDN(pyArrayXPC  , coordxPC  );
    RELEASESHAREDN(pyArrayYPC  , coordyPC  );
    RELEASESHAREDN(pyArrayZPC  , coordzPC  );
    RELEASESHAREDN(pyArrayXPW  , coordxPW  );
    RELEASESHAREDN(pyArrayYPW  , coordyPW  );
    RELEASESHAREDN(pyArrayZPW  , coordzPW  );
    if (oki1 != 0) { RELEASESHAREDN(pyArrayXPI  , coordxPI  );}
    if (oki2 != 0) { RELEASESHAREDN(pyArrayYPI  , coordyPI  );}
    if (oki3 != 0) { RELEASESHAREDN(pyArrayZPI  , coordzPI  );}
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfersD: coordinates of interpolated points are invalid.");
    return NULL;
  }
  if (okD*okP*okU*okY == 0 )
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    RELEASESHAREDN(pyIndDonor  , donorPtsI  );
    RELEASESHAREDN(pyArrayTypes, typesI     );
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    RELEASESHAREDN(pyArrayXPC  , coordxPC  );
    RELEASESHAREDN(pyArrayYPC  , coordyPC  );
    RELEASESHAREDN(pyArrayZPC  , coordzPC  );
    RELEASESHAREDN(pyArrayXPW  , coordxPW  );
    RELEASESHAREDN(pyArrayYPW  , coordyPW  );
    RELEASESHAREDN(pyArrayZPW  , coordzPW  );
    RELEASESHAREDN(pyArrayXPI  , coordxPI  );
    RELEASESHAREDN(pyArrayYPI  , coordyPI  );
    RELEASESHAREDN(pyArrayZPI  , coordzPI  );
    if (okD != 0) { RELEASESHAREDN(pyArrayDens    , densF   );}
    if (okP != 0) { RELEASESHAREDN(pyArrayPressure, pressF  );}
    if (okU != 0) { RELEASESHAREDN(pyArrayUtau    , utauF   );}
    if (okY != 0) { RELEASESHAREDN(pyArrayYplus   , yplusF  );}
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfersD: Post array are invalid.");
    return NULL;
  }

  
  // Transferts 
  // Types valides: 2, 3, 4, 5 
  PyObject* tpl;
  if (resr == 1) 
  {
    tpl = K_ARRAY::buildArray(fr->getNfld(), varStringR, imr, jmr, kmr);
  }
  else // unstructured 
  {
    E_Int crsize = cnr->getSize()*cnr->getNfld(); 
    tpl = K_ARRAY::buildArray(fr->getNfld(), varStringR,
                              fr->getSize(), cnr->getSize(),
                              -1, eltTypeR, false, crsize);
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cnr->begin(), cnr->getSize()*cnr->getNfld());
  }
  E_Float* ptrFieldOut = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fieldROut(fr->getSize(), fr->getNfld(), ptrFieldOut, true);
  fieldROut = *fr; 

  vector<E_Float*> vectOfDnrFields;  
  vector<E_Float*> vectOfRcvFields;
  E_Int posvr, posvd;
  // Extrait les positions des variables a transferer
  E_Int nfoundvar = 0;
  if (PyList_Check(pyVariables) != 0)
  {
    int nvariables = PyList_Size(pyVariables);
    if (nvariables > 0)
    {
      for (int i = 0; i < nvariables; i++)
      {
        PyObject* tpl0 = PyList_GetItem(pyVariables, i);
        if (PyString_Check(tpl0) == 0)
          PyErr_Warn(PyExc_Warning, "setIBCTransfers: variable must be a string. Skipped.");
        else 
        {
          char* varname = PyString_AsString(tpl0);        
          posvd = K_ARRAY::isNamePresent(varname, varStringD);      
          posvr = K_ARRAY::isNamePresent(varname, varStringR);      
          if (posvd != -1 && posvr != -1) 
          {
            vectOfDnrFields.push_back(fd->begin(posvd+1));
            vectOfRcvFields.push_back(fieldROut.begin(posvr+1));
            nfoundvar += 1;
          }
        }
      }
    }// nvariables > 0
  }
  if ( nfoundvar != nvars )
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr); 
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    BLOCKRELEASEMEM;
    BLOCKRELEASEMEM2;
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfers: number of variables is inconsistent with varType.");
    return NULL;
  }

////
////
//  Interpolation parallele
////  
////  
# include "commonInterpTransfers_indirect.h"


    E_Int threadmax_sdm  = __NUMTHREADS__;

    E_Int size = (nbRcvPts/threadmax_sdm)+1; // on prend du gras pour gerer le residus
    E_Int    r =  size % 8;
          size = size + r;                  // on rajoute du bas pour alignememnt 64bits

    FldArrayF  tmp(size*13*threadmax_sdm);
    E_Float* ipt_tmp=  tmp.begin();

#pragma omp parallel default(shared)
  {
   //indice loop pour paralelisation omp
   E_Int ideb, ifin;
#ifdef _OPENMP
   E_Int  ithread           = omp_get_thread_num()+1;
   E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
   E_Int ithread = 1;
   E_Int Nbre_thread_actif = 1;
#endif
   // Calcul du nombre de champs a traiter par chaque thread
   E_Int chunk = nbRcvPts/Nbre_thread_actif;
   E_Int r = nbRcvPts - chunk*Nbre_thread_actif;
   // pts traitees par thread
   if (ithread <= r) 
        { ideb = (ithread-1)*(chunk+1); ifin = ideb + (chunk+1); }  
   else { ideb = (chunk+1)*r+(ithread-r-1)*chunk; ifin = ideb + chunk; }

   if (varType == 1 || varType == 11) 
     setIBCTransfersCommonVar1(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread, 
			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, density, pressure, utau, yplus,
                              ipt_tmp, size,
                              gamma, cv, muS, Cs, Ts,
	      	              vectOfDnrFields, vectOfRcvFields);
  
   else if (varType == 2 || varType == 21) 
     setIBCTransfersCommonVar2(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, density, pressure, utau, yplus,
                              ipt_tmp, size,
                              gamma, cv, muS, Cs, Ts,
		              vectOfDnrFields, vectOfRcvFields);

   else if (varType == 3 || varType == 31)
     setIBCTransfersCommonVar3(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
    			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, density, pressure, utau, yplus,
                              ipt_tmp, size,
                              gamma, cv, muS, Cs, Ts,
  			      vectOfDnrFields, vectOfRcvFields);
 
  } // Fin zone // omp


  // sortie
  RELEASESHAREDB(resr, arrayR, fr, cnr);
  RELEASESHAREDB(resd, arrayD, fd, cnd);
  BLOCKRELEASEMEM;
  BLOCKRELEASEMEM2;
  return tpl;
}
//=============================================================================
/* Effectue les transfers IBC en in-place */
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR, *zoneD;
  PyObject *pyVariables; 
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayXPC, *pyArrayXPI, *pyArrayXPW;
  PyObject *pyArrayYPC, *pyArrayYPI, *pyArrayYPW;
  PyObject *pyArrayZPC, *pyArrayZPI, *pyArrayZPW;
  PyObject *pyArrayDens, *pyArrayPressure, *pyArrayUtau, *pyArrayYplus ;
  E_Int bctype, loc, vartype, compact;
  E_Float gamma, cv, muS, Cs, Ts;

  if (!PYPARSETUPLE(args,
                    "OOOOOOOOOOOOOOOOOOOOlllldddddsss", "OOOOOOOOOOOOOOOOOOOOiiiidddddsss",
                    "OOOOOOOOOOOOOOOOOOOOllllfffffsss", "OOOOOOOOOOOOOOOOOOOOiiiifffffsss",
                    &zoneR, &zoneD, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, 
                    &pyArrayXPC, &pyArrayYPC, &pyArrayZPC,
                    &pyArrayXPW, &pyArrayYPW, &pyArrayZPW,
                    &pyArrayXPI, &pyArrayYPI, &pyArrayZPI, &pyArrayDens,  &pyArrayPressure,  &pyArrayUtau,  &pyArrayYplus,
                    &bctype    , &loc       , &vartype   , &compact   ,&gamma, &cv, &muS, &Cs, &Ts,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
  {
      return NULL;
  }

  vector<PyArrayObject*> hook;
  
  E_Int bcType = E_Int(bctype); // 0 : glissement; 1: adherence; 2: loi log; 3: loi de Musker


  /* varType : 
     1  : conservatives, 
     11 : conservatives + ronutildeSA 
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA 
     3  : (ro,u,v,w,p)     
     31 : (ro,u,v,w,p) + ronutildeSA   */
  E_Int nvars;
  E_Int varType = E_Int(vartype);
  if      ( varType ==  1 || varType ==  2 || varType ==  3 ) 
    nvars = 5;
  else if ( varType == 11 || varType == 21 || varType == 31 )
    nvars = 6;
  else 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "_setIBCTransfers:  varType value is not valid.");
    return NULL;
  }



  // recupere les champs du donneur (nodes)
  E_Int imdjmd, imd, jmd, kmd, cnNfldD, ndimdxR, ndimdxD, meshtype;;
  E_Float* iptroD; E_Float* iptroR;


# include "extract_interpD.h"
  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI, true);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();  

# include "IBC/extract_IBC.h"


  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Float*> vectOfDnrFields(nvars); vector<E_Float*> vectOfRcvFields(nvars);
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;
  //codage general (lent ;-) )
  if( compact==0)
    {// recupere les champs du donneur (nodes)
      E_Int cnSizeD;
      char* varStringD;
      vector<E_Int> locsD;
      vector<E_Int*> cnd;
      E_Int resd = K_PYTREE::getFromZone(zoneD, 0, 0, varStringD,
                                         fieldsD, locsD, imd, jmd, kmd,
                                         cnd, cnSizeD, cnNfldD, 
                                         eltTypeD, hook,
                                         GridCoordinates, 
                                         FlowSolutionNodes, FlowSolutionCenters);
      if (cnd.size() > 0) ptrcnd = cnd[0];

      meshtype = resd; // 1: structure, 2: non structure
      // recupere les champs du receveur 
      E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
      char* varStringR; vector<E_Int> locsR;
      vector<E_Int*> cnr;
      K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                            fieldsR, locsR, imr, jmr, kmr,
                            cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                            GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

      // -- no check (perfo) --
      // Transferts 
      // Types valides: 2, 3, 4, 5
      E_Int posvr, posvd;

      // Extrait les positions des variables a transferer
      E_Int nfoundvar = 0;
      if (PyList_Check(pyVariables) != 0)
      {
        int nvariables = PyList_Size(pyVariables);
        if (nvariables > 0)
        {
          for (int i = 0; i < nvariables; i++)
          {
            PyObject* tpl0 = PyList_GetItem(pyVariables, i);
            if (PyString_Check(tpl0) == 0)
              PyErr_Warn(PyExc_Warning, "_setIBCTransfers: variable must be a string. Skipped.");
            else 
            {
              char* varname = PyString_AsString(tpl0);        
              posvd = K_ARRAY::isNamePresent(varname, varStringD);      
              posvr = K_ARRAY::isNamePresent(varname, varStringR);      
              if (posvd != -1 && posvr != -1) 
              {
                
                vectOfRcvFields[nfoundvar]= fieldsR[posvr];
                vectOfDnrFields[nfoundvar]= fieldsD[posvd];
                //vectOfDnrFields.push_back(fieldsD[posvd]);
                //vectOfRcvFields.push_back(fieldsR[posvr]);
                nfoundvar += 1;
              }         
            }
          }
        }
      }
      delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

      if (nfoundvar != nvars)
      {
        RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
        BLOCKRELEASEMEM;
        BLOCKRELEASEMEM2;
        PyErr_SetString(PyExc_TypeError, 
                     "_setIBCTransfers: number of variables is inconsistent with varType.");
        return NULL;
      }
    }  

  // les variables a transferes sont compactes: on recuperes uniquement la premiere et la taille 
  else
    {
# include "getfromzonecompact.h"
      for (E_Int eq = 0; eq < nvars; eq++)
         {
            vectOfRcvFields[eq]= iptroR + eq*ndimdxR;
            vectOfDnrFields[eq]= iptroD + eq*ndimdxD;
         }
    }

////
////
////
//  Interpolation parallele
////  
////  
////  
# include "commonInterpTransfers_indirect.h"

    E_Int threadmax_sdm  = __NUMTHREADS__;

    E_Int size = (nbRcvPts/threadmax_sdm)+1; // on prend du gras pour gerer le residus
    E_Int    r =  size % 8;
    if (r != 0) size  = size + 8 - r;        // on rajoute du bas pour alignememnt 64bits
    //printf("size %d %d \n", size, r);
    if (bctype <=1 ) size = 0;               // tableau inutile

    FldArrayF  tmp(size*13*threadmax_sdm);
    E_Float* ipt_tmp=  tmp.begin();

#    pragma omp parallel default(shared)
     {

   //indice loop pour paralelisation omp
   E_Int ideb, ifin;
#ifdef _OPENMP
   E_Int  ithread           = omp_get_thread_num()+1;
   E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
   E_Int ithread = 1;
   E_Int Nbre_thread_actif = 1;
#endif
   // Calcul du nombre de champs a traiter par chaque thread
   E_Int chunk = nbRcvPts/Nbre_thread_actif;
   E_Int r = nbRcvPts - chunk*Nbre_thread_actif;
   // pts traitees par thread
   if (ithread <= r) 
        { ideb = (ithread-1)*(chunk+1); ifin = ideb + (chunk+1); }  
   else { ideb = (chunk+1)*r+(ithread-r-1)*chunk; ifin = ideb + chunk; }

  if (varType == 1 || varType == 11)
    setIBCTransfersCommonVar1(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, density, pressure, utau, yplus,
                              ipt_tmp, size,
                              gamma, cv, muS, Cs, Ts,
                              vectOfDnrFields, vectOfRcvFields);
  else if (varType == 2 || varType == 21)
    setIBCTransfersCommonVar2(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
	             	      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, density, pressure, utau, yplus,
                              ipt_tmp, size,
                              gamma, cv, muS, Cs, Ts,
			      vectOfDnrFields, vectOfRcvFields);
  else if (varType == 3 || varType == 31)
    setIBCTransfersCommonVar3(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, density, pressure, utau, yplus,
                              ipt_tmp, size,
                              gamma, cv, muS, Cs, Ts,
			      vectOfDnrFields, vectOfRcvFields);

     } // Fin zone // omp
  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  BLOCKRELEASEMEM;
  BLOCKRELEASEMEM2;
  Py_INCREF(Py_None);
  return Py_None;
}
