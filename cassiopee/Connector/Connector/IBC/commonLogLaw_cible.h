     err = 1; skip =0;
     //initialisation Newton SA  + vitesse cible
#ifdef _OPENMP4
     #pragma omp simd
#endif 
     for (E_Int noind = 0; noind < ifin-ideb; noind++)
       {

        yplus   = utau_vec[noind]*yplus_vec[noind];
        yplus_vec[noind] = yplus ;//* yext/yibc;//yplus au pt image

        umod = utau_vec[noind]*( kappainv*log(yplus) + cc );// U au pt IBC ;
        //umod = K_FUNC::E_abs(umod);

        ucible0 = sign_vec[noind] * umod;

        ucible_vec[noind] += ucible0 * ut_vec[noind]; // vitesse tangentielle pour le pt IBC
        vcible_vec[noind] += ucible0 * vt_vec[noind];
        wcible_vec[noind] += ucible0 * wt_vec[noind];


        nutcible_vec[noind] = kappa * alpha_vec[noind]* utau_vec[noind]; //negatif si pt ibc interieur aux corps

        nutilde             = K_FUNC::E_abs(  nutcible_vec[noind] );
        aa_vec[noind]       = nutilde;
#       include "IBC/fnutilde_vec.h"
        // out= mut  et err 
        // utau_vec[noind]    = nutilde; // uatu ne peut plus etre ecraser, car variable de post
         ut_vec[noind]    = nutilde; // Sauvegarde du nutilde, au cas ou newton non convergent.
       }
