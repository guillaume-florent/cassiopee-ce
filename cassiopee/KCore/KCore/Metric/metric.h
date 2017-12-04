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

#ifndef _KCORE_METRIC_H
#define _KCORE_METRIC_H
# include "Def/DefTypes.h"
# include "Fld/FldArray.h"
# include "Def/DefFunction.h"

namespace K_METRIC
{
  /* Calcul les volumes pour des elements NGon
     IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
     IN: cn: connectivite NGon
     OUT: volp: pointeur sur le tableau des volumes calcules aux centres des elements
  */
  E_Int CompNGonVol(E_Float* xt, E_Float* yt, E_Float* zt, 
                    K_FLD::FldArrayI& cn, E_Float* volp); 
  
  E_Int CompNGonVolOfElement(
    E_Float* xt,E_Float* yt,E_Float* zt, 
    K_FLD::FldArrayI& cn, E_Int indE, std::vector<std::vector<E_Int> > cnEV, 
    K_FLD::FldArrayI& posElt, K_FLD::FldArrayI& posFace, 
    K_FLD::FldArrayI& dimElt, E_Float& vol);

  E_Int compNGonSurf(E_Float* xt, E_Float* yt, E_Float* zt, 
                     K_FLD::FldArrayI& cn, 
                     E_Float* sxp, E_Float* syp,  E_Float* szp); 
  
  /* Calcul des surfaces orientees des faces et la norme associee
     On suppose que le NGON est deja correctement oriente
     IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
     IN:  cnp: pointeur sur la connectivite NGon
     OUT: sxp, syp, szp, snp: surface orientee calculee pour les faces et norme associee
     Return 0 (OK), 1 (Failed)
  */
  E_Int compNGonFacesSurf(
    E_Float* xt, E_Float* yt, E_Float* zt, K_FLD::FldArrayI& cn,
    E_Float* sxp, E_Float* syp,  E_Float* szp, E_Float* snp, 
    K_FLD::FldArrayI* cFE=NULL);  

  /* Calcule l aire d'une cellule d un maillage surfacique nk=1. N'est pas n�cessairement dans le plan */
  E_Float compVolOfStructCell2D(E_Int ni, E_Int nj, E_Int indcell, E_Float* xt, E_Float* yt, E_Float* zt);
}
#endif
