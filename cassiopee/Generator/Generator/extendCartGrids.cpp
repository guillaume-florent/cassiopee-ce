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

#include "generator.h"

using namespace K_FLD;
using namespace std;
using namespace K_CONST;

//=============================================================================
/* Extend Cartesian grids w.r.t depth with a minimum overlapping */
//=============================================================================
PyObject* K_GENERATOR::extendCartGrids(PyObject* self, PyObject* args)
{
  PyObject *arrays;
  E_Int ext, optimized;
#ifdef E_DOUBLEINT 
  if (!PyArg_ParseTuple(args, "Oll", &arrays, &ext, &optimized)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Oii", &arrays, &ext, &optimized)) return NULL;
#endif
  if (ext < 0) 
  {
     PyErr_SetString(PyExc_TypeError, 
                    "extendCartGrids: ext must not be negative.");
     return NULL;
  }
  if (ext == 0) return arrays;
  if (optimized != 0 && optimized != 1)
  {printf("Warning: extendCartGrids: optimized is set to 1.\n"); optimized = 1;}

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit;
  vector<E_Int> njt; 
  vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true;
  E_Boolean skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured);
  if ( isOk == -1 ) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "extendCartGrids: invalid list of arrays.");
     return NULL;   
  }
  /* Verification des positions de x,y,z */
  E_Int nzones = structF.size();
  E_Int posxi, posyi, poszi;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  for (E_Int i = 0; i < nzones; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarString[i]);
    posyi = K_ARRAY::isCoordinateYPresent(structVarString[i]);
    poszi = K_ARRAY::isCoordinateZPresent(structVarString[i]);
    if ( posxi == -1 || posyi == -1 || poszi == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "extendCartGrids: arrays must contain coordinates.");
      K_ARRAY::cleanStructFields(structF); 
      return NULL;
    }
    posxi++; posyi++; poszi++;
    posxt.push_back(posxi); posyt.push_back(posyi); poszt.push_back(poszi);
  }
  E_Int dim = 3;
  if ( nkt[0] == 1 ) dim = 2;

  // determination des elts dont les bbox intersectent l elt courant
  FldArrayF bbox(nzones, 6);// xmin, ymin, zmin, xmax, ymax, zmax
  E_Float minB[3];  E_Float maxB[3];
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  E_Float tol = 1.e-6; E_Float tol2 = tol*tol;
  for (E_Int v = 0; v < nzones; v++)
  {
    K_COMPGEOM::boundingBox(
      nit[v], njt[v], nkt[v], posxt[v], posyt[v], poszt[v], 
      *structF[v], xminp[v], yminp[v], zminp[v], xmaxp[v], ymaxp[v], zmaxp[v]);
    minB[0] = xminp[v]; minB[1] = yminp[v]; minB[2] = zminp[v];
    maxB[0] = xmaxp[v]; maxB[1] = ymaxp[v]; maxB[2] = zmaxp[v];
    if (dim == 2) { minB[2] = 0.; maxB[2] = 1.; }
  }

  // Determination des extensions pour chq zone a partir de l'octree
  E_Int extg = ext; E_Int extf = ext;
  if ( optimized == 1 ) {extg = ext-1;}

  FldArrayI extension(nzones, 6); extension.setAllValuesAtNull();
  vector<E_Int> indicesBB;
  E_Int indA1, indB1, indC1, indD1, indE1, indF1, indG1, indH1;
  E_Int indA2, indB2, indC2, indD2, indE2, indF2, indG2, indH2, v2;
  E_Int* ext1 = extension.begin(1);
  E_Int* ext2 = extension.begin(2);
  E_Int* ext3 = extension.begin(3);
  E_Int* ext4 = extension.begin(4);
  E_Int* ext5 = extension.begin(5);
  E_Int* ext6 = extension.begin(6);
  E_Int ret;
  E_Float xp, yp, zp, dx, dy, dz;
  E_Float s1, s2; // surface de la facette initiale et surface de facette projetee
  E_Int nbboxes;
  E_Int ni1, nj1, nk1, ni2, nj2, nk2;
  E_Int shift1, shift2;
  E_Float *xt1, *yt1, *zt1, *xt2, *yt2, *zt2;
  E_Float dhmax;
  E_Int found1, found2, found3, found4;
  vector< vector<E_Int> > dejaVu(nzones); 
  if (dim == 2) 
  {
    FldArrayF face(2,3);//3D
    FldArrayI cnf(1,2); // cn BAR de la facette
    cnf(0,1) = 1; cnf(0,2) = 2;  
    for (E_Int v1 = 0; v1 < nzones; v1++)
    { 
      vector<E_Int>& dejaVu1 = dejaVu[v1];
      xt1 = structF[v1]->begin(posxt[v1]);
      yt1 = structF[v1]->begin(posyt[v1]);
      zt1 = structF[v1]->begin(poszt[v1]);
      ni1 = nit[v1]; nj1 = njt[v1]; nk1 = nkt[v1]; 
      indA1 = 0; indB1 = ni1-1; indD1 = (nj1-1)*ni1; indC1 = indB1 + indD1;

      /* facette AD ou i = 1 */
      s1 = (ymaxp[v1]-yminp[v1])/(njt[v1]-1);
      minB[0] = xt1[indA1]; minB[1] = yt1[indA1]; minB[2] = zt1[indA1];
      maxB[0] = xt1[indD1]; maxB[1] = yt1[indD1]; maxB[2] = zt1[indD1];
      indicesBB.clear(); getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
      // facette opposee en i = imax : B'C'

      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; 
        v2 = indicesBB[noe];
        if (v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end() ) 
          goto finimin2;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        face(0,1) = xt2[indB2]; face(0,2) = yt2[indB2]; face(0,3) = zt2[indB2];
        face(1,1) = xt2[indC2]; face(1,2) = yt2[indC2]; face(1,3) = zt2[indC2];

        if (K_FUNC::fEqualZero(xt1[indA1]-xt2[indB2],tol) == false ) goto finimin2;
        s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1); 
        
        // projeter le pt A sur la facette opposee
        ret = K_COMPGEOM::projectOrtho(xt1[indA1], yt1[indA1], zt1[indA1], 
                                       face.begin(1), face.begin(2), 
                                       face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indA1]; dy = yp-yt1[indA1]; dz = zp-zt1[indA1];
        if (ret > -1 && dx*dx + dy*dy + dz*dz <= tol2) 
        {found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2); }
        if (found1 == 1) goto finimin2;

        // projeter le pt D sur la facette opposee
        ret = K_COMPGEOM::projectOrtho(xt1[indD1], yt1[indD1], zt1[indD1], 
                                       face.begin(1), face.begin(2), 
                                       face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indD1]; dy = yp-yt1[indD1]; dz = zp-zt1[indD1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1) goto finimin2;

        finimin2:;
        if ( found1+found2 > 0 ) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);

          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext1[v1] = extf; ext2[v2] = K_FUNC::E_max(ext2[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext1[v1] = extf; ext2[v2] = K_FUNC::E_max(ext2[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext2[v2] = extf; ext1[v1] = K_FUNC::E_max(ext1[v1],extg);
          }
          goto faceimax2;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette i = 1
      faceimax2:;

      /* facette BC ou i = imax */
      minB[0] = xt1[indB1]; minB[1] = yt1[indB1]; minB[2] = zt1[indB1];
      maxB[0] = xt1[indC1]; maxB[1] = yt1[indC1]; maxB[2] = zt1[indC1];
      indicesBB.clear(); getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
      // facette opposee en i = imax : A'D'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; 
        v2 = indicesBB[noe]; 
        if ( v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end()) 
          goto finimax2;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        face(0,1) = xt2[indA2]; face(0,2) = yt2[indA2]; face(0,3) = zt2[indA2];
        face(1,1) = xt2[indD2]; face(1,2) = yt2[indD2]; face(1,3) = zt2[indD2];

        if ( K_FUNC::fEqualZero(xt1[indB1]-xt2[indA2],tol) == false ) goto finimax2;
        s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1); 
        
        // projeter le pt B sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indB1], yt1[indB1], zt1[indB1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indB1]; dy = yp-yt1[indB1]; dz = zp-zt1[indB1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        {found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2); }
        if ( found1 == 1) goto finimax2;

        // projeter le pt C sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indC1], yt1[indC1], zt1[indC1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indC1]; dy = yp-yt1[indC1]; dz = zp-zt1[indC1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1) goto finimax2;

        finimax2:;
        if ( found1+found2 > 0) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext2[v1] = extf; ext1[v2] = K_FUNC::E_max(ext1[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext2[v1] = extf; ext1[v2] = K_FUNC::E_max(ext1[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext1[v2] = extf; ext2[v1] = K_FUNC::E_max(ext2[v1],extg);
          }
          goto facejmin2;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette i = imax      

      facejmin2:;
      /* facette AB ou i = 1 */
      s1 = (xmaxp[v1]-xminp[v1])/(nit[v1]-1);
      minB[0] = xt1[indA1]; minB[1] = yt1[indA1]; minB[2] = zt1[indA1];
      maxB[0] = xt1[indB1]; maxB[1] = yt1[indB1]; maxB[2] = zt1[indB1];
      indicesBB.clear(); getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
      // facette opposee en i = imax : C'D'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; 
        v2 = indicesBB[noe]; 
        if ( v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end()) 
          goto finjmin2;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        face(0,1) = xt2[indC2]; face(0,2) = yt2[indC2]; face(0,3) = zt2[indC2];
        face(1,1) = xt2[indD2]; face(1,2) = yt2[indD2]; face(1,3) = zt2[indD2];

        if ( K_FUNC::fEqualZero(yt1[indA1]-yt2[indC2],tol) == false ) goto finjmin2;
        s2 = (xmaxp[v2]-xminp[v2])/(nit[v2]-1); 
        
        // projeter le pt A sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indA1], yt1[indA1], zt1[indA1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indA1]; dy = yp-yt1[indA1]; dz = zp-zt1[indA1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        {found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2); }
        if ( found1 == 1 ) goto finjmin2;

        // projeter le pt B sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indB1], yt1[indB1], zt1[indB1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indB1]; dy = yp-yt1[indB1]; dz = zp-zt1[indB1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1) goto finjmin2;
        finjmin2:;

        if ( found1+found2 > 0 ) 
        {
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext3[v1] = extf; ext4[v2] = K_FUNC::E_max(ext4[v2],extg);
          }
          else if (s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext3[v1] = extf; ext4[v2] = K_FUNC::E_max(ext4[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext4[v2] = extf; ext3[v1] = K_FUNC::E_max(ext3[v1],extg);
          }
          goto facejmax2;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette i = 1

      facejmax2:;      
      /* facette DC ou i = imax */
      minB[0] = xt1[indD1]; minB[1] = yt1[indD1]; minB[2] = zt1[indD1];
      maxB[0] = xt1[indC1]; maxB[1] = yt1[indC1]; maxB[2] = zt1[indC1];
      indicesBB.clear(); getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
      // facette opposee en i = imax : A'B'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; 
        v2 = indicesBB[noe]; 
        if (v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end() ) 
          goto finjmax2;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        face(0,1) = xt2[indA2]; face(0,2) = yt2[indA2]; face(0,3) = zt2[indA2];
        face(1,1) = xt2[indB2]; face(1,2) = yt2[indB2]; face(1,3) = zt2[indB2];

        if ( K_FUNC::fEqualZero(yt1[indC1]-yt2[indA2],tol) == false ) goto finjmax2;
        s2 = (xmaxp[v2]-xminp[v2])/(nit[v2]-1); 
        
        // projeter le pt B sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indC1], yt1[indC1], zt1[indC1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indC1]; dy = yp-yt1[indC1]; dz = zp-zt1[indC1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        {found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2); }
        if ( found1 == 1 ) goto finjmax2;

        // projeter le pt C sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indD1], yt1[indD1], zt1[indD1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indD1]; dy = yp-yt1[indD1]; dz = zp-zt1[indD1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1) goto finjmax2;

        finjmax2:;
        if ( found1+found2 > 0 ) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext4[v1] = extf; ext3[v2] = K_FUNC::E_max(ext3[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext4[v1] = extf; ext3[v2] = K_FUNC::E_max(ext3[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext3[v2] = extf; ext4[v1] = K_FUNC::E_max(ext4[v1],extg);
          }
          goto end2;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette j = jmax      
      end2:;
    }
  }
  else //( dim == 3 )
  {    
    FldArrayF face(4,3);//facette = 2 TRI
    FldArrayI cnf(2,3);// cn TRI de la facette de l HEXA
    cnf(0,1) = 1; cnf(0,2) = 2; cnf(0,3) = 4;
    cnf(1,1) = 2; cnf(1,2) = 3; cnf(1,3) = 4;

    for (E_Int v1 = 0; v1 < nzones; v1++)
    {
      xt1 = structF[v1]->begin(posxt[v1]);
      yt1 = structF[v1]->begin(posyt[v1]);
      zt1 = structF[v1]->begin(poszt[v1]);
      ni1 = nit[v1]; nj1 = njt[v1]; nk1 = nkt[v1]; shift1 = (nk1-1)*ni1*nj1;
      indA1 = 0; indB1 = ni1-1; indD1 = (nj1-1)*ni1; indC1 = indB1 + indD1;
      indE1 = indA1+shift1; indF1 = indB1+shift1; indH1 = indD1+shift1; indG1 = indC1+shift1;
      s1 = (ymaxp[v1]-yminp[v1])/(njt[v1]-1);

      /* facette ADHE ou i = 1 */
      minB[0] = xt1[indA1]; minB[1] = yt1[indA1]; minB[2] = zt1[indA1];
      maxB[0] = xt1[indH1]; maxB[1] = yt1[indH1]; maxB[2] = zt1[indH1];

      // elts intersectant la facette i = 1
      indicesBB.clear();
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
      vector<E_Int>& dejaVu1 = dejaVu[v1];

      // facette opposee en i = imax : B'C'G'F'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; found3 = 0; found4 = 0;
        v2 = indicesBB[noe]; 
        if (v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end() ) 
          goto finimin;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        indE2 = indA2+shift2; indF2 = indB2+shift2; indH2 = indD2+shift2; indG2 = indC2+shift2;
        face(0,1) = xt2[indB2]; face(0,2) = yt2[indB2]; face(0,3) = zt2[indB2];
        face(1,1) = xt2[indC2]; face(1,2) = yt2[indC2]; face(1,3) = zt2[indC2];
        face(2,1) = xt2[indG2]; face(2,2) = yt2[indG2]; face(2,3) = zt2[indG2];
        face(3,1) = xt2[indF2]; face(3,2) = yt2[indF2]; face(3,3) = zt2[indF2];

        if ( K_FUNC::fEqualZero(xt1[indA1]-xt2[indB2],tol) == false ) goto finimin;
        s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1); 
        
        // projeter le pt A sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indA1], yt1[indA1], zt1[indA1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indA1]; dy = yp-yt1[indA1]; dz = zp-zt1[indA1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        {found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2); }
        if ( found1 == 1 ) goto finimin;

        // projeter le pt D sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indD1], yt1[indD1], zt1[indD1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indD1]; dy = yp-yt1[indD1]; dz = zp-zt1[indD1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1 ) goto finimin;

        // projeter le pt H sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indH1], yt1[indH1], zt1[indH1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indH1]; dy = yp-yt1[indH1]; dz = zp-zt1[indH1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found3 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found3 == 1 ) goto finimin;

        // projeter le pt E sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indE1], yt1[indE1], zt1[indE1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indE1]; dy = yp-yt1[indE1]; dz = zp-zt1[indE1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 )
        { found4 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found4 == 1 ) goto finimin;

        finimin:;
        if ( found1+found2+found3+found4 > 0) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext1[v1] = extf; ext2[v2] = K_FUNC::E_max(ext2[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext1[v1] = extf; ext2[v2] = K_FUNC::E_max(ext2[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext2[v2] = extf; ext1[v1] = K_FUNC::E_max(ext1[v1],extg);
          }
          goto faceimax;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette i = 1

      /* facette BCGF ou i = imax */
      faceimax:;
      minB[0] = xt1[indB1]; minB[1] = yt1[indB1]; minB[2] = zt1[indB1];
      maxB[0] = xt1[indG1]; maxB[1] = yt1[indG1]; maxB[2] = zt1[indG1];

      // elts intersectant la facette i = imax
      indicesBB.clear(); 
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes

      // facette opposee en i = 1: A'D'H'E'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; found3 = 0; found4 = 0;
        v2 = indicesBB[noe]; 
        if ( v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end()) 
          goto finimax;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        indE2 = indA2+shift2; indF2 = indB2+shift2; indH2 = indD2+shift2; indG2 = indC2+shift2;
        face(0,1) = xt2[indA2]; face(0,2) = yt2[indA2]; face(0,3) = zt2[indA2];
        face(1,1) = xt2[indD2]; face(1,2) = yt2[indD2]; face(1,3) = zt2[indD2];
        face(2,1) = xt2[indH2]; face(2,2) = yt2[indH2]; face(2,3) = zt2[indH2];
        face(3,1) = xt2[indE2]; face(3,2) = yt2[indE2]; face(3,3) = zt2[indE2];

        if ( K_FUNC::fEqualZero(xt1[indB1]-xt2[indA2],tol) == false ) goto finimax;
        s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1); 
        
        // projeter le pt B sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indB1], yt1[indB1], zt1[indB1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indB1]; dy = yp-yt1[indB1]; dz = zp-zt1[indB1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found1 == 1 ) goto finimax;

        // projeter le pt C sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indC1], yt1[indC1], zt1[indC1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indC1]; dy = yp-yt1[indC1]; dz = zp-zt1[indC1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1 ) goto finimax;

        // projeter le pt G sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indG1], yt1[indG1], zt1[indG1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indG1]; dy = yp-yt1[indG1]; dz = zp-zt1[indG1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found3 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found3 == 1 ) goto finimax;

        // projeter le pt F sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indF1], yt1[indF1], zt1[indF1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indF1]; dy = yp-yt1[indF1]; dz = zp-zt1[indF1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 )
        { found4 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found4 == 1 ) goto finimax;

        finimax:;
        if ( found1+found2+found3+found4 > 0) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext2[v1] = extf; ext1[v2] = K_FUNC::E_max(ext1[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext2[v1] = extf; ext1[v2] = K_FUNC::E_max(ext1[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext1[v2] = extf; ext2[v1] = K_FUNC::E_max(ext2[v1],extg);
          }
          goto facejmin;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette i = imax

      s1 = (xmaxp[v1]-xminp[v1])/(nit[v1]-1);
      /* facette j =jmin ou ABFE */
      facejmin:;    
      minB[0] = xt1[indA1]; minB[1] = yt1[indA1]; minB[2] = zt1[indA1];
      maxB[0] = xt1[indF1]; maxB[1] = yt1[indF1]; maxB[2] = zt1[indF1];

      // elts intersectant la facette j=1
      indicesBB.clear(); 
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes

      // facette opposee en j = jmax: D'C'G'H'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; found3 = 0; found4 = 0;
        v2 = indicesBB[noe]; 
        if ( v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end()) 
          goto finjmin;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        indE2 = indA2+shift2; indF2 = indB2+shift2; indH2 = indD2+shift2; indG2 = indC2+shift2;
        face(0,1) = xt2[indD2]; face(0,2) = yt2[indD2]; face(0,3) = zt2[indD2];
        face(1,1) = xt2[indC2]; face(1,2) = yt2[indC2]; face(1,3) = zt2[indC2];
        face(2,1) = xt2[indG2]; face(2,2) = yt2[indG2]; face(2,3) = zt2[indG2];
        face(3,1) = xt2[indH2]; face(3,2) = yt2[indH2]; face(3,3) = zt2[indH2];

        if ( K_FUNC::fEqualZero(yt1[indA1]-yt2[indD2],tol) == false ) goto finjmin;
        s2 = (xmaxp[v2]-xminp[v2])/(nit[v2]-1); 
        
        // projeter le pt A sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indA1], yt1[indA1], zt1[indA1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indA1]; dy = yp-yt1[indA1]; dz = zp-zt1[indA1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found1 == 1 ) goto finjmin;

        // projeter le pt B sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indB1], yt1[indB1], zt1[indB1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indB1]; dy = yp-yt1[indB1]; dz = zp-zt1[indB1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1 ) goto finjmin;

        // projeter le pt F sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indF1], yt1[indF1], zt1[indF1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indF1]; dy = yp-yt1[indF1]; dz = zp-zt1[indF1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found3 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found3 == 1 ) goto finjmin;

        // projeter le pt E sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indE1], yt1[indE1], zt1[indE1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indE1]; dy = yp-yt1[indE1]; dz = zp-zt1[indE1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 )
        { found4 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found4 == 1 ) goto finjmin;

        finjmin:;
        if ( found1+found2+found3+found4 > 0 ) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext3[v1] = extf; ext4[v2] = K_FUNC::E_max(ext4[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext3[v1] = extf; ext4[v2] = K_FUNC::E_max(ext4[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext4[v2] = extf; ext3[v1] = K_FUNC::E_max(ext3[v1],extg);
          }
          goto facejmax;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette j=1
      
      facejmax:;
      //facette j =jmax ou DCGH
      minB[0] = xt1[indD1]; minB[1] = yt1[indD1]; minB[2] = zt1[indD1];
      maxB[0] = xt1[indG1]; maxB[1] = yt1[indG1]; maxB[2] = zt1[indG1];

      // elts intersectant la facette j=jmax
      indicesBB.clear(); 
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes

      // facette opposee en i = 1: A'B'F'E'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; found3 = 0; found4 = 0;
        v2 = indicesBB[noe]; 
        if ( v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end()) 
         goto finjmax;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        indE2 = indA2+shift2; indF2 = indB2+shift2; indH2 = indD2+shift2; indG2 = indC2+shift2;
        face(0,1) = xt2[indA2]; face(0,2) = yt2[indA2]; face(0,3) = zt2[indA2];
        face(1,1) = xt2[indB2]; face(1,2) = yt2[indB2]; face(1,3) = zt2[indB2];
        face(2,1) = xt2[indF2]; face(2,2) = yt2[indF2]; face(2,3) = zt2[indF2];
        face(3,1) = xt2[indE2]; face(3,2) = yt2[indE2]; face(3,3) = zt2[indE2];

        if ( K_FUNC::fEqualZero(yt1[indD1]-yt2[indA2],tol) == false ) goto finjmax;
        s2 = (xmaxp[v2]-xminp[v2])/(nit[v2]-1); 

        // projeter le pt D sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indD1], yt1[indD1], zt1[indD1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indD1]; dy = yp-yt1[indD1]; dz = zp-zt1[indD1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found1 == 1 ) goto finjmax;

        // projeter le pt C sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indC1], yt1[indC1], zt1[indC1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indC1]; dy = yp-yt1[indC1]; dz = zp-zt1[indC1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1 ) goto finjmax;

        // projeter le pt G sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indG1], yt1[indG1], zt1[indG1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indG1]; dy = yp-yt1[indG1]; dz = zp-zt1[indG1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found3 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found3 == 1 ) goto finjmax;

        // projeter le pt H sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indH1], yt1[indH1], zt1[indH1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indH1]; dy = yp-yt1[indH1]; dz = zp-zt1[indH1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 )
        { found4 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found4 == 1 ) goto finjmax;

        finjmax:;
        if ( found1+found2+found3+found4 > 0) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext4[v1] = extf; ext3[v2] = K_FUNC::E_max(ext3[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext4[v1] = extf; ext3[v2] = K_FUNC::E_max(ext3[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext3[v2] = extf; ext4[v1] = K_FUNC::E_max(ext4[v1],extg);
          }
          goto facekmin;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette j=jmax

      facekmin:; 
      /* facette ABCD */
      minB[0] = xt1[indA1]; minB[1] = yt1[indA1]; minB[2] = zt1[indA1];
      maxB[0] = xt1[indC1]; maxB[1] = yt1[indC1]; maxB[2] = zt1[indC1];
      // elts intersectant la facette k = 1
      indicesBB.clear();
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
       // facette opposee en i = imax : E'F'G'H'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; found3 = 0; found4 = 0;
        v2 = indicesBB[noe]; 
        if ( v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end()) 
          goto finkmin;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        indE2 = indA2+shift2; indF2 = indB2+shift2; indH2 = indD2+shift2; indG2 = indC2+shift2;
        face(0,1) = xt2[indE2]; face(0,2) = yt2[indE2]; face(0,3) = zt2[indE2];
        face(1,1) = xt2[indF2]; face(1,2) = yt2[indF2]; face(1,3) = zt2[indF2];
        face(2,1) = xt2[indG2]; face(2,2) = yt2[indG2]; face(2,3) = zt2[indG2];
        face(3,1) = xt2[indH2]; face(3,2) = yt2[indH2]; face(3,3) = zt2[indH2];

        if ( K_FUNC::fEqualZero(zt1[indA1]-zt2[indE2],tol) == false ) goto finkmin;
        s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1); 
        
        // projeter le pt A sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indA1], yt1[indA1], zt1[indA1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indA1]; dy = yp-yt1[indA1]; dz = zp-zt1[indA1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        {found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2); }
        if ( found1 == 1 ) goto finkmin;

        // projeter le pt B sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indB1], yt1[indB1], zt1[indB1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indB1]; dy = yp-yt1[indB1]; dz = zp-zt1[indB1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1 ) goto finkmin;

        // projeter le pt C sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indC1], yt1[indC1], zt1[indC1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indC1]; dy = yp-yt1[indC1]; dz = zp-zt1[indC1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found3 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found3 == 1 ) goto finkmin;

        // projeter le pt D sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indD1], yt1[indD1], zt1[indD1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indD1]; dy = yp-yt1[indD1]; dz = zp-zt1[indD1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 )
        { found4 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found4 == 1 ) goto finkmin;

        finkmin:;
        if ( found1+found2+found3+found4 > 0 ) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext5[v1] = extf; ext6[v2] = K_FUNC::E_max(ext6[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext5[v1] = extf; ext6[v2] = K_FUNC::E_max(ext6[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext6[v2] = extf; ext5[v1] = K_FUNC::E_max(ext5[v1],extg);
          }
          goto facekmax;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette k = 1     

      facekmax:; 
      /* facette EFGH */
      minB[0] = xt1[indE1]; minB[1] = yt1[indE1]; minB[2] = zt1[indE1];
      maxB[0] = xt1[indG1]; maxB[1] = yt1[indG1]; maxB[2] = zt1[indG1];
      // elts intersectant la facette k = kmax
      indicesBB.clear();
      getBlocksIntersecting(v1, minB, maxB, bbox, tol, indicesBB);
      nbboxes = indicesBB.size();
      dhmax = 0.;// dh max des grilles adjacentes
       // facette opposee en i = imax : A'B'C'D'
      for (E_Int noe = 0; noe < nbboxes; noe++)
      {
        found1 = 0; found2 = 0; found3 = 0; found4 = 0;
        v2 = indicesBB[noe]; 
        if ( v2 == v1 || find(dejaVu1.begin(), dejaVu1.end(), v2) != dejaVu1.end()) 
          goto finkmax;
        xt2 = structF[v2]->begin(posxt[v2]);
        yt2 = structF[v2]->begin(posyt[v2]);
        zt2 = structF[v2]->begin(poszt[v2]);
        ni2 = nit[v2]; nj2 = njt[v2]; nk2 = nkt[v2];
        shift2 = (nk2-1)*ni2*nj2;
        indA2 = 0; indB2 = ni2-1; indD2 = (nj2-1)*ni2; indC2 = indB2 + indD2;
        indE2 = indA2+shift2; indF2 = indB2+shift2; indH2 = indD2+shift2; indG2 = indC2+shift2;
        face(0,1) = xt2[indA2]; face(0,2) = yt2[indA2]; face(0,3) = zt2[indA2];
        face(1,1) = xt2[indB2]; face(1,2) = yt2[indB2]; face(1,3) = zt2[indB2];
        face(2,1) = xt2[indC2]; face(2,2) = yt2[indC2]; face(2,3) = zt2[indC2];
        face(3,1) = xt2[indD2]; face(3,2) = yt2[indD2]; face(3,3) = zt2[indD2];

        if ( K_FUNC::fEqualZero(zt1[indE1]-zt2[indA2],tol) == false ) goto finkmax;
        s2 = (ymaxp[v2]-yminp[v2])/(njt[v2]-1); 
        
        // projeter le pt E sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indE1], yt1[indE1], zt1[indE1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indE1]; dy = yp-yt1[indE1]; dz = zp-zt1[indE1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        {found1 = 1; dhmax = K_FUNC::E_max(dhmax,s2); }
        if ( found1 == 1 ) goto finkmax;

        // projeter le pt F sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indF1], yt1[indF1], zt1[indF1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indF1]; dy = yp-yt1[indF1]; dz = zp-zt1[indF1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found2 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found2 == 1 ) goto finkmax;

        // projeter le pt G sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indG1], yt1[indG1], zt1[indG1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indG1]; dy = yp-yt1[indG1]; dz = zp-zt1[indG1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 ) 
        { found3 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found3 == 1 ) goto finkmax;

        // projeter le pt H sur la facette opposee
        ret = K_COMPGEOM::projectOrtho( xt1[indH1], yt1[indH1], zt1[indH1], 
                                        face.begin(1), face.begin(2), 
                                        face.begin(3), cnf, xp, yp, zp);
        dx = xp-xt1[indH1]; dy = yp-yt1[indH1]; dz = zp-zt1[indH1];
        if ( ret > -1 && dx*dx + dy*dy + dz*dz <= tol2 )
        { found4 = 1; dhmax = K_FUNC::E_max(dhmax,s2);}
        if ( found4 == 1 ) goto finkmax;

        finkmax:;
        if ( found1+found2+found3+found4 > 0 ) 
        { 
          vector<E_Int>& dejaVu1 = dejaVu[v1]; dejaVu1.push_back(v2);
          vector<E_Int>& dejaVu2 = dejaVu[v2]; dejaVu2.push_back(v1);
          if ( K_FUNC::fEqualZero(s1-dhmax,tol2) == true ) 
          {
            ext6[v1] = extf; ext5[v2] = K_FUNC::E_max(ext5[v2],extg);
          }
          else if ( s1 < dhmax - tol2) // niveau fin etendu de ext+1
          {
            ext6[v1] = extf; ext5[v2] = K_FUNC::E_max(ext5[v2],extg);
          }
          else // niveau fin etendu de ext+1
          {
            ext5[v2] = extf; ext6[v1] = K_FUNC::E_max(ext6[v1],extg);
          }
          goto end;
        }
      }// fin parcours de ts les elts intersectant 
      // fin test facette k = kmax           
      end:;
    }// pour ts les elts
  }

  E_Int nio, njo, nko, ni, nj,nk, npts, ind;
 
  for (E_Int v = 0; v < nzones; v++)
  {
    ni = nit[v]; nj = njt[v]; nk = nkt[v]; 
    E_Float* xp = structF[v]->begin(posxt[v]);
    E_Float* yp = structF[v]->begin(posyt[v]);
    E_Float* zp = structF[v]->begin(poszt[v]);
    E_Float dh = xp[1]-xp[0];
    E_Float xxor = xp[0]-ext1[v]*dh;
    E_Float yyor = yp[0]-ext3[v]*dh;
    E_Float zzor = zp[0]-ext5[v]*dh;
    nio = ni+ext1[v]+ext2[v]; njo = nj+ext3[v]+ext4[v]; nko = nk+ext5[v]+ext6[v];
    npts = nio*njo*nko;
    FldArrayF* newcoords = new FldArrayF(npts,3);
    E_Float* xn = newcoords->begin(1);
    E_Float* yn = newcoords->begin(2);
    E_Float* zn = newcoords->begin(3);
    E_Int nionjo = nio*njo;
    
    for (E_Int k = 0; k < nko; k++)    
      for (E_Int j = 0; j < njo; j++)
        for (E_Int i = 0; i < nio; i++)
        {
          ind = i + j*nio + k*nionjo; 
          xn[ind] = xxor + i*dh;
          yn[ind] = yyor + j*dh;
          zn[ind] = zzor + k*dh;
        }
    delete structF[v];
    nit[v] = nio; njt[v] = njo; nkt[v] = nko; structF[v] = newcoords;  
  }
  
  PyObject* l = PyList_New(0); 
  for (E_Int v = 0; v < nzones; v++)
  {
    PyObject* tpl = K_ARRAY::buildArray(*structF[v], "x,y,z", nit[v], njt[v], nkt[v]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  //nettoyage
  for (E_Int v = 0; v < nzones; v++) delete structF[v];
  return l;
}

//=============================================================================
/* Intersection des bbox des elements */
//=============================================================================
void K_GENERATOR::getBlocksIntersecting(E_Int noz1, 
                                        E_Float* minB, E_Float* maxB,
                                        FldArrayF& bbox, E_Float tol, 
                                        vector<E_Int>& listOfZones)
{
  listOfZones.clear();
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  E_Float maxB0 = maxB[0]+tol;
  E_Float minB0 = minB[0]-tol;
  E_Float maxB1 = maxB[1]+tol;
  E_Float minB1 = minB[1]-tol;
  E_Float maxB2 = maxB[2]+tol;
  E_Float minB2 = minB[2]-tol;

  for (E_Int noz = 0; noz < bbox.getSize(); noz++)
  {
    if ( noz != noz1 ) 
    {
       if ( xminp[noz] <= maxB0 && xmaxp[noz] >= minB0 &&
            yminp[noz] <= maxB1 && ymaxp[noz] >= minB1 &&
            zminp[noz] <= maxB2 && zmaxp[noz] >= minB2 ) 
         listOfZones.push_back(noz);
    }
  }
}
