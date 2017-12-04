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

#include "connector.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Search for the fringe of interpolated nodes near blanked points; depth is 
   the number of layers of interpolated nodes. 
   IN: blankedCells: -1, point masque, 0: point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedNodesUnstr(
  E_Int depth,  FldArrayI& cnEV,
  FldArrayI& blankedCells,
  FldArrayI& cellN)
{
  E_Int nvert = blankedCells.getSize();
  std::vector< std::vector<E_Int> > cVN(nvert);
  K_CONNECT::connectEV2VNbrs(cnEV, cVN);
  E_Int nvoisins;

  for (E_Int ind = 0; ind < nvert; ind++)
  {
    if (blankedCells[ind]  == -1) 
    {
      std::vector<E_Int>& voisins = cVN[ind];
      nvoisins = voisins.size();
      for (E_Int nov = 0; nov < nvoisins; nov++)
      {
        E_Int indv = voisins[nov]-1;
        cellN[indv] = K_FUNC::E_min(0,cellN[indv]);
      }
    }
  }
  FldArrayIS tag(nvert); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int ind = 0; ind<nvert; ind++)
    {
      if (cellN[ind] == 0)// pt interpole 
      {
        std::vector<E_Int>& voisins = cVN[ind];
        nvoisins = voisins.size();
        for (E_Int nov = 0; nov < nvoisins; nov++)
        {
          E_Int indv = voisins[nov]-1;
          if (cellN[indv] == 1) tag[indv] = 1;
        }
      }
    }
  }
  for (E_Int ind = 0; ind < nvert; ind++)
  {if (tag[ind] == 1) cellN[ind] = 0;}
  return;
}

//=============================================================================
/* Search for the fringe of interpolated cells near blanked points; depth is 
   the number of layers of interpolated cells. 
   IN: blankedCells: -1, point masque, 0 : point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsNGON(E_Int depth, FldArrayI& cNG,
                                                  FldArrayI& blankedCells,
                                                  FldArrayI& cellN)
{
  FldArrayI cFE;
  E_Int* cnp = cNG.begin();       
  E_Int sizeFN = cnp[1];         // taille de la connectivite face/noeuds
  E_Int nelts = cnp[sizeFN+2];         // nombre d elements       
  std::vector< std::vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectNG2FE(cNG, cFE);
  K_CONNECT::connectFE2EENbrs(cFE, cEEN);
  E_Int nvoisins;

  //1st layer, depth = 1
  for (E_Int et = 0; et < nelts; et++)
  {
    if (blankedCells[et] == -1)// pt masque 
    {
      std::vector<E_Int>& voisins = cEEN[et];
      nvoisins = voisins.size();
      for (E_Int noev = 0; noev < nvoisins; noev++)
      {
        E_Int et2 = voisins[noev];
        cellN[et2] = K_FUNC::E_min(0,cellN[et2]);
      }
    }
  }

  FldArrayIS tag(nelts); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      if (cellN[et] == 0)// pt interpole 
      {
        std::vector<E_Int>& voisins = cEEN[et];
        nvoisins = voisins.size();
        for (E_Int noev = 0; noev < nvoisins; noev++)
        {
          E_Int et2 = voisins[noev];
          if (cellN[et2] == 1) tag[et2] = 1;
        }
      }
    }
  }
  for (E_Int et = 0; et < nelts; et++)
  { if ( tag[et] == 1) cellN[et] = 0; }
}
//=============================================================================
/* Search for the fringe of interpolated cells near blanked points; depth is 
   the number of layers of interpolated cells. 
   IN: blankedCells: -1, point masque, 0 : point interpole, 1, point normal.
   IN/OUT: cellN: -1, point masque, 0, point interpole, 1, point normal.*/
//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsUnstr(char* eltType, 
                                                   E_Int depth, FldArrayI& cnEV,
                                                   FldArrayI& blankedCells,
                                                   FldArrayI& cellN)
{
  E_Int nelts = cnEV.getSize();
  E_Int nvert = nelts*cnEV.getNfld();
  std::vector< std::vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs(eltType, nvert, cnEV, cEEN); 
                       
  E_Int nvoisins;

  //1st layer, depth = 1
  for (E_Int et = 0; et < nelts; et++)
  {
    if (blankedCells[et] == -1)// pt masque 
    {
      std::vector<E_Int>& voisins = cEEN[et];
      nvoisins = voisins.size();
      for (E_Int noev = 0; noev < nvoisins; noev++)
      {
        E_Int et2 = voisins[noev];
        cellN[et2] = K_FUNC::E_min(0,cellN[et2]);
      }
    }
  }

  FldArrayIS tag(nelts); tag.setAllValuesAtNull();
  for (E_Int d = 2; d<= depth; d++)
  {
    for (E_Int et = 0; et < nelts; et++)
    {
      if (cellN[et] == 0)// pt interpole 
      {
        std::vector<E_Int>& voisins = cEEN[et];
        nvoisins = voisins.size();
        for (E_Int noev = 0; noev < nvoisins; noev++)
        {
          E_Int et2 = voisins[noev];
          if (cellN[et2] == 1) tag[et2] = 1;
        }
      }
    }
  }
  for (E_Int et = 0; et < nelts; et++)
  { if ( tag[et] == 1) cellN[et] = 0; }
}

//=============================================================================
void K_CONNECTOR::searchMaskInterpolatedCellsStruct(E_Int imc, E_Int jmc, E_Int kmc, E_Int depth, E_Int dir,
                                                    FldArrayI& blankedCells, FldArrayI& cellN)
{
  E_Int imjmc = imc*jmc;
  E_Int imjmkmc = imjmc*kmc;
  E_Int i, j, k, sensor, unmsensor, ind2;
  E_Int im1, ip1, jm1, jp1, km1, kp1;
  E_Int km1imjmc, kimjmc, kp1imjmc;
  E_Int nindices;
  //On n etend que les points masques (blankedcells = -1)
  if (dir == 0) //directionnel
  {
    if (kmc == 1) 
    {
      nindices = 4;      
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          j = ind/imc;
          i = ind-j*imc;
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;
              
          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          indices[0] = im1 + j*imc;
          indices[1] = ip1 + j*imc;
          indices[2] = i + jm1*imc;
          indices[3] = i + jp1*imc;      
          
          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }        
      }
    }// fin 2D
    else 
    {
      nindices = 6;      
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          k = ind/imjmc;
          j = ( ind-k*imjmc )/imc;
          i = ind-k*imjmc-j*imc;  
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;
              
          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          km1 = K_FUNC::E_max(0,k-d); kp1 = K_FUNC::E_min(k+d,kmc-1);

          indices[0] = im1 + j*imc + k*imjmc;
          indices[1] = ip1 + j*imc + k*imjmc;
          indices[2] = i + jm1*imc + k*imjmc;
          indices[3] = i + jp1*imc + k*imjmc;      
          indices[4] = i + j*imc + km1*imjmc;      
          indices[5] = i + j*imc + kp1*imjmc;

          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }        
      }
    }//fin 3D dir = 0
  }//dir = 0
  else 
  {
    if (kmc == 1) 
    {
      nindices = 8;      
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmc; ind++)
        {
          j = ind/imc;
          i = ind-j*imc;
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;
              
          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          indices[0] = im1 + jm1*imc;
          indices[1] = i + jm1*imc;
          indices[2] = ip1 + jm1*imc;
          
          indices[3] = im1 + j*imc;
          indices[4] = ip1 + j*imc;

          indices[5] = im1 + jp1*imc;
          indices[6] = i  +  jp1*imc;
          indices[7] = ip1 + jp1*imc;
          
          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }        
      }
    }// 2D dir = 1
    else // 3D 
    {
      nindices = 26;
      vector<E_Int> indices(nindices);
      for (E_Int d = 1; d <= depth; d++)
      {
        for (E_Int ind = 0; ind < imjmkmc; ind++)
        {
          k = ind/imjmc;
          j = ( ind-k*imjmc )/imc;
          i = ind-k*imjmc-j*imc;  
          sensor = (2+blankedCells[ind])/2;
          unmsensor = 1-sensor;
          
          im1 = K_FUNC::E_max(0,i-d); ip1 = K_FUNC::E_min(i+d,imc-1);
          jm1 = K_FUNC::E_max(0,j-d); jp1 = K_FUNC::E_min(j+d,jmc-1);
          km1 = K_FUNC::E_max(0,k-d); kp1 = K_FUNC::E_min(k+d,kmc-1);

          km1imjmc= km1*imjmc;
          kp1imjmc= kp1*imjmc;
          kimjmc= k*imjmc;
          
          indices[0] = im1 + jm1*imc + km1imjmc;
          indices[1] = i   + jm1*imc + km1imjmc;
          indices[2] = ip1 + jm1*imc + km1imjmc;
          
          indices[3] = im1 + j*imc + km1imjmc;
          indices[4] = i   + j*imc + km1imjmc;
          indices[5] = ip1 + j*imc + km1imjmc;

          indices[6] = im1 + jp1*imc + km1imjmc;
          indices[7] = i  +  jp1*imc + km1imjmc;
          indices[8] = ip1 + jp1*imc + km1imjmc;


          indices[9]  = im1 + jm1*imc + kimjmc;
          indices[10] = i   + jm1*imc + kimjmc;
          indices[11] = ip1 + jm1*imc + kimjmc;
          
          indices[12] = im1 + j*imc + kimjmc;
          indices[13] = ip1 + j*imc + kimjmc;
          
          indices[14] = im1 + jp1*imc + kimjmc;
          indices[15] = i  +  jp1*imc + kimjmc;
          indices[16] = ip1 + jp1*imc + kimjmc;

          indices[17] = im1 + jm1*imc + kp1imjmc;
          indices[18] = i   + jm1*imc + kp1imjmc;
          indices[19] = ip1 + jm1*imc + kp1imjmc;
          
          indices[20] = im1 + j*imc + kp1imjmc;
          indices[21] = i   + j*imc + kp1imjmc;
          indices[22] = ip1 + j*imc + kp1imjmc;

          indices[23] = im1 + jp1*imc + kp1imjmc;
          indices[24] = i  +  jp1*imc + kp1imjmc;
          indices[25] = ip1 + jp1*imc + kp1imjmc;

          for (E_Int noi = 0; noi < nindices; noi++)
          {
            ind2 = indices[noi];
            cellN[ind2] = sensor*cellN[ind2] + unmsensor*K_FUNC::E_min(cellN[ind2],0);
          }
        }
      }
    }
  }//dir = 1
}

//=============================================================================
/* Modifie le celln des centres des cellules interpolees
   a partir de la fenetre [i1,i2,j1,j2,k1,k2] de la BCOverlap */
//=============================================================================
PyObject* K_CONNECTOR::getBCOverlapInterpCellCenters(PyObject* self, PyObject* args)
{
  PyObject *coordArray, *range1;
  E_Int depth;
  if (!PYPARSETUPLEI(args,
                    "OOl", "OOi",
                    &coordArray, &range1, &depth))
  {
      return NULL;
  }
  if (PyList_Check(range1) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "getBCOverlapInterpCellCenters: 2nd argument must be a list.");
    return NULL;
  }
  if (depth != 1 && depth != 2) 
  {
     PyErr_SetString(PyExc_TypeError, 
                    "getBCOverlapInterpCellCenters: depth must be equal to 1 or 2.");
    return NULL;  
  }
  // Check: range de la BC
  FldArrayI range(6);
  E_Int sizeRange = PyList_Size(range1);
  for (int i = 0; i <  sizeRange; i++)
  {
    PyObject* tpl = PyList_GetItem(range1, i);
    if ( PyLong_Check(tpl) == 0 &&
         PyInt_Check(tpl) == 0 )
    {
      PyErr_SetString(PyExc_TypeError,
                      "getBCOverlapInterpCellCenters: range value must be an integer.");
      return NULL;
    }
    range[i] = PyLong_AsLong(tpl);
    if ( range[i] < 1 ) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "getBCOverlapInterpCellCenters: range value must be positive.");
      return NULL;
    }    
  }
  // Check: coordonnees en centres de z1
  E_Int im, jm, km;
  FldArrayF* f;
  FldArrayI* cn;
  char* varString;
  char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(coordArray, varString, f, im, jm, km, cn, eltType); 
  if ( res != 1 ) 
  {
    if ( res == 2 ) {delete f; delete cn;}
    PyErr_SetString(PyExc_TypeError, 
                    "getBCOverlapInterpCellCenters: array must be structured.");
    return NULL;   
  }
  // verification de la coherence entre le range de la paroi en centres et 
  // des dimensions de la zone
  if ( range[1] > im || range[3] > jm || range[5] > km ) 
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "getBCOverlapInterpCellCenters: wrong range max.");
    return NULL; 
  }
  E_Int i1 = range[0]; E_Int j1 = range[2]; E_Int k1 = range[4];
  E_Int i2 = range[1]; E_Int j2 = range[3]; E_Int k2 = range[5];

  E_Int posc = K_ARRAY::isCellNatureField2Present(varString);
  if ( posc == -1 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "getBCOverlapInterpCellCenters: array must contain cellN variable.");
    delete f; return NULL;
  }
  posc++;
  /* fin verifs*/
  E_Float* cellN = f->begin(posc);

  E_Int imjm = im*jm;
  E_Int inc = 0;
  switch ( depth ) 
  {
    case 1:
      for (E_Int k = k1; k <= k2; k++)
        for (E_Int j = j1; j <= j2; j++)
          for (E_Int i = i1; i <= i2; i++)
          {
            E_Int ind = (i-1) + (j-1)* im + (k-1)*imjm;
            cellN[ind] = 2.;
          }
      break;
  
    case 2:
      if ( i1 == i2 && i1 == 1 ) inc = 1;
      else if ( j1 == j2 && j1 == 1 ) inc = im;
      else if ( k1 == k2 && k1 == 1 ) inc = imjm;
      else if ( i1 == i2 && i1 > 1 ) inc = -1;
      else if ( j1 == j2 && j1  > 1 ) inc =-im;
      else if ( k1 == k2 && k2  > 1 ) inc =-imjm;
      else {printf("Fatal: getBCOverlapInterpCellCenters: invalid increment.\n"); exit(0);}  
      for (E_Int k = k1; k <= k2; k++)
        for (E_Int j = j1; j <= j2; j++)
          for (E_Int i = i1; i <= i2; i++)
          {
            E_Int ind = (i-1) + (j-1)* im + (k-1)*imjm;
            cellN[ind] = 2.; cellN[ind+inc] = 2.;
          }    
      break;
  }

  PyObject* tpl =  K_ARRAY::buildArray(*f, varString, im, jm, km);
  delete f;
  return tpl;
}
//=============================================================================
/* Determine les noeuds interpoles a partir du cellN en noeuds
   Si le celln contient des pts masques, alors les points interpoles autour 
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::getOversetHolesInterpNodes(PyObject* self, PyObject* args)
{
  PyObject *array;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLEI(args,
                    "Olls", "Oiis",
                    &array, &depth, &dir, &cellNName))
  {
      return NULL;
  }
  if (dir != 0 && dir != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: dir must be 0 or 1.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, 
                                    field, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {    
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: first argument is not recognized");
    return NULL;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);

  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: array must contain cellN variable.");
    delete field; if (res == 2) delete cn; return NULL;
  }
  posc++;

  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int npts = field->getSize();
  FldArrayI blankedCells(npts); blankedCells.setAllValuesAt(1);
  FldArrayI cellNatFld(npts); cellNatFld.setAllValuesAt(1);
  for (E_Int ind = 0; ind < npts; ind++)
  {
    if (cellNp[ind] == 2.){ blankedCells[ind] = 0; cellNatFld[ind] = 0;}
    else if (cellNp[ind] == 0.){ blankedCells[ind] = -1; cellNatFld[ind] = -1;}
  }
  if (res == 1) 
  {
    searchMaskInterpolatedCellsStruct(im, jm, km, depth, dir, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
    }
    
    PyObject* tpl =  K_ARRAY::buildArray(*field, varString, im, jm, km);
    delete field; return tpl;
  }
  else 
  {
    if ( K_STRING::cmp(eltType,"NGON")==0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getOversetHolesInterpNodes: not implemented for NGON zones.");
      delete field; delete cn; return NULL;
    }
    searchMaskInterpolatedNodesUnstr(depth, *cn, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
    }
    
    PyObject* tpl =  K_ARRAY::buildArray(*field, varString, *cn, -1, eltType);
    delete field; delete cn; return tpl;
  }
}
//=============================================================================
/* Determine les centres interpoles a partir du cellN 
   Si le celln contient des pts masques, alors les points interpoles autour 
   sont construits */
//=============================================================================
PyObject* K_CONNECTOR::getOversetHolesInterpCellCenters(PyObject* self, PyObject* args)
{
  PyObject *centersArray;
  E_Int depth; E_Int dir;
  char* cellNName;
  if (!PYPARSETUPLEI(args,
                    "Olls", "Oiis",
                    &centersArray, &depth, &dir, &cellNName))
  {
      return NULL;
  }

  if (dir != 0 && dir != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpNodes: dir must be 0 or 1.");
    return NULL;
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(centersArray, varString, 
                                    field, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {    
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpCellCenters:  first argument is not recognized");
    return NULL;
  }

  E_Int posc;
  if (strcmp(cellNName, "cellN") == 0)
    posc = K_ARRAY::isCellNatureField2Present(varString);
  else posc = K_ARRAY::isNamePresent(cellNName, varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOversetHolesInterpCellCenters: array must contain cellN variable.");
    delete field; if ( res == 2) delete cn; return NULL;
  }
  posc++;
  E_Float* cellNp = field->begin(posc);
  /* Fin des verifs */
  E_Int ncells = field->getSize();
  FldArrayI blankedCells(ncells); blankedCells.setAllValuesAt(1);
  FldArrayI cellNatFld(ncells); cellNatFld.setAllValuesAt(1);
  for (E_Int ind = 0; ind < ncells; ind++)
  {
    if (cellNp[ind] == 2.){ blankedCells[ind] = 0; cellNatFld[ind] = 0;}
    else if (cellNp[ind] == 0.){ blankedCells[ind] = -1; cellNatFld[ind] = -1;}
  }
  if (res == 1) 
  {
    searchMaskInterpolatedCellsStruct(im, jm, km, depth, dir, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
    }
    
    PyObject* tpl =  K_ARRAY::buildArray(*field, varString, im, jm, km);
    delete field; return tpl;
  }
  else 
  {
    if ( K_STRING::cmp(eltType,"NGON*")==0)
      searchMaskInterpolatedCellsNGON(depth, *cn, blankedCells, cellNatFld);
    else
      searchMaskInterpolatedCellsUnstr(eltType, depth, *cn, blankedCells, cellNatFld);
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if (cellNatFld[ind] == 0) cellNp[ind] = 2.; 
      else if (cellNatFld[ind] == -1) cellNp[ind] = 0.; 
    }
    
    PyObject* tpl =  K_ARRAY::buildArray(*field, varString, *cn, -1, eltType);
    delete field; delete cn; return tpl;
  }
}

//=============================================================================
PyObject* K_CONNECTOR::getInterpolatedPoints(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array))
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPoints: wrong arguments.");
    return NULL;
  }
  // Check: 
  E_Int im, jm, km;
  FldArrayF* f;
  FldArrayI* cn;
  char* varString;
  char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true); 
  if ( res != 1 && res != 2) 
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "getInterpolatedPoints: invalid array.");
    return NULL;   
  }
  E_Int posc = K_ARRAY::isCellNatureField2Present(varString);
  if ( posc == -1 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "getInterpolatedPoints: array must contain cellN.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posc++;

  /*fin verifs*/
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();
  char varStringOut[K_ARRAY::VARSTRINGLENGTH]; varStringOut[0] = '\0';
  E_Int nfldOut = nfld+1;
  strcpy(varStringOut,varString); strcat(varStringOut,",indcell");
  E_Float* cellnp = f->begin(posc);
  FldArrayF* fout = new FldArrayF(npts,nfldOut);
  E_Int c=0;
  for (E_Int ind=0; ind < npts; ind++)
  {
    if ( cellnp[ind] == 2.)
    { 
      for (E_Int eq = 1; eq <= nfld; eq++)
        (*fout)(c,eq) = (*f)(ind,eq);
      (*fout)(c,nfldOut) = E_Float(ind);
      c++;
    }
  }
  fout->reAllocMat(c,nfldOut);
  RELEASESHAREDB(res, array, f, cn);
  FldArrayI* cnl = new FldArrayI(0);
  PyObject* tpl = K_ARRAY::buildArray(*fout, varStringOut, *cnl, -1, 
                                      "NODE", false);
  delete fout; delete cnl;
  return tpl;
}
