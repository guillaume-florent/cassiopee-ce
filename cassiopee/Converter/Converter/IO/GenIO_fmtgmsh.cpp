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

// Formated gmsh (2.2) file support

# include <string.h>
# include <stdio.h>
# include "GenIO.h"
# include "Array/Array.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* gmshread
*/
//=============================================================================
E_Int K_IO::GenIO::gmshread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  E_Int res;
  E_Float t; E_Int ti; E_Int type;
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: gmshread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de l'entete
  char buf[BUFSIZE];
  res = readWord(ptrFile, buf);
  if (res == -1) { fclose(ptrFile); return 1; }
  if (strcmp(buf, "$MeshFormat") != 0) { fclose(ptrFile); return 1; }
  res = readDouble(ptrFile, t, -1);
  //printf("version %f\n", t);
  res = readInt(ptrFile, ti, -1); type = ti;
  if (type == 1) return 1; // c'est un binary file
  //printf("file type %d\n", ti);
  res = readInt(ptrFile, ti, -1);
  //printf("data size %d\n", ti);
  res = readGivenKeyword(ptrFile, "$ENDMESHFORMAT");
  
  // Lecture Vertices (Global)
  res = readGivenKeyword(ptrFile, "$NODES");
  res = readInt(ptrFile, ti, -1);
  E_Int nb = E_Int(ti);
  //printf("Number of nodes %d\n", nb);
  FldArrayF f(nb, 3);
  E_Float* f1 = f.begin(1); E_Float* f2 = f.begin(2); E_Float* f3 = f.begin(3);
  FldArrayI indirNodes(nb);
  for (E_Int i = 0; i < nb; i++)
  {
    res = readInt(ptrFile, ti, -1); indirNodes[i] = ti;
    res = readDouble(ptrFile, t, -1); f1[i] = t;
    res = readDouble(ptrFile, t, -1); f2[i] = t;
    res = readDouble(ptrFile, t, -1); f3[i] = t;
    //printf("%f %f %f\n", f(i,1), f(i,2), f(i,3));
  }
  res = readGivenKeyword(ptrFile, "$ENDNODES");

  /* Elements by zone type */
  res = readGivenKeyword(ptrFile, "$ELEMENTS");
  res = readInt(ptrFile, ti, -1);
  E_Int ne = E_Int(ti); // Global
  //printf("Number of elements %d\n", ne);
  FldArrayI indirElements(ne);
  E_Int nNODE = 0; FldArrayI* indNODE = NULL; // type = 15
  E_Int nBAR = 0; FldArrayI* cnBAR = NULL; // type = 1
  E_Int nTRI = 0; FldArrayI* cnTRI = NULL; // type = 2
  E_Int nQUAD = 0; FldArrayI* cnQUAD = NULL; // type = 3
  E_Int nTETRA = 0; FldArrayI* cnTETRA = NULL; // type = 4
  E_Int nHEXA = 0; FldArrayI* cnHEXA = NULL; // type = 5
  E_Int nPENTA = 0; FldArrayI* cnPENTA = NULL; // type = 6
  E_Int nPYRA = 0; FldArrayI* cnPYRA = NULL; // type = 7
  E_Int nDiscard = 0;

  E_Int tagl, ind;
  /* Compte les elements par type */
  E_LONG pos = KFTELL(ptrFile);
#define READI readInt(ptrFile, ti, -1)
#include "GenIO_gmsh1.h"
  //printf("Elements BAR=%d TRI=%d QUAD=%d TETRA=%d HEXA=%d NODES=%d\n", 
  //       nBAR, nTRI, nQUAD, nTETRA, nHEXA, nNODE);

  /* Allocations */
  if (nBAR > 0) {
    cnBAR = new FldArrayI(nBAR, 2);
    connect.push_back(cnBAR);
    eltType.push_back(1); nBAR = 0; }
  if (nTRI > 0) {
    cnTRI = new FldArrayI(nTRI, 3);
    connect.push_back(cnTRI);
    eltType.push_back(2); nTRI = 0; }
  if (nQUAD > 0) {
    cnQUAD = new FldArrayI(nQUAD, 4);
    connect.push_back(cnQUAD);
    eltType.push_back(3); nQUAD = 0; }
  if (nTETRA > 0) {
    cnTETRA = new FldArrayI(nTETRA, 4);
    connect.push_back(cnTETRA);
    eltType.push_back(4); nTETRA = 0; }
  if (nHEXA > 0) {
    cnHEXA = new FldArrayI(nHEXA, 8);
    connect.push_back(cnHEXA);
    eltType.push_back(7); nHEXA = 0; }
  if (nPENTA > 0) {
    cnPENTA = new FldArrayI(nPENTA, 6);
    connect.push_back(cnPENTA);
    eltType.push_back(6); nPENTA = 0; }
  if (nPYRA > 0) {
    cnPYRA = new FldArrayI(nPYRA, 5);
    connect.push_back(cnPYRA);
    eltType.push_back(5); nPYRA = 0; }
  if (nNODE > 0) {
    indNODE = new FldArrayI(nNODE);
    connect.push_back(new FldArrayI());
    eltType.push_back(0);
    nNODE = 0; }

  /* Lecture reelle des elements par type */
  KFSEEK(ptrFile, pos, SEEK_SET);
#include "GenIO_gmsh2.h"

  /* Nodes duplications */
  if (nBAR > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nTRI > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nQUAD > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nTETRA > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nHEXA > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nPENTA > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nPYRA > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  //printf("node %d\n", nNODE);
  if (nNODE > 0)
  {
    FldArrayF* an = new FldArrayF(nNODE,3);
    for (E_Int i = 0; i < nNODE; i++) 
    { ind = (*indNODE)[i]-1; ind = indirNodes[ind]-1;
      (*an)(i,1) = f1[ind]; (*an)(i,2) = f2[ind]; (*an)(i,3) = f3[ind]; }
    unstructField.push_back(an);
    delete indNODE;
  }

  // Cree le nom des zones
  //printf("Number of zones %d\n", unstructField.size());
  for (unsigned int i=0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%d", i);
    zoneNames.push_back(zoneName);
  }
  //printf("sizes: %d %d\n", unstructField.size(), connect.size());

  varString = new char [8];
  strcpy(varString, "x,y,z");
  fclose(ptrFile);

  return 0;
}

//=============================================================================
// Only write NODE, BAR, TRI, QUAD, TETRA, HEXA, PYRA, PENTA.
//=============================================================================
E_Int K_IO::GenIO::gmshwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  E_Int nzone = unstructField.size();
  E_Int nvalidZones = 0;
  for (E_Int zone = 0; zone < nzone; zone++)
  {
    // NODE, BAR, TRI, QUADS, TETRA, HEXA , PENTA, PYRA, supported
    if (eltType[zone] < 8) nvalidZones++;
    else
      printf("Warning: gmshwrite: zone %d not written (not a valid elements in zone).", zone);
  }

  if (nvalidZones == 0) return 1;

  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1)
  {
    printf("Warning: gmshwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;
    
  char format1[40]; char format2[40]; char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format for data
  strcpy(format1, "%d ");
  sprintf(format2,"%s%s%s", dataFmt, dataFmt, dataFmtl);
  strcat(format1, format2);
  strcat(format1, "\n");

  // Concatenate all vertices in one field
  E_Int size = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    size += unstructField[i]->getSize();
  }
  FldArrayF* vertices = new FldArrayF(size, 3);
  FldArrayF& v = *vertices;
  E_Int c = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    FldArrayF& field = *unstructField[i];
    for (E_Int n = 0; n < field.getSize(); n++)
    {
      v(c, 1) = field(n, posx);
      v(c, 2) = field(n, posy);
      v(c, 3) = (posz>0) ? field(n, posz) : 0.; c++;
    }
  }

  // Ecriture
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL)
  {
    printf("Warning: gmshwrite: I can't open file %s.\n", file);
    return 1;
  }

  fprintf(ptrFile, "$MeshFormat\n");
  fprintf(ptrFile, "2.2 0 8\n");
  fprintf(ptrFile, "$EndMeshFormat\n");

  fprintf(ptrFile, "$Nodes\n");
  fprintf(ptrFile, "%d\n", v.getSize());
  for (E_Int i = 0; i < v.getSize(); i++)
    fprintf(ptrFile, format1, i+1, v(i,1), v(i,2), v(i,3));
  fprintf(ptrFile, "$EndNodes\n");

  fprintf(ptrFile, "$Elements\n");
  size = 0;
  for (E_Int i = 0; i < nzone; i++)
  {
    size += connect[i]->getSize();
  }
  fprintf(ptrFile, "%d\n", size);

  // Connectivite par elts par rapport a la definition globale des vertices
  E_Int shift = 0;
  E_Int ind = 1;
  for (E_Int i = 0; i < nzone; i++)
  {
    FldArrayI& cn = *connect[i];
    E_Int elt = eltType[i];
    switch(elt)
    {
      case 0: // NODE
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, "%d 15 2 0 14 %d\n", ind, ind);
          ind++;
        }
      }
      break;
      case 1: // BAR
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, "%d 1 2 0 14 %d %d\n", ind, cn(n,1)+shift, cn(n,2)+shift);
          ind++;
        }
      }
      break;
      case 2: // TRI
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, "%d 2 2 0 14 %d %d %d\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift);
          ind++;
        }
      }
      break;
      case 3: // QUAD
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, "%d 3 2 0 14 %d %d %d %d\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,cn(n,4)+shift);
          ind++;
        }
      }
      break;
      case 4: // TETRA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, "%d 4 2 0 14 %d %d %d %d\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,cn(n,4)+shift);
          ind++;
        }
      }
      break;
      case 5: // PYRA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, "%d 7 2 0 14 %d %d %d %d %d\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,
                  cn(n,4)+shift,cn(n,5)+shift);
          ind++;
        }
      }
      break;
      case 6: // PENTA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, "%d 6 2 0 14 %d %d %d %d %d %d\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,
                  cn(n,4)+shift,cn(n,5)+shift,cn(n,6)+shift);
          ind++;
        }
      }
      break;
      case 7: // HEXA
      {
        for (E_Int n = 0; n < cn.getSize(); n++)
        {
          fprintf(ptrFile, "%d 5 2 0 14 %d %d %d %d %d %d %d %d\n", ind, 
                  cn(n,1)+shift, cn(n,2)+shift,cn(n,3)+shift,
                  cn(n,4)+shift,cn(n,5)+shift,cn(n,6)+shift,
                  cn(n,7)+shift,cn(n,8)+shift);
          ind++;
        }
      }
      break;
    }
    shift += unstructField[i]->getSize();
  }
  fprintf(ptrFile, "$EndElements\n");

  fclose(ptrFile);
  return 0;
}
