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

# include "kcore.h"
# include "Data.h"
# include "cplot.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Cree une zone structuree pour CPlot a partir d'un FldArray 
   IN: structF: coord + fields
   IN: varString: chaine des variables
   IN: posx, posy, posz: position des coords dans structF
   IN: ni,nj,nk: nbre de noeuds
   IN: zoneName: nom de la zone
   IN: zoneTags: tags de la zone (peut etre NULL)
   IN: referenceZone: sert de reference pour les variables
*/
//=============================================================================
StructZone* Data::createStructZone(FldArrayF* structF, char* varString,
                                   E_Int posx, E_Int posy, E_Int posz,
                                   E_Int ni, E_Int nj, E_Int nk,
                                   char* zoneName, char* zoneTags,
                                   Zone* referenceZone)
{
  StructZone* sz = new StructZone( ptrState, createZoneImpl() );
  StructZone& z = *sz;
  strcpy(z.zoneName, zoneName);
 
  z.ni = ni;
  z.nj = nj;
  z.nk = nk;
  z.npts = z.ni * z.nj * z.nk;
  if (referenceZone != NULL) z.nfield = referenceZone->nfield;
  else z.nfield = structF->getNfld()-3;
  if (z.nk != 1) z.dim = 3;
  else if (z.nj != 1) z.dim = 2;
  else z.dim = 1;
  z.x = new E_Float[z.npts];
  memcpy(z.x, structF->begin(posx), z.npts*sizeof(E_Float));
  z.y = new E_Float[z.npts];
  memcpy(z.y, structF->begin(posy), z.npts*sizeof(E_Float));
  z.z = new E_Float[z.npts];
  memcpy(z.z, structF->begin(posz), z.npts*sizeof(E_Float));

  vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int varsSize = vars.size();

  // Allocation of var fields
  reallocNFieldArrays(z.nfield);
  z.f = new double* [z.nfield];
  for (E_Int i = 0; i < z.nfield; i++) z.f[i] = NULL;
  z.varnames = new char* [z.nfield];
  for (E_Int n = 0; n < z.nfield; n++) 
    z.varnames[n] = new char [MAXSTRINGLENGTH];
  z.minf = new double [z.nfield];
  z.maxf = new double [z.nfield];

  if (referenceZone != NULL)
  {
    for (E_Int n = 0; n < referenceZone->nfield; n++)
    {
      for (E_Int p = 0; p < varsSize; p++)
      {
        if (K_STRING::cmp(vars[p], referenceZone->varnames[n]) == 0)
        {
          z.f[n] = new E_Float[z.npts];
          memcpy(z.f[n], structF->begin(p+1), z.npts*sizeof(E_Float));
          strcpy(z.varnames[n], vars[p]); break;
        }
      }
    }
    for (E_Int n = 0; n < referenceZone->nfield; n++)
    {
      if (z.f[n] == NULL)
      {
        z.f[n] = new E_Float[z.npts];
        if (K_STRING::cmp(referenceZone->varnames[n], "cellN") == 0)
        { for (int i = 0; i < z.npts; i++) z.f[n][i] = 1.; }
        else { for (int i = 0; i < z.npts; i++) z.f[n][i] = 0.; }
        strcpy(z.varnames[n], referenceZone->varnames[n]);
      }
    }
  }
  else // no reference zone
  {
    E_Int p = 0;
    for (E_Int n = 0; n < structF->getNfld(); n++)
    {
      if (n != posx-1 && n != posy-1 && n != posz-1)
      {
        z.f[p] = new E_Float[z.npts];
        memcpy(z.f[p], structF->begin(n+1), z.npts*sizeof(E_Float));
        strcpy(z.varnames[p], vars[n]); p++;
      }
    }
  }

  for (E_Int n = 0; n < varsSize; n++) delete [] vars[n];

  if (zoneTags != NULL)
  {
    strcpy(z.renderTag, zoneTags);
    codeFromRenderTag(z, z.renderTag, z.colorR, z.colorG, z.colorB, 
                      z.material, z.blending, z.meshOverlay, 
                      z.shaderParam1, z.shaderParam2);
  }
  else
  {
    strcpy(z.renderTag, "None:None:None:None");
    z.colorR = -1; z.colorG = -1; z.colorB = -1; z.material = -1;
    z.blending = -1.; z.meshOverlay = 0; z.shaderParam1 = 1.;
    z.shaderParam2 = 1.;
  }

  z.surf = NULL; z.compNorm();

  z.activePlane = 0;
  z.iPlane = -1;
  z.jPlane = -1;
  z.kPlane = -1;
  z.iLine = 0;
  z.jLine = 0;
  z.kLine = 0;
  z.blank = 0;
  z.active = 1;
  z.selected = 0;
  findMinMax(&z); findFMinMax(&z);

  return sz;
}

//=============================================================================
/* Cree une zone non-structuree pour CPlot a partir d'un FldArray 
   IN: unstrF: coord + fields
   IN: varString: chaine des variables
   IN: posx, posy, posz: position des coords dans unstrF
   IN: cn: connectivite
   IN: eltType: type d'element
   IN: zoneName: nom de la zone
   IN: zoneTags: tags de la zone (peut etre NULL)
*/
//=============================================================================
UnstructZone* Data::createUnstrZone(FldArrayF* unstrF, char* varString,
                                    E_Int posx, E_Int posy, E_Int posz,
                                    FldArrayI* cn, char* eltType,
                                    char* zoneName, char* zoneTags,
                                    Zone* referenceZone)
{
  UnstructZone* uz = new UnstructZone( ptrState, createZoneImpl() );
  UnstructZone& z = *uz;
  strcpy(z.zoneName, zoneName);
  z.npts = unstrF->getSize();
  z.np = z.npts;
  if (referenceZone != NULL) z.nfield = referenceZone->nfield;
  else z.nfield = unstrF->getNfld()-3;
  
  z.x = new E_Float[z.npts];
  memcpy(z.x, unstrF->begin(posx), z.npts*sizeof(E_Float));
  z.y = new E_Float[z.npts];
  memcpy(z.y, unstrF->begin(posy), z.npts*sizeof(E_Float));
  z.z = new E_Float[z.npts];
  memcpy(z.z, unstrF->begin(posz), z.npts*sizeof(E_Float));
  z.ne = cn->getSize();
  vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int varsSize = vars.size();

  // Allocation of var fields
  reallocNFieldArrays(z.nfield);
  z.f = new double* [z.nfield];
  for (E_Int i = 0; i < z.nfield; i++) z.f[i] = NULL;
  z.varnames = new char* [z.nfield];
  for (E_Int n = 0; n < z.nfield; n++) 
    z.varnames[n] = new char [MAXSTRINGLENGTH];
  z.minf = new double [z.nfield];
  z.maxf = new double [z.nfield];

  if (referenceZone != NULL)
  {
    for (E_Int n = 0; n < referenceZone->nfield; n++)
    {
      for (E_Int p = 0; p < varsSize; p++)
      {
        if (K_STRING::cmp(vars[p], referenceZone->varnames[n]) == 0)
        {
          z.f[n] = new E_Float[z.npts];
          memcpy(z.f[n], unstrF->begin(p+1), z.npts*sizeof(E_Float));
          strcpy(z.varnames[n], vars[p]); break;
        }
      }
    }
    for (E_Int n = 0; n < referenceZone->nfield; n++)
    {
      if (z.f[n] == NULL)
      {
        z.f[n] = new E_Float[z.npts];
        if (K_STRING::cmp(referenceZone->varnames[n], "cellN") == 0)
        { for (int i = 0; i < z.npts; i++) z.f[n][i] = 1.; }
        else { for (int i = 0; i < z.npts; i++) z.f[n][i] = 0.; }
        strcpy(z.varnames[n], referenceZone->varnames[n]);
      }
    }
  }
  else // no reference zone
  {
    E_Int p = 0;
    for (E_Int n = 0; n < unstrF->getNfld(); n++)
    {
      if (n != posx-1 && n != posy-1 && n != posz-1)      
      {
        z.f[p] = new E_Float[z.npts];
        memcpy(z.f[p], unstrF->begin(n+1), z.npts*sizeof(E_Float));
        strcpy(z.varnames[p], vars[n]); p++;
      }
    }
  }

  for (E_Int n = 0; n < varsSize; n++) delete [] vars[n];
  
  if (zoneTags != NULL)
  {
    strcpy(z.renderTag, zoneTags);
    codeFromRenderTag(z, z.renderTag, z.colorR, z.colorG, z.colorB, 
                      z.material, z.blending, z.meshOverlay, 
                      z.shaderParam1, z.shaderParam2);
  }
  else
  {
    strcpy(z.renderTag, "None:None:None:None");
    z.colorR = -1.; z.colorG = -1.; z.colorB = -1.; z.material = -1;
    z.blending = -1.; z.meshOverlay = 0; z.shaderParam1 = 1.;
    z.shaderParam2 = 1.;
  }

  if (K_STRING::cmp(eltType, "NODE") == 0)
  {
    z.eltType = 0;
    z.eltSize = 0;
    z.dim = 0;
  }
  else if (K_STRING::cmp(eltType, "BAR") == 0)
  {
    z.eltType = 1;
    z.eltSize = 2;
    z.dim = 1;
  }
  else if (K_STRING::cmp(eltType, "TRI") == 0)
  {
    z.eltType = 2;
    z.eltSize = 3;
    z.dim = 2;
  }
  else if (K_STRING::cmp(eltType, "QUAD") == 0)
  {
    z.eltType = 3;
    z.eltSize = 4;
    z.dim = 2;
  }
  else if (K_STRING::cmp(eltType, "TETRA") == 0)
  {
    z.eltType = 4;
    z.eltSize = 4;
    z.dim = 3;
  }
  else if (K_STRING::cmp(eltType, "PENTA") == 0)
  {
    z.eltType = 5;
    z.eltSize = 6;
    z.dim = 3;
  }
  else if (K_STRING::cmp(eltType, "PYRA") == 0)
  {
    z.eltType = 6;
    z.eltSize = 5;
    z.dim = 3;
  }
  else if (K_STRING::cmp(eltType, "HEXA") == 0)
  {
    z.eltType = 7;
    z.eltSize = 8;
    z.dim = 3;
  }
  else if (K_STRING::cmp(eltType, "NGON") == 0)
  {
    z.ne = (*cn)[2+(*cn)[1]];
    z.eltType = 10;
    z.eltSize = 1;
    z.dim = 3;
  }
  else
  {
    printf("Warning: element type is unknown. Set to TRI.\n");
    z.eltType = 2;
    z.eltSize = 3;
    z.dim = 2;
  }
  E_Int size = cn->getSize() * cn->getNfld();
  z.connect = new E_Int[size];
  memcpy(z.connect, cn->begin(), size*sizeof(E_Int));
 
  z.posFaces = NULL;

  if (z.eltType == 10) // NGONS
  {
    // calcul posFaces (position des faces dans connect)
    int nfaces = NFACES(z.connect);
    z.posFaces = new E_Int[nfaces];
    int c = POSFACES(z.connect); int l;
    for (int i = 0; i < nfaces; i++)
    {
      z.posFaces[i] = c; l = z.connect[c]; c += l+1;
    }
    
    // calcul le nombre d'elements 1D et 2D
    int nelts = NELTS(z.connect); c = POSELTS(z.connect);
    int dim, s, c1, c2;
    z.nelts1D = 0; z.nelts2D = 0;
    for (int i = 0; i < nelts; i++)
    {
      l = z.connect[c]; // nbre de faces
      dim = 0;
      for (int j = 0; j < l; j++)
      {
        s = z.posFaces[z.connect[c+j+1]-1];
        dim = max(dim, z.connect[s]);
      }
      if (dim == 1) { z.nelts1D++; }
      else if (dim == 2) { z.nelts2D++; }
      c += l+1;
    }
    
    //printf("1D: %d, 2D: %d\n", z.nelts1D, z.nelts2D);
    if (z.nelts1D > 0) z.posElts1D = new E_Int[z.nelts1D];
    if (z.nelts2D > 0) z.posElts2D = new E_Int[z.nelts2D];
    c = POSELTS(z.connect);
    c1 = 0; c2 = 0;
    for (int i = 0; i < nelts; i++)
    {
      l = z.connect[c]; // nbre de faces
      dim = 0;
      for (int j = 0; j < l; j++)
      {
        s = z.posFaces[z.connect[c+j+1]-1];
        dim = max(dim, z.connect[s]);
      }
      if (dim == 1) { z.posElts1D[c1] = c; c1++; }
      else if (dim == 2) { z.posElts2D[c2] = c; c2++; }
      c += l+1;
    }
  }

  z.surf = NULL; z.compNorm();
  z.blank = 0;
  z.active = 1;
  z.selected = 0;
  findMinMax(&z); findFMinMax(&z);

  return uz;
}
