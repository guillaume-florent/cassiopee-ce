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

// Routines d'identification geometrique

# include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Identifie les noeuds de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: zone a identfier
   Retourne la liste des noeuds dans la numerotation du KDT. */
// ============================================================================
PyObject* K_CONVERTER::identifyNodes(PyObject* self, PyObject* args)
{  
  PyObject* array; PyObject* hook; E_Float tol;
  if (!PYPARSETUPLEF(args, "OOd", "OOf", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3 &&
      *type != 100 && *type != 102 && *type != 103)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
  //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cnl, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: array is invalid.");
    return NULL;
  }

  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res,array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Cree le numpy de sortie
  E_Int npts = f->getSize();
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);

  PyObject* ac = K_NUMPY::buildNumpyArray(npts, 1, 1);
  E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

  // Remplissage
#pragma omp parallel default(shared)
  {
  E_Float pt[3];
  E_Float xf,yf,zf,dx,dy,dz;
  E_Int ind;
#pragma omp for schedule(dynamic)
  for (E_Int i = 0; i < npts; i++)
  {
    xf = xp[i]; yf = yp[i]; zf = zp[i];
    pt[0] = xf; pt[1] = yf; pt[2] = zf;
    ind = globalKdt->getClosest(pt); // closest pt
    dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
    //d = dx*dx + dy*dy + dz*dz;
    //if (d < tol) nptr[i] = ind+1;
    //else nptr[i] = -1;
    if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
    else nptr[i] = -1;
  }
  }
  RELEASESHAREDB(res, array, f, cnl);
  return ac;
}

// ============================================================================
/* Identifie les centres des faces de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: NGON
   Retourne la liste des faces. */
// ============================================================================
PyObject* K_CONVERTER::identifyFaces(PyObject* self, PyObject* args)
{ 
  PyObject* array; PyObject* hook;
  E_Float tol;
  if (!PYPARSETUPLEF(args, "OOd", "OOf", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3&&
      *type != 100 && *type != 102 && *type != 103)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
  //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cnl, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must be a NGON.");
    return NULL; 
  }
  if (strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must be a NGON.");
    return NULL; 
  }

  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Cree le numpy de sortie
  E_Int nfaces = (*cnl)[0];
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);
  E_Int* cnp = cnl->begin();

  PyObject* ac = K_NUMPY::buildNumpyArray(nfaces, 1, 1);
  E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

  // Remplissage
  FldArrayI posFace;
  K_CONNECT::getPosFaces(*cnl, posFace);
  E_Int* posFacep = posFace.begin();

#pragma omp parallel default(shared)
  {
#pragma omp for schedule(dynamic)
  for (E_Int i = 0; i < nfaces; i++)
  {
    E_Int pos = posFacep[i];
    E_Int* ptr = &cnp[pos];
    E_Int nv = ptr[0];
    E_Float xf = 0.,yf = 0.,zf = 0.;

#ifdef SORTHOOK
    std::vector<E_Float> sxp; sxp.reserve(1024);
    std::vector<E_Float> syp; syp.reserve(1024);
    std::vector<E_Float> szp; szp.reserve(1024);
    std::vector<E_Int> vertices; vertices.reserve(1024);

    for (E_Int n = 1; n <= nv; n++)
    {
      E_Int ind = ptr[n]-1; vertices.push_back(ind);
    }
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end() );
    for (E_Int nov = 0; nov < vertices.size(); nov++)
    {
      E_Int indv = vertices[nov];
      sxp.push_back(xp[indv]); syp.push_back(yp[indv]); szp.push_back(zp[indv]); 
    }
    std::sort(sxp.begin(), sxp.end());
    std::sort(syp.begin(), syp.end());
    std::sort(szp.begin(), szp.end());
    for (E_Int n = 0; n < sxp.size(); n++) {xf += sxp[n]; yf+= syp[n]; zf += szp[n];}
    E_Float inv = 1./E_Float(sxp.size()); xf *= inv; yf *= inv; zf *= inv;  
#else
    for (E_Int n = 1; n <= nv; n++)
    {
      E_Int ind = ptr[n]-1; 
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
    }
    E_Float inv = 1./E_Float(nv); xf *= inv; yf *= inv; zf *= inv;
#endif
    E_Float pt[3];
    pt[0] = xf; pt[1] = yf; pt[2] = zf;
    E_Int indk = globalKdt->getClosest(pt); // closest pt
    E_Float dx = xt[indk]-xf; E_Float dy = yt[indk]-yf; E_Float dz = zt[indk]-zf;
    //d = dx*dx + dy*dy + dz*dz;
    //if (d < tol) nptr[i] = ind+1; 
    //else nptr[i] = -1;
    if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = indk+1; 
    else nptr[i] = -1;
  }
  }

  RELEASESHAREDU(array, f, cnl);
  return ac;
}

// ============================================================================
/* Identifie les centres des elements de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: NGON
   Retourne les indices des points du kdtree correspondant. */
// ============================================================================
PyObject* K_CONVERTER::identifyElements(PyObject* self, PyObject* args)
{ 
  PyObject* array; PyObject* hook;
  E_Float tol;
  if (!PyArg_ParseTuple(args, "OOd", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3&&
      *type != 100 && *type != 102 && *type != 103) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
  //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cnl, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: array must be a NGON.");
    return NULL; 
  }
  
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;
  PyObject* ac = NULL;

  if (strcmp(eltType, "NGON") == 0)
  {  
    // Cree le numpy de sortie
    E_Int sizef = (*cnl)[1];
    E_Int nelts = (*cnl)[sizef+2];
    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);

    ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

    // Remplissage
    E_Int* ptr = cnl->begin();
    FldArrayI posFace; K_CONNECT::getPosFaces(*cnl, posFace);
    E_Int* posFacep = posFace.begin();
    FldArrayI posElts; K_CONNECT::getPosElts(*cnl, posElts);
    E_Int* posEltsp = posElts.begin();

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(dynamic)
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int* ptrElt = &ptr[posEltsp[i]];
      E_Int nf = ptrElt[0];
      E_Float xf = 0.; E_Float yf = 0.; E_Float zf = 0.; E_Int c = 0;
  #ifdef SORTHOOK
      std::vector<E_Float> sxp; sxp.reserve(1024);
      std::vector<E_Float> syp; syp.reserve(1024);
      std::vector<E_Float> szp; szp.reserve(1024);
      std::vector<E_Int> vertices; vertices.reserve(1024);

      for (E_Int n = 1; n <= nf; n++)
      { 
        E_Int ind = ptrElt[n]-1;
        E_Int pos = posFacep[ind];
        E_Int* ptrFace = &ptr[pos];
        E_Int nv = ptrFace[0];
      
        for (E_Int p = 1; p <= nv; p++)
        { 
          E_Int indv = ptrFace[p]-1; 
          vertices.push_back(indv);
        }
      }        

      std::sort(vertices.begin(), vertices.end());
      vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end() );

      for (E_Int nov = 0; nov < vertices.size(); nov++)
      {
        E_Int indv = vertices[nov];
        sxp.push_back(xp[indv]); syp.push_back(yp[indv]); szp.push_back(zp[indv]); 
      }
      std::sort(sxp.begin(), sxp.end());
      std::sort(syp.begin(), syp.end());
      std::sort(szp.begin(), szp.end());

      for (E_Int n = 0; n < sxp.size(); n++) { xf += sxp[n]; yf+= syp[n]; zf += szp[n];}

      E_Float inv = 1./E_Float(sxp.size()); xf *= inv; yf *= inv; zf *= inv;       

#else
      for (E_Int n = 1; n <= nf; n++)
      { 
        E_Int ind = ptrElt[n]-1;
        E_Int pos = posFacep[ind];
        E_Int* ptrFace = &ptr[pos];
        E_Int nv = ptrFace[0];
        for (E_Int p = 1; p <= nv; p++)
        {
          E_Int indp = ptrFace[p]-1; xf += xp[indp]; yf += yp[indp]; zf += zp[indp]; c++;
        }
      }
      E_Float inv = 1./E_Float(c); xf *= inv; yf *= inv; zf *= inv;
#endif
      E_Float pt[3];
      pt[0] = xf; pt[1] = yf; pt[2] = zf;
      E_Int ind = globalKdt->getClosest(pt); // closest pt
      E_Float dx = xt[ind]-xf; E_Float dy = yt[ind]-yf; E_Float dz = zt[ind]-zf;
      //d = (xt[ind]-xf)*(xt[ind]-xf)+(yt[ind]-yf)*(yt[ind]-yf)+(zt[ind]-zf)*(zt[ind]-zf);
      //if (d < tol) nptr[i] = ind+1;
      //else nptr[i] = -1;
      if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
      else nptr[i] = -1;
      //if (K_FUNC::E_abs(xf-3.)<1.e-1 && K_FUNC::E_abs(yf-9.5)<1.e-1 && K_FUNC::E_abs(zf-4.5)<1.e-1)
      //if (nptr[i] != -1)
      //printf("HERE %d: %f %f %f -> %f %f %f -> %d\n",i,xf,yf,zf,dx,dy,dz,nptr[i]);
    }
    }
  }
  else if (strcmp(eltType,"BAR")==0 || strcmp(eltType,"TRI")==0 || 
           strcmp(eltType,"QUAD")==0 || strcmp(eltType,"TETRA")==0 || 
           strcmp(eltType,"HEXA")==0 || strcmp(eltType,"PENTA")==0 ||
           strcmp(eltType,"PYRA")==0)
  {
    E_Int nelts = cnl->getSize();
    E_Int nvert = cnl->getNfld();
    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);
    ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);
    E_Float pt[3];
    E_Float inv = 1./E_Float(nvert);

#pragma omp parallel for default(shared) private(pt) schedule(dynamic)
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Float xf,yf,zf,dx,dy,dz;
      E_Int ind;
      xf = 0.; yf = 0.; zf = 0.;

#ifdef SORTHOOK
      std::vector<E_Float> sxp; sxp.reserve(1024);
      std::vector<E_Float> syp; syp.reserve(1024);
      std::vector<E_Float> szp; szp.reserve(1024);
      std::vector<E_Int> vertices; vertices.reserve(1024);

      for (E_Int n = 1; n <= nvert; n++)
      {
        ind = (*cnl)(i,n)-1; vertices.push_back(ind);
      }
      std::sort(vertices.begin(), vertices.end());
      vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end() );

      for (E_Int nov = 0; nov < vertices.size(); nov++)
      {
        E_Int indv = vertices[nov];
        sxp.push_back(xp[indv]); syp.push_back(yp[indv]); szp.push_back(zp[indv]); 
      }
      std::sort(sxp.begin(), sxp.end());
      std::sort(syp.begin(), syp.end());
      std::sort(szp.begin(), szp.end());
      for (E_Int j = 0; j < sxp.size(); j++)
      {xf+= sxp[j]; yf += syp[j]; zf += szp[j];}
      E_Float inv0 = 1./(E_Float(sxp.size()));
      xf *= inv0; yf *= inv0; zf *= inv0;
#else
      for (E_Int n = 1; n <= nvert; n++)
      {
        ind = (*cnl)(i,n)-1; 
        xf += xp[ind]; yf += yp[ind]; zf += zp[ind];        
      }
      xf*=inv; yf*=inv; zf*=inv;
#endif
      pt[0] = xf; pt[1] = yf; pt[2] = zf;

      ind = globalKdt->getClosest(pt); // closest pt
      dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
      //d = dx*dx + dy*dy + dz*dz;
      //if (d < tol) nptr[i] = ind+1;
      //else nptr[i] = -1;
      if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
      else nptr[i] = -1;
    }
  }
  else
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: invalid type of array.");
    return NULL; 
  }
  RELEASESHAREDU(array, f, cnl);
  return ac;
}
