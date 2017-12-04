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
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Effectue les transferts par interpolation */
//=============================================================================
PyObject* K_CONNECTOR::setInterpTransfers(PyObject* self, PyObject* args)
{
  PyObject *arrayR, *arrayD;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyVariables; 
  E_Float AngleX, AngleY, AngleZ; 
  if (!PYPARSETUPLEF(args, "OOOOOOO(ddd)", "OOOOOOO(fff)",
                     &arrayR, &arrayD,  &pyVariables, &pyIndRcv, &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, 
                     &AngleX, &AngleY, &AngleZ))
  {
      return NULL;
  }
  
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
                    "setInterpTransfers: 1st arg is not a valid array.");
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
  if (resd != 2 && resd != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpTransfers: 2nd arg is not a valid array.");
    RELEASESHAREDB(resr, arrayR, fr, cnr); 
    return NULL; 
  }

  E_Int meshtype = resd; // 1: structure, 2: non structure

  if (resd == 2)
  {
    if (K_STRING::cmp(eltTypeD, "TETRA") != 0 && 
        K_STRING::cmp(eltTypeD, "NGON") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "setInterpTransfers: unstructured donor zone must only be TETRA or NGON.");
      RELEASESHAREDB(resr, arrayR, fr, cnr);
      RELEASESHAREDB(resd, arrayD, fd, cnd); 
      return NULL; 
    }
  }
  /*---------------------------------------------------*/
  /*  Extrait les positions des variables a transferer */
  /*---------------------------------------------------*/
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  E_Int posvr, posvd;
  if (PyList_Check(pyVariables) != 0)
  {
    int nvariables = PyList_Size(pyVariables);
    if ( nvariables > 0 )
    {
      for (int i = 0; i < nvariables; i++)
      {
        PyObject* tpl0 = PyList_GetItem(pyVariables, i);
        if (PyString_Check(tpl0) == 0)
          PyErr_Warn(PyExc_Warning, "setInterpTransfers: variable must be a string. Skipped.");
        else 
        {
          char* varname = PyString_AsString(tpl0);        
          posvd = K_ARRAY::isNamePresent(varname, varStringD);      
          posvr = K_ARRAY::isNamePresent(varname, varStringR);      
          if (posvd != -1 && posvr != -1) 
          {
            posvarsD.push_back(posvd+1);
            posvarsR.push_back(posvr+1);
          }
        }
      }
    }
    else // toutes les variables communes sont transferees sauf le cellN
    {
      E_Int poscd = K_ARRAY::isCellNatureField2Present(varStringD)+1;
      E_Int poscr = K_ARRAY::isCellNatureField2Present(varStringR)+1;
      char* varStringC; // chaine de caractere commune
      E_Int l = strlen(varStringR);
      varStringC = new char [l+1];
      K_ARRAY::getPosition(varStringR, varStringD, 
                           posvarsR, posvarsD, varStringC);
      
      delete [] varStringC; 
      if (poscd != 0) 
        posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poscd), posvarsD.end());
      if (poscr != 0)
        posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poscr), posvarsR.end());    
    }
  }
  else 
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr); 
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpTransfers: name of transfered variables must be defined by a list.");
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
    PyErr_SetString(PyExc_TypeError,"setInterpTransfers: 4th to 6th arg must be a numpy of integers. 7th arg a numpy floats ");
    return NULL;
  }

  // Extraction de l angle de rotation
  E_Int dirR = 0; E_Float theta = 0.;
  if ( K_FUNC::E_abs(AngleX) > 0.) {dirR = 1; theta=AngleX;}
  else if (K_FUNC::E_abs(AngleY) > 0.){dirR=2; theta=AngleY;}
  else if (K_FUNC::E_abs(AngleZ) > 0.) {dirR=3; theta=AngleZ;} 
  E_Int posvx=-1, posvy=-1, posvz=-1, posmx=-1, posmy=-1, posmz=-1;
  if ( dirR > 0 )
  {
    posvx = K_ARRAY::isNamePresent("VelocityX", varStringR);
    posvy = K_ARRAY::isNamePresent("VelocityY", varStringR);
    posvz = K_ARRAY::isNamePresent("VelocityZ", varStringR);
    posmx = K_ARRAY::isNamePresent("MomentumX", varStringR);
    posmy = K_ARRAY::isNamePresent("MomentumY", varStringR);
    posmz = K_ARRAY::isNamePresent("MomentumZ", varStringR);
  }
  //
  E_Int nvars   = fr->getNfld();//nb de champs a interpoler
  E_Int* ptrcnd = cnd->begin(); 
  E_Int cnNfldD = cnd->getNfld();

  PyObject* tpl;
  if (resr == 1) 
  {
    tpl = K_ARRAY::buildArray(nvars, varStringR, imr, jmr, kmr);
  }
  else // unstructured 
  {
    E_Int crsize = cnr->getSize()*cnr->getNfld(); 
    tpl = K_ARRAY::buildArray(nvars, varStringR,
                              fr->getSize(), cnr->getSize(),
                              -1, eltTypeR, false, crsize);
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cnr->begin(), cnr->getSize()*cnr->getNfld());
  }

  E_Float* frp  = K_ARRAY::getFieldPtr(tpl);
  //FldArrayF fieldROut(nbRcvPts, nvars, frp, true);
  FldArrayF fieldROut( fr->getSize() , nvars, frp, true);
  fieldROut = *fr;

  // Transferts 
  // Types valides: 0,1, 2, 3, 4, 5 

  nvars = posvarsR.size();
  
  vector<E_Float*> vectOfRcvFields(nvars);
  vector<E_Float*> vectOfDnrFields(nvars);

  for (E_Int eq = 0; eq < nvars; eq++)
   {
    vectOfRcvFields[eq] = fieldROut.begin(posvarsR[eq]);
    vectOfDnrFields[eq] = fd->begin(posvarsD[eq]);
   }

////
////
//  Interpolation parallele
////  
////  
# include "commonInterpTransfers_indirect.h" 

  // Prise en compte de la periodicite par rotation
  if ( dirR != 0 )
  {
    # include "includeTransfers.h"
  }

  // sortie
  RELEASESHAREDB(resr, arrayR, fr, cnr); 
  RELEASESHAREDB(resd, arrayD, fd, cnd); 
  RELEASESHAREDN(pyIndRcv    , rcvPtsI    );
  RELEASESHAREDN(pyIndDonor  , donorPtsI  );
  RELEASESHAREDN(pyArrayTypes, typesI     );
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  return tpl;
}

//=============================================================================
// Idem: in place + from zone
//=============================================================================
PyObject* K_CONNECTOR::_setInterpTransfers(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;

  PyObject *zoneR, *zoneD;
  PyObject *pyIndRcv, *pyIndDonor; PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs; PyObject *pyVariables;
  E_Int loc, vartype, compact;
  char* cellNVariable;
  E_Float AngleX, AngleY, AngleZ; 

  if (!PYPARSETUPLE(args,
                    "OOOOOOOlllssssddd", "OOOOOOOiiissssddd",
                    "OOOOOOOlllssssfff", "OOOOOOOiiissssfff",
                    &zoneR, &zoneD, &pyVariables, &pyIndRcv, 
                    &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, &loc , &vartype, &compact, &cellNVariable,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters, 
                    &AngleX, &AngleY, &AngleZ))
  {
      return NULL;
  }
  // Extraction de l angle de rotation
  E_Int dirR = 0; E_Float theta = 0.;
  if ( K_FUNC::E_abs(AngleX) > 0.) {dirR = 1; theta=AngleX;}
  else if (K_FUNC::E_abs(AngleY) > 0.){dirR=2; theta=AngleY;}
  else if (K_FUNC::E_abs(AngleZ) > 0.) {dirR=3; theta=AngleZ;} 

  vector<PyArrayObject*> hook;
  E_Int imdjmd, imd,jmd,kmd, cnNfldD, nvars,ndimdxR, ndimdxD,meshtype;
  E_Float* iptroD; E_Float* iptroR; 

# include "extract_interpD.h"

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI, true);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts = rcvPtsI->getSize();

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;

  //codage general (lent ;-) )
  if (compact==0)
  {
    // recupere les champs du donneur (nodes)
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
    // recupere les champs du receveur (centers)
    E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
    char* varStringR;  vector<E_Int> locsR;
    vector<E_Int*> cnr;
    K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                          fieldsR, locsR, imr, jmr, kmr,
                          cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                          GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
    if (varStringD == NULL)
    {
      RELEASESHAREDZ(hook, varStringD, eltTypeD);
      RELEASESHAREDN(pyIndRcv, rcvPtsI);
      RELEASESHAREDN(pyIndDonor, donorPtsI);
      RELEASESHAREDN(pyArrayTypes, typesI);
      RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
      PyErr_SetString(PyExc_TypeError, 
                      "_setInterpTransfers: no field found in donor zones.");
      return NULL;
    }
    if (varStringR == NULL) 
    {
      RELEASESHAREDZ(hook, varStringR, eltTypeR);
      RELEASESHAREDN(pyIndRcv, rcvPtsI);
      RELEASESHAREDN(pyIndDonor, donorPtsI);
      RELEASESHAREDN(pyArrayTypes, typesI);
      RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
      PyErr_SetString(PyExc_TypeError, 
                      "_setInterpTransfers: no field found in receiver zones.");
      return NULL;
    }

    // Extrait les positions des variables a transferer
    E_Int posvr, posvd;
    E_Int initAll = false;
    E_Int poscd = K_ARRAY::isNamePresent(cellNVariable, varStringD);      
    E_Int poscr = K_ARRAY::isNamePresent(cellNVariable, varStringR);

    if (PyList_Check(pyVariables) != 0)
    {
      int nvariables = PyList_Size(pyVariables);
      if (nvariables > 0)
      {
        for (int i = 0; i < nvariables; i++)
        {
          PyObject* tpl0 = PyList_GetItem(pyVariables, i);
          if (PyString_Check(tpl0) == 0)
            PyErr_Warn(PyExc_Warning, "_setInterpTransfers: variable must be a string. Skipped.");
          else 
          {
            char* varname = PyString_AsString(tpl0);        
            posvd = K_ARRAY::isNamePresent(varname, varStringD);      
            posvr = K_ARRAY::isNamePresent(varname, varStringR);      
            if (posvd != -1 && posvr != -1) 
            {
              posvarsD.push_back(posvd);
              posvarsR.push_back(posvr);
            }
          }
        }
      }
      else initAll = true;
    } // ckeck list of variables non empty
    else { initAll = true; }

    E_Int posvx=-1, posvy=-1, posvz=-1, posmx=-1, posmy=-1, posmz=-1;
    if ( dirR > 0 )
    {
      posvx = K_ARRAY::isNamePresent("VelocityX", varStringR);
      posvy = K_ARRAY::isNamePresent("VelocityY", varStringR);
      posvz = K_ARRAY::isNamePresent("VelocityZ", varStringR);
      posmx = K_ARRAY::isNamePresent("MomentumX", varStringR);
      posmy = K_ARRAY::isNamePresent("MomentumY", varStringR);
      posmz = K_ARRAY::isNamePresent("MomentumZ", varStringR);
    }
    if (initAll == true)// all common variables are transfered
    {
      posvd = K_ARRAY::isCellNatureField2Present(varStringD);
      posvr = K_ARRAY::isCellNatureField2Present(varStringR);
      char* varStringC; // chaine de caractere commune
      E_Int l = strlen(varStringR);
      varStringC = new char [l+1];
      // les positions demarrent a 1 
      K_ARRAY::getPosition(varStringR, varStringD, 
                           posvarsR, posvarsD, varStringC);
      delete [] varStringC;
      E_Int sizeVarsD = posvarsD.size();
      E_Int sizeVarsR = posvarsR.size();

      for (E_Int i = 0; i < sizeVarsD; i++) posvarsD[i] -= 1;
      for (E_Int i = 0; i < sizeVarsR; i++) posvarsR[i] -= 1;
    
      if (posvd != -1) 
        posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), posvd), posvarsD.end());
      if (posvr != -1)
        posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), posvr), posvarsR.end());
    }
    if (poscd > -1 && poscr > -1) // cellNVariable exists : do not interpolate but specific update
    {
      posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poscd), posvarsD.end());
      posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poscr), posvarsR.end());
    }
    delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

    // -- no check (perfo) --
    // Transferts 
    // Types valides: 0, 1, 2, 3, 4, 5 
    nvars = posvarsR.size();

    vector<E_Float*> vectOfRcvFields(nvars);
    vector<E_Float*> vectOfDnrFields(nvars);

    for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq] = fieldsR[posvarsR[eq]];
      vectOfDnrFields[eq] = fieldsD[posvarsD[eq]];
    }
    // interpolation of all fields
# include "commonInterpTransfers_indirect.h" 

    // transfer of cellN variable
    if ( poscd > -1 && poscr > -1) // cellNVariable exists
    {  
# include "commonCellNTransfersStrict.h"    
    }    
    // Prise en compte de la periodicite par rotation
    if ( dirR != 0 )
    {
    # include "includeTransfers.h"
    }
  }// end of  compact = 0
  
  else  // compacted fields
  {
    // les variables a transferer sont compactees: on recupere uniquement la premiere et la taille 
#include "getfromzonecompact.h"
    if( vartype <= 3 &&  vartype >= 1) nvars =5;
    else                               nvars =6;

    vector<E_Float*> vectOfRcvFields(nvars);
    vector<E_Float*> vectOfDnrFields(nvars);
    
    for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq] = iptroR + eq*ndimdxR;
      vectOfDnrFields[eq] = iptroD + eq*ndimdxD;
    }
    //Parallel interpolation of all fields
# include "commonInterpTransfers_indirect.h" 

    // Prise en compte de la periodicite par rotation
    if ( dirR != 0 )
    {
      printf("Warning: _setInterpTransfers: no correction for periodicity by rotation can be applied in compacted version.\n");
     //# include "includeTransfers.h"
    }
  }// end of compacted field transfers

  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
  RELEASESHAREDN(pyIndRcv    , rcvPtsI    );
  RELEASESHAREDN(pyIndDonor  , donorPtsI  );
  RELEASESHAREDN(pyArrayTypes, typesI     );
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
// Idem: in place + from zone + tc compact au niveau zone donneuse
//=============================================================================
PyObject* K_CONNECTOR::__setInterpTransfers(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zoneD;
  PyObject *pyVariables;
  PyObject *pyParam_int, *pyParam_real;
  E_Int loc, vartype, compact, flagibc, bctype;
  E_Float gamma, cv, muS, Cs, Ts;

  if (!PYPARSETUPLE(args,
                    "OOOOOllllddddd", "OOOOOiiiiddddd",
                    "OOOOOllllfffff", "OOOOOiiiifffff",
                    &zonesR, &zoneD, &pyVariables, &pyParam_int,  &pyParam_real, &vartype, &compact,
                    &flagibc, &bctype, &gamma, &cv, &muS, &Cs, &Ts))
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
  E_Int flagIbc = E_Int(flagibc); 
  vector<PyArrayObject*> hook;
  E_Int imdjmd, imd,jmd,kmd, cnNfldD, nvars,ndimdxR, ndimdxD,meshtype;
  E_Float* iptroD; E_Float* iptroR; 

  E_Int nidom = PyList_Size(zonesR);

  E_Int* ipt_ndimdxR; E_Float** ipt_roR;
  ipt_ndimdxR = new E_Int[nidom];
  ipt_roR     = new E_Float*[nidom];

  /*-------------------------------------*/
  /* Extraction tableau int et real      */
  /*-------------------------------------*/
  FldArrayI* param_int;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real, true);
  E_Float* ipt_param_real = param_real->begin();


  E_Int nracID = ipt_param_int[0];
  E_Int nracIBC= ipt_param_int[1];
  E_Int nractot= nracID+nracIBC;

  E_Int nrac, shift_rac;
  if( flagIbc == 0) {nrac = nracID ; shift_rac = 2       ;}
  else              {nrac = nracIBC; shift_rac = 2+nracID;}

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;
  

  //Loc restrictif, tous les zones receveuses doivent avoir la meme loc
  loc = ipt_param_int[ 1+2*nractot + 6 ];
# include "getfromzoneDcompact.h"

  if( vartype <= 3 &&  vartype >= 1) nvars =5;
  else                               nvars =6;

  for  (E_Int irac=0; irac< nrac; irac++)
    { E_Int pos_rac   = ipt_param_int[ shift_rac +irac];
      E_Int no_zR     = ipt_param_int[ pos_rac   +  6 ];
      //printf("pos_rac %d %d %d %d \n", pos_rac, no_zR, shift_rac+irac, irac); 
      PyObject* zoneR = PyList_GetItem(zonesR, no_zR); // domaine i
#     include "getfromzoneRcompact.h"
      ipt_ndimdxR[irac] = ndimdxR;
      ipt_roR[irac]     = iptroR;
    }
 

  vector<E_Float*> vectOfRcvFields(nvars);
  vector<E_Float*> vectOfDnrFields(nvars);

  for  (E_Int irac=0; irac< nrac; irac++)
  {
    for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq] = ipt_roR[irac] + eq*ipt_ndimdxR[irac];
      vectOfDnrFields[eq] = iptroD        + eq*ndimdxD;
    }

    ////
    //  Interpolation parallele
    ////  
    //# include "commonInterpTransfers_indirect.h"
    ////  
    imdjmd = imd*jmd;
    E_Int max_thread = min(nvars , __NUMTHREADS__);

    # pragma omp parallel default(shared) num_threads(max_thread)
    {

#ifdef _OPENMP
     E_Int  ithread           = omp_get_thread_num()+1;
     E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
     E_Int ithread = 1;
     E_Int Nbre_thread_actif = 1;
#endif
     // Calcul du nombre de champs a traiter par chaque thread
     E_Int chunk = nvars/Nbre_thread_actif;
     E_Int r = nvars - chunk*Nbre_thread_actif;
     E_Int eq_deb, eq_fin;
     // equations traitees par thread
     if (ithread <= r) { eq_deb = (ithread-1)*(chunk+1)          ; eq_fin = eq_deb + (chunk+1); }  
     else              { eq_deb = (chunk+1)*r+(ithread-r-1)*chunk; eq_fin = eq_deb +  chunk;    }

     E_Int pos         = ipt_param_int[ shift_rac + nractot + irac];
     E_Float* ptrCoefs = ipt_param_real + pos;

     pos            = ipt_param_int[ shift_rac + irac];
     E_Int nbRcvPts = ipt_param_int[  pos + 1        ];
     E_Int nbDonPts = ipt_param_int[  pos            ];
  
     E_Int* types    = ipt_param_int +  pos + 7 + nbRcvPts + nbDonPts;
     E_Int* donorPts = ipt_param_int +  pos + 7             ;
     E_Int* rcvPts   = ipt_param_int +  pos + 7 +   nbDonPts;// donor et receveur inverser car storage donor
 
     //printf("nbRcvPts %d %d %d %d %d \n", nbRcvPts , types[0], irac, pos, pos + 7 + nbRcvPts + nbDonPts );
   
     E_Int indR, type;
     E_Int indD0, indD, i, j, k, ncfLoc, nocf;
     E_Int noi = 0; // compteur sur le tableau d indices donneur
     E_Int sizecoefs = 0;

     indR = rcvPts[0];
     //printf(" indR00 %d %d %d %d %d \n",nbRcvPts , nbRcvPts, nbDonPts, imd, jmd);
     //printf(" ndimdxR %d  %d  \n", ipt_ndimdxR[irac],ndimdxD );

     for (E_Int noind = 0; noind < nbRcvPts; noind++)
     {  
      //
      // adressage indirect pour indR
      //
       indR = rcvPts[noind];
#      include "commonInterpTransfers.h" 
       ptrCoefs += sizecoefs;
     } 
    }// omp

    if(flagIbc== 1)
    {
      E_Int threadmax_sdm  = __NUMTHREADS__;

      E_Int pos      = ipt_param_int[ shift_rac + irac];
      E_Int nbRcvPts = ipt_param_int[  pos + 1        ];
      E_Int nbDonPts = ipt_param_int[  pos            ];
      E_Int nbInterpD= ipt_param_int[  pos + 2        ];

      E_Int* rcvPts  = ipt_param_int +  pos + 7 +   nbRcvPts;// donor et receveur inverser car storage donor

            pos         = ipt_param_int[ shift_rac + nractot + irac];
      E_Float* ptrCoefs = ipt_param_real + pos;

      E_Int size = (nbRcvPts/threadmax_sdm)+1; // on prend du gras pour gerer le residus
      E_Int    r =  size % 8;
      if (r != 0) size  = size + 8 - r;        // on rajoute du bas pour alignememnt 64bits
      if (bctype <=1 ) size = 0;               // tableau inutile

      FldArrayF  tmp(size*13*threadmax_sdm);
      E_Float* ipt_tmp=  tmp.begin();

      E_Float* xPC     = ptrCoefs + nbInterpD;
      E_Float* xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
      E_Float* xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
      E_Float* densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;
#     pragma omp parallel default(shared)
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
        //  creer 2 zone  para pour threder au Max les loi de paroi
        if (varType == 1 || varType == 11)
		setIBCTransfersCommonVar1(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread, 
                             xPC, xPC+nbRcvPts, xPC     +nbRcvPts*2,
                                        xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                        xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                        densPtr, densPtr +nbRcvPts, densPtr +nbRcvPts*2, densPtr +nbRcvPts*3,
					ipt_tmp, size,
					gamma, cv, muS, Cs, Ts,
					vectOfDnrFields, vectOfRcvFields);
	else if (varType == 2 || varType == 21)
		setIBCTransfersCommonVar2(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
					xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                        xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                        xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                        densPtr, densPtr +nbRcvPts, densPtr +nbRcvPts*2, densPtr +nbRcvPts*3,
					ipt_tmp, size,
					gamma, cv, muS, Cs, Ts,
					vectOfDnrFields, vectOfRcvFields);
	else if (varType == 3 || varType == 31)
		setIBCTransfersCommonVar3(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
					xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                        xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                        xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                        densPtr, densPtr +nbRcvPts, densPtr +nbRcvPts*2, densPtr +nbRcvPts*3,
					ipt_tmp, size,
					gamma, cv, muS, Cs, Ts,
					vectOfDnrFields, vectOfRcvFields);

      } // Fin zone // omp
    } //ibc

  }//irac


 delete [] ipt_ndimdxR; delete [] ipt_roR;

 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );

 Py_INCREF(Py_None);
 return Py_None;
}
//=============================================================================
// Idem: in place + from zone + tc compact au niveau base 
//=============================================================================
PyObject* K_CONNECTOR::___setInterpTransfers(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *pyVariables;
  PyObject *pyParam_int, *pyParam_real;
  E_Int /*loc,*/ vartype, bctype, type_transfert, no_transfert, nitrun;
  E_Float gamma, cv, muS, Cs, Ts;

  if (!PYPARSETUPLE(args,
                    "OOOOOlllllddddd", "OOOOOiiiiiddddd",
                    "OOOOOlllllfffff", "OOOOOiiiiifffff",
                    &zonesR, &zonesD, &pyVariables, &pyParam_int,  &pyParam_real, &nitrun, &vartype,
                    &bctype, &type_transfert, &no_transfert, &gamma, &cv, &muS, &Cs, &Ts))
  {
      return NULL;
  }

  E_Int bcType = E_Int(bctype); // 0 : wallslip; 1: noslip; 2: log law of wall; 3: Musker law of wall
  E_Int NitRun = E_Int(nitrun);
  /* varType : 
     1  : conservatives, 
     11 : conservatives + ronutildeSA 
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA 
     3  : (ro,u,v,w,p)     
     31 : (ro,u,v,w,p) + ronutildeSA */
  E_Int varType = E_Int(vartype); 

  //gestion nombre de pass pour ID et/ou IBC
  E_Int TypeTransfert = E_Int(type_transfert);
  E_Int pass_deb, pass_fin;
  if     (TypeTransfert==0) { pass_deb =1; pass_fin =2; }//ID
  else if(TypeTransfert==1) { pass_deb =0; pass_fin =1; }//IBCD
  else                      { pass_deb =0; pass_fin =2; }//ALL


  E_Int NoTransfert   = E_Int(no_transfert);

  vector<PyArrayObject*> hook;
  //E_Int kmd, cnNfldD, nvars,ndimdxR, ndimdxD,meshtype;
  E_Int /*imd, jmd, imdjmd,*/ kmd, cnNfldD, nvars,/*ndimdxR, ndimdxD,*/meshtype;
  //E_Float* iptroD; E_Float* iptroR; 

  if( vartype <= 3 &&  vartype >= 1) nvars =5;
  else                               nvars =6;

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //pointeur pour stocker solution au centre ET au noeud 
  E_Int* ipt_ndimdxR; E_Int* ipt_ndimdxD; E_Int** ipt_cnd;
  E_Float** ipt_roR; E_Float** ipt_roD; E_Float** ipt_roR_vert;  E_Float** ipt_roD_vert;

  ipt_ndimdxR      = new E_Int[nidomR*2];   // on stocke ndimdx  en centre et vertexe
  ipt_roR          = new E_Float*[nidomR*2];
  ipt_roR_vert     = ipt_roR + nidomR;

  ipt_ndimdxD      = new E_Int[nidomD*8];  //on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  ipt_cnd          = new E_Int*[nidomD];
  ipt_roD          = new E_Float*[nidomD*2];
  ipt_roD_vert     = ipt_roD + nidomD;

  /*-------------------------------------*/
  /* Extraction tableau int et real      */
  /*-------------------------------------*/
  FldArrayI* param_int;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real, true);
  E_Float* ipt_param_real = param_real->begin();

  

  //On recupere le nom de la 1ere variable a recuperer 
   PyObject* tpl0= PyList_GetItem(pyVariables, 0); char* varname = PyString_AsString(tpl0);

  //on recupere sol et solcenter ainsi que connectivite et taille zones Donneuses (tc)
  for (E_Int nd = 0; nd < nidomD; nd++)
  {  
     PyObject* zoneD = PyList_GetItem(zonesD, nd);
#    include "getfromzoneDcompact_all.h"
  }

  //on recupere sol et solcenter et taille zones receuveuses (t)
  for (E_Int nd = 0; nd < nidomR; nd++)
  {  
     PyObject* zoneR = PyList_GetItem(zonesR, nd);
#    include "getfromzoneRcompact_all.h"
  }


  E_Int threadmax_sdm  = __NUMTHREADS__;
  E_Int ech            = ipt_param_int[ NoTransfert ];
  E_Int nrac           = ipt_param_int[ ech +1 ];          //nb total de raccord
  E_Int nrac_inst      = ipt_param_int[ ech +2 ];          //nb total de raccord instationnaire
  E_Int timelevel      = ipt_param_int[ ech +3 ];          //nb de pas de temps stocker pour chaque raccord instationnaire

  E_Int it_target = 0;
  if (timelevel !=0) { it_target = NitRun%timelevel; }//selection du pas de temps dans la base instationnaire 

  E_Int nrac_steady = nrac - nrac_inst;                     //nb total de raccord stationnaire

  // printf("nrac = %d, nrac_inst = %d, level= %d, it_target= %d  \n",  nrac, nrac_inst, timelevel,it_target);
  //on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx =0;
  for  (E_Int ipass_inst=0; ipass_inst< 2; ipass_inst++)
  {
    E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
    if(ipass_inst == 1){ irac_deb = ipt_param_int[ ech + 4 + it_target             ];
                         irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];  }

    for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
    { 
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
      if( ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1]; 
    }
  }

  E_Int size = (nbRcvPts_mx/threadmax_sdm)+1; // on prend du gras pour gerer le residus
  E_Int r =  size % 8;
  if (r != 0) size  = size + 8 - r;           // on rajoute du bas pour alignememnt 64bits
  if (bctype <=1 ) size = 0;                  // tableau inutile

  FldArrayF  tmp(size*13*threadmax_sdm);
  E_Float* ipt_tmp=  tmp.begin();

    //# pragma omp parallel default(shared)  num_threads(1)
    # pragma omp parallel default(shared)
    {

#ifdef _OPENMP
     E_Int  ithread           = omp_get_thread_num()+1;
     E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
     E_Int ithread = 1;
     E_Int Nbre_thread_actif = 1;
#endif

     E_Int indR, type;
     E_Int indD0, indD, i, j, k, ncfLoc/*, nocf*/, indCoef, noi, sizecoefs, /*Nbchunk,*/ imd, jmd, imdjmd;

     vector<E_Float*> vectOfRcvFields(nvars);
     vector<E_Float*> vectOfDnrFields(nvars);

  //1ere pass_typ: IBC
  //2eme pass_typ: transfert
  //
  for  (E_Int ipass_typ=pass_deb; ipass_typ< pass_fin; ipass_typ++)
  {
   // printf("ipass_typ = %d \n",  ipass_typ );
   //1ere pass_inst: les raccord fixe
   //2eme pass_inst: les raccord instationnaire
   for  (E_Int ipass_inst=0; ipass_inst< 2; ipass_inst++)
   {

    //printf("ipass_inst = %d, level= %d \n",  ipass_inst, nrac_inst_level );
    E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
    if(ipass_inst == 1){ irac_deb = ipt_param_int[ ech + 4 + it_target             ];
                         irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];  }

    //printf("iracdeb=  %d, iracfin= %d \n", irac_deb, irac_fin  );
    for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
    {
      //E_Int shift_rac =  ech + 4 + irac;
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

      //printf("ipass_typ = %d, ipass_inst= %d, irac=  %d, ithread= %d \n", ipass_typ,ipass_inst,irac , ithread );
      E_Int ibc =  ipt_param_int[ shift_rac + nrac*3     ];
      if( 1-ibc != ipass_typ)  continue;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
      E_Int loc      =  ipt_param_int[ shift_rac + nrac*9  +1 ]; //+1 a cause du nrac mpi
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ];
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ]; //neq fonction raccord rans/LES
      E_Int rotation =  ipt_param_int[ shift_rac + nrac*14 +1 ]; //flag pour periodicite azymuthal

      E_Int meshtype = ipt_ndimdxD[NoD + nidomD*6];
      E_Int cnNfldD  = ipt_ndimdxD[NoD + nidomD*7];
      E_Int* ptrcnd  = ipt_cnd[    NoD           ];

      if(loc == 0)
      {
        for (E_Int eq = 0; eq < nvars_loc; eq++)
        {
         vectOfRcvFields[eq] = ipt_roR_vert[ NoR] + eq*ipt_ndimdxR[ NoR + nidomR  ];
         vectOfDnrFields[eq] = ipt_roD_vert[ NoD] + eq*ipt_ndimdxD[ NoD + nidomD*3];
        }
         imd= ipt_ndimdxD[ NoD+ nidomD*4]; jmd= ipt_ndimdxD[ NoD + nidomD*5];
      }
      else
      {
        for (E_Int eq = 0; eq < nvars_loc; eq++)
        {
         vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_ndimdxR[ NoR ];
         vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_ndimdxD[ NoD ];
        }
         imd= ipt_ndimdxD[ NoD+ nidomD  ]; jmd= ipt_ndimdxD[ NoD+ nidomD*2];
      }

      imdjmd = imd*jmd;

      ////
      //  Interpolation parallele
      ////  
      ////  

       E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 + 1 ];
       //E_Int nbDonPts = ipt_param_int[ shift_rac                ];
  
       E_Int pos;
       pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int +  pos;
       pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int +  pos;
       pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts   = ipt_param_int +  pos;
       pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int +  pos;   // donor et receveur inverser car storage donor
       pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;

       E_Int nbInterpD = ipt_param_int[ shift_rac +  nrac ]; E_Float* xPC=NULL; E_Float* xPI=NULL; E_Float* xPW=NULL; E_Float* densPtr=NULL;
       if(ibc ==1)
       { 
        xPC     = ptrCoefs + nbInterpD;
        xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
        xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
        densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;
       }

       E_Int ideb      = 0;
       E_Int ifin      = 0;
       E_Int shiftCoef = 0;
       E_Int shiftDonnor  = 0;
       
       for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++)
       { 
        type      = types[ifin];

        SIZECF(type, meshtype, sizecoefs);
        ifin =  ifin + ntype[ 1 + ndtyp];

/*
// *      New school: meilleur equilibrage, mais gestion looop dynamique rame...
// *
        E_Int size_bc =  ifin-ideb;
        E_Int size_min=   16;
        //E_Int chunk = size_bc/Nbre_thread_actif;
        //if     (chunk < size_min && size_bc >= size_min) { chunk = size_min;}
        //else if(chunk < size_min && size_bc <  size_min) { chunk = size_bc ;}
        E_Int chunk = size_min;
        if(size_bc <  size_min) { chunk = size_bc ;}

        if      ( type == 0 ||  chunk <= 0) { Nbchunk = 1;                }
        else if ( chunk > 0)                { Nbchunk = size_bc/chunk;}

        chunk = size_bc/Nbchunk;

        E_Int r = size_bc - chunk*Nbchunk;

        #pragma omp for nowait schedule(dynamic,1)
        for (E_Int nd = 0; nd < Nbchunk; nd++)
        { 
*/  
        E_Int pt_deb, pt_fin;

/// oldschool
        // Calcul du nombre de champs a traiter par chaque thread
        E_Int size_bc =  ifin-ideb;
        E_Int chunk   =  size_bc/Nbre_thread_actif;
        E_Int r       =  size_bc - chunk*Nbre_thread_actif;
        // pts traitees par thread
        if (ithread <= r)
             { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); }
        else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk; } 

        //Si type 0, calcul sequentiel
        if      ( type == 0 )
          { if (ithread ==1 ){ pt_deb = ideb; pt_fin = ifin;}
            else             { pt_deb = ideb; pt_fin = ideb;}
          }

/// newschool suite
//        if (nd  <  r) { pt_deb = ideb + nd*(chunk+1)               ; pt_fin = pt_deb + (chunk+1); }  
//        else          { pt_deb = ideb +    (chunk+1)*r+(nd-r)*chunk; pt_fin = pt_deb +  chunk;    }

      //printf(" irac= %d, NoR= %d, nvar=  %d, NoD= %d, Rans=  %d, rot= %d, fin= %d, ithread= %d \n", irac, NoR, nvars_loc, NoD, ipass_inst ,rotation, pt_fin , ithread );
      //if(ithread <=8 && NoD==83 )  printf(" shift %d  %d %d %d  %d %d %d  %d \n", irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  , type, ithread );
      //if(ithread <=8 && NoR==114 )  printf(" new   %d  %d %d %d  %d %d %d  %d \n", irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  , type, ithread );

          noi       = shiftDonnor;                             // compteur sur le tableau d indices donneur
          indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;

          if     (nvars_loc==5)
          {
#           include "commonInterpTransfers_reorder_5eq.h" 
          }
          else if(nvars_loc==6)
          {
#           include "commonInterpTransfers_reorder_6eq.h" 
          }
          else
          {
#           include "commonInterpTransfers_reorder_neq.h" 
          }
           

          // Prise en compte de la periodicite par rotation
          if ( rotation == 1 )
          {
           E_Float* angle = ptrCoefs + nbInterpD;
#          include "includeTransfers_rotation.h"
          }

          // ibc    
          if(ibc== 1)
          {
            if (varType == 1 || varType == 11)
              setIBCTransfersCommonVar1(bcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                        xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                        xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                        xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                        densPtr, densPtr +nbRcvPts, densPtr +nbRcvPts*2, densPtr +nbRcvPts*3,
                                        ipt_tmp, size,
                                        gamma, cv, muS, Cs, Ts,
                                        vectOfDnrFields, vectOfRcvFields);
            else if (varType == 2 || varType == 21)
              setIBCTransfersCommonVar2(bcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                        xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                        xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                        xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                        densPtr, densPtr +nbRcvPts, densPtr +nbRcvPts*2, densPtr +nbRcvPts*3,
                                        ipt_tmp, size,
                                        gamma, cv, muS, Cs, Ts,
                                        vectOfDnrFields, vectOfRcvFields);
            else if (varType == 3 || varType == 31)
              setIBCTransfersCommonVar3(bcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                        xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                        xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                        xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                        densPtr, densPtr +nbRcvPts, densPtr +nbRcvPts*2, densPtr +nbRcvPts*3,
                                        ipt_tmp, size,
                                        gamma, cv, muS, Cs, Ts,
                                        vectOfDnrFields, vectOfRcvFields);
          }//ibc          
  //*
  //        } //chunk
  //*/
          ideb       = ideb + ifin;
          shiftCoef  = shiftCoef    +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
          shiftDonnor= shiftDonnor  +  ntype[1+ndtyp];           //shift donnor entre 2 types successif
       }// type 
    }//irac
   }//ipass_inst
  #pragma omp barrier 
  }//ipass
  }// omp

 delete [] ipt_ndimdxR; delete [] ipt_roR; delete [] ipt_ndimdxD; delete [] ipt_roD; delete [] ipt_cnd;

 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );

 Py_INCREF(Py_None);
 return Py_None;
}
