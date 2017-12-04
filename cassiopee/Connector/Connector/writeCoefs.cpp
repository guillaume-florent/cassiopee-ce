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
# include <stdio.h>
#include <map>
# include "connector.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Ecriture des coefficients d'interpolation dans un fichier relisible 
   par elsA */
//=============================================================================
PyObject* K_CONNECTOR::writeCoefs(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject *RcvIndexMap, *DonorIndexMap;  // indices des pts interpoles et des cellules donneuses
  PyObject *EXDirectionMap;
  PyObject *DonorInterpCoefMap;           // coefficients d interpolation
  PyObject *DonorCellNMap;                // champs cellN des zones d interpolations
  PyObject *InterpTypesMap;               // types des interpolations
  PyObject *ZoneDimMap;                   // dimensions des grilles aux centres des cellules pour chaque zone
  PyObject *BlockRcvIdMap;                // Id des zones interpolees
  PyObject *Npts;                         // nombre total de points interpoles (dans tous les blocs) par ne zone d interpolation
  char* PrefixFile;                       // prefixe pour le nommage des fichiers elsA
  E_Int isEX = 0;                      // isEX = 0 ou 1 : Chimere avec 2 ou 1 rangees de cellules fictives
  E_Int NZones;                        // nombre de points interpoles
  E_Int Solver;                        // solveur pour lequel les fichiers sont ecrits
  E_Int NGhostCells;
  if (!PYPARSETUPLEI(args,
                    "lOOOOOOOOOslll", "iOOOOOOOOOsiii",
                    &NZones, &BlockRcvIdMap, &RcvIndexMap, 
                    &EXDirectionMap, &DonorIndexMap, &DonorInterpCoefMap, &InterpTypesMap,
                    &DonorCellNMap, &ZoneDimMap,
                    &Npts, &PrefixFile, &isEX, &Solver, &NGhostCells))
  {
      return NULL;
  }

  /*-------------------*/
  /* Variables locales */
  /*-------------------*/
  E_Int intKey;

  /*-------------------------------*/
  /* Extraction du solveur */
  /*-------------------------------*/
  E_Int solver = Solver; //1: elsA, 2: Cassiopee

  /*-------------------------------*/
  /* Extraction du nombre de zones */
  /*-------------------------------*/
  E_Int nzones = NZones;

  /*-------------------------------*/
  /* Nb de ghost cells */
  /*-------------------------------*/
  E_Int nbOfGc = NGhostCells;

  /*-------------------------------------------------------------*/
  /* Extraction du nombre de points interpoles par zone donneuse */
  /*-------------------------------------------------------------*/
  vector<E_Int> nbInterpCells;
  PyObject* tpl;
  E_Int ind;
  E_Int size = PyList_Size(Npts);
  for (E_Int v  = 0 ; v < size; v++)
  {
    tpl = PyList_GetItem(Npts, v);
    if (PyInt_Check(tpl) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable %d. Skipped...\n", v);
    }
    else
    {
      ind = PyInt_AsLong(tpl);
      nbInterpCells.push_back(ind);
    }
  }

  /*------------------------------------------------*/
  /* Extraction des listes d'Id des blocs receveurs */
  /*------------------------------------------------*/
  map<E_Int, vector<E_Int> > blockRcvIdMap;
  PyObject* pyListValues;
  if (PyDict_Check (BlockRcvIdMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: BlockRcvIdMap must be a dictionary.");
    return NULL;
  }
  PyObject * key_list = PyDict_Keys(BlockRcvIdMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable %d. Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int>
    vector<E_Int> tmpListValue;
    pyListValues = PyDict_GetItem(BlockRcvIdMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: BlockRcvIdMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* pyIntValue = PyList_GetItem(pyListValues, v);
      if (PyInt_Check(pyIntValue) == 0)
      {
        printf("Warning: writeCoefs: invalid int value in  BlockRcvIdMap\n");
      }
      else
      {
        E_Int intValue = PyInt_AsLong(pyIntValue);
        tmpListValue.push_back(intValue);
      }
    }
   blockRcvIdMap[intKey] = tmpListValue;
  }

  /*--------------------------------------------------------*/
  /* Extraction des listes de tableaux d indices            */
  /*--------------------------------------------------------*/
  map<E_Int,vector<E_Int*> > rcvIndexMap;
  map<E_Int,vector<E_Int*> > donorIndexMap;

  // rcvIndexMap : indices du bloc interpole
  if (PyDict_Check (RcvIndexMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: RcvIndexMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(RcvIndexMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable %d. Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int*>
    vector<E_Int*> tmpListValue;
    pyListValues = PyDict_GetItem(RcvIndexMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: RcvIndexMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* intArray = PyList_GetItem(pyListValues, v);
      PyArrayObject* a = (PyArrayObject*)intArray;
      if (PyArray_Check(a) == 0)
      {
      printf("Warning: writeCoefs: RcvIndexMap must contain a list of arrays\n");
      }
      else
      {
        E_Int* indArray = (E_Int*)PyArray_DATA(a); 
        tmpListValue.push_back(indArray);
      }
    }
   rcvIndexMap[intKey] = tmpListValue;
  }

  // DonorIndexArray  : indices du bloc donneur
  if (PyDict_Check (DonorIndexMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: DonorIndexMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(DonorIndexMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable %d. Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int*>
    vector<E_Int*> tmpListValue;
    pyListValues = PyDict_GetItem(DonorIndexMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: DonorIndexMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* intArray = PyList_GetItem(pyListValues, v);
      PyArrayObject* a = (PyArrayObject*)intArray;
      if (PyArray_Check(a) == 0)
      {
      printf("Warning: writeCoefs: DonorIndexMap must contain a list of arrays\n");
      }
      else
      { 
        E_Int* indArray = (E_Int*)PyArray_DATA(a); tmpListValue.push_back(indArray);
        //int* indArray = (int*)PyArray_DATA(a); tmpListValue.push_back(indArray);
      }
    }
   donorIndexMap[intKey] = tmpListValue;
  }

  /*----------------------------------------------------------*/
  /* Extraction des tableaux d indirection pour les points EX */
  /*----------------------------------------------------------*/
  map<E_Int,vector<E_Int*> >directionEXMap;
  if (isEX)
  {    
    if (PyDict_Check (EXDirectionMap) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: EXDirectionMap must be a dictionary.");
      return NULL;
    }
    key_list = PyDict_Keys(EXDirectionMap);
    size = PyList_Size(key_list);
    //Get keys and corresponding values from the dictionary
    for (E_Int i  = 0 ; i < size; i++)
    {
      PyObject * pyKey = 0;
      pyKey = PyList_GetItem(key_list, i);
      if (PyInt_Check(pyKey) == 0)
      {
        printf("Warning: writeCoefs: invalid int for variable %d. Skipped...\n", i);
      }
      else
        intKey = PyInt_AsLong(pyKey);
      //Convert to a C++ vector<E_Int*>
      vector<E_Int*> tmpListValue;
      pyListValues = PyDict_GetItem(EXDirectionMap, pyKey);
      // check if pyListValue is a list
      if (PyList_Check (pyListValues) == 0)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "writeCoefs: EXDirectionMap must contain lists.");
        return NULL;
      }
      // fill C++ map
      E_Int tmpListValueSize = PyList_Size(pyListValues);
      for (E_Int v  = 0 ; v < tmpListValueSize; v++)
      {
        PyObject* intArray = PyList_GetItem(pyListValues, v);
        PyArrayObject* a = (PyArrayObject*)intArray;
        if (PyArray_Check(a) == 0)
        {
          printf("Warning: writeCoefs: EXDirectionMap must contain a list of arrays\n");
        }
        else
        {
          E_Int* indArray = (E_Int*)PyArray_DATA(a); 
          tmpListValue.push_back(indArray);
        }
      }
      directionEXMap[intKey] = tmpListValue;
    }
  }

  /*--------------------------------------------------------*/
  /* Extraction des listes de coefficients d interpolation  */
  /*--------------------------------------------------------*/
  // DonorInterpCoef : coefficient d'interpolation
  map<E_Int,vector<FldArrayF> > donorCoefMap;
  if (PyDict_Check (DonorInterpCoefMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: DonorInterpCoefMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(DonorInterpCoefMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable %d. Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<FldArrayF>
    vector<FldArrayF> tmpListValue;
    pyListValues = PyDict_GetItem(DonorInterpCoefMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: DonorInterpCoefMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* pyArrayValue = PyList_GetItem(pyListValues, v);
      PyArrayObject* a = (PyArrayObject*)pyArrayValue;
      if (PyArray_Check(a) == 0)
      {
      printf("Warning: writeCoefs: DonorInterpCoefMap must contain a list of arrays\n");
      }
      else
      {
#ifdef NPY_1_7_API_VERSION
        E_Int isFortran = PyArray_IS_F_CONTIGUOUS(a);
#else
        E_Int isFortran = PyArray_CHKFLAGS(a, NPY_F_CONTIGUOUS);
#endif
        E_Int adim1, adim2;
        if (isFortran == 1) {adim1 = PyArray_DIMS(a)[0]; adim2 = PyArray_DIMS(a)[1];}
        else  {adim1 = PyArray_DIMS(a)[1]; adim2 = PyArray_DIMS(a)[0];}
        FldArrayF coefArray(adim1, adim2, (E_Float*)PyArray_DATA(a));
        
        tmpListValue.push_back(coefArray);
      }
    }
    donorCoefMap[intKey] = tmpListValue;
  }

  /*--------------------------------------------------------*/
  /* Extraction des types d interpolation                   */
  /*--------------------------------------------------------*/
  map<E_Int,vector<E_Int*> > interpTypesMap;
  if (PyDict_Check (InterpTypesMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: InterpTypesMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(InterpTypesMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i=0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
      printf("Warning: writeCoefs: invalid integer for variable %d. Skipped...\n", i);
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int*>
    vector<E_Int*> tmpListValue;
    pyListValues = PyDict_GetItem(InterpTypesMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: InterpTypesMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* intArray = PyList_GetItem(pyListValues, v);
      PyArrayObject* a = (PyArrayObject*)intArray;
      if (PyArray_Check(a) == 0)
      {
        printf("Warning: writeCoefs: InterpTypesMap must contain a list of arrays.\n");
      }
      else
      {
        E_Int* indArray = (E_Int*)PyArray_DATA(a); 
        tmpListValue.push_back(indArray);
      }
    }
    interpTypesMap[intKey] = tmpListValue;
  }
  /*--------------------------------------------------------*/
  /* Extraction des listes de champs cellN                  */
  /*--------------------------------------------------------*/
  // DonorCellN : champs cellNatureField
  map<E_Int,FldArrayF> lDonorCellN;
  if (PyDict_Check (DonorCellNMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: DonorCellNMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(DonorCellNMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable %d. Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    // fill C++ map
    PyObject* pyArrayValue = PyDict_GetItem(DonorCellNMap, pyKey);
    PyArrayObject* a = (PyArrayObject*)pyArrayValue;
    if (PyArray_Check(a) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: lDonorCellN must contain a list of arrays.");
      return NULL;
    }
    FldArrayF cellnArray(PyArray_DIMS(a)[1], PyArray_DIMS(a)[0], (E_Float*)PyArray_DATA(a));
    lDonorCellN[intKey] = cellnArray;
  }

  /*--------------------------------------------------------*/
  /* Extraction des listes de dimensions des blocs donneurs */
  /*--------------------------------------------------------*/
  map<E_Int, vector<E_Int> > zoneDimMap;
  if (PyDict_Check (ZoneDimMap) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "writeCoefs: ZoneDimMap must be a dictionary.");
    return NULL;
  }
  key_list = PyDict_Keys(ZoneDimMap);
  size = PyList_Size(key_list);
  //Get keys and corresponding values from the dictionary
  for (E_Int i  = 0 ; i < size; i++)
  {
    PyObject * pyKey = 0;
    pyKey = PyList_GetItem(key_list, i);
    if (PyInt_Check(pyKey) == 0)
    {
      printf("Warning: writeCoefs: invalid int for variable %d. Skipped...\n", i);
    }
    else
      intKey = PyInt_AsLong(pyKey);
    //Convert to a C++ vector<E_Int>
    vector<E_Int> tmpListValue;
    pyListValues = PyDict_GetItem(ZoneDimMap, pyKey);
    // check if pyListValue is a list
    if (PyList_Check (pyListValues) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "writeCoefs: ZoneDimMap must contain lists.");
      return NULL;
    }
    // fill C++ map
    E_Int tmpListValueSize = PyList_Size(pyListValues);
    for (E_Int v  = 0 ; v < tmpListValueSize; v++)
    {
      PyObject* pyIntValue = PyList_GetItem(pyListValues, v);
      if (PyInt_Check(pyIntValue) == 0)
      {
        printf("Warning: writeCoefs: invalid int value in  ZoneDimMap\n");
      }
      else
      {
        E_Int intValue = PyInt_AsLong(pyIntValue);
        tmpListValue.push_back(intValue);
      }
    }
   zoneDimMap[intKey] = tmpListValue;
  }

  /*--------------------------------------------------------*/
  /* Ecriture des fichiers d interpolation                  */
  /*--------------------------------------------------------*/
  E_Int typeInterp;
  // Parcours des zones d interpolation
  for (E_Int noz = 0; noz < nzones; noz++)
  {

    // Id du bloc donneur
    E_Int BlockDonorId = noz;

    // Donnees pour la zone noz
    E_Int nbRcvZones = rcvIndexMap[noz].size();

    // 1- Fichier de cellN (necessaire pour une relecture de fichiers des coefficients dans elsA)
    // REMARQUES : - le fichier cellN doit contenir les cellules fictives!
    //             - on peut mettre cellN a 1 sur ces cellules fictives, car la valeur sera ensuite  
    //             ecrasee dans le traitement des raccords
    if (!isEX)
    {
      E_Int ni = zoneDimMap[noz][0]; E_Int nj = zoneDimMap[noz][1]; E_Int nk = zoneDimMap[noz][2];
      E_Int nig = ni+2*nbOfGc; E_Int njg = nj+2*nbOfGc; E_Int nkg = nk+2*nbOfGc;
      E_Int nignjg = nig*njg; E_Int ninj = ni*nj;
      FldArrayI cellNgc;
      E_Float eps = K_CONST::E_GEOM_CUTOFF;
      if (solver == 1) // elsA : prise en compte des cellules fictives pour l indicage
      {
        cellNgc.malloc(nig*njg*nkg); cellNgc.setAllValuesAt(1);
        for (E_Int k=0; k<nk; k++)
          for (E_Int j=0; j<nj; j++)
            for (E_Int i=0; i<ni; i++)
            {
              E_Int ind = k*ninj + j*ni +i;
              E_Int indg = (k+nbOfGc)*nignjg + (j+nbOfGc)*nig +i+nbOfGc;
              if (K_FUNC::fEqualZero(lDonorCellN[noz][ind],eps) == true)
                cellNgc[indg] = -1;
              else if (K_FUNC::fEqualZero(lDonorCellN[noz][ind] - 2.,eps) == true)
                cellNgc[indg] = 0;
            }
      }
      else // Cassiopee : pas de cellules fictives
      {
        cellNgc.malloc(ni*nj*nk); cellNgc.setAllValuesAt(1);
         for (E_Int k=0; k<nk; k++)
          for (E_Int j=0; j<nj; j++)
            for (E_Int i=0; i<ni; i++)
            {
              E_Int ind = k*ninj + j*ni +i;
              if (K_FUNC::fEqualZero(lDonorCellN[noz][ind],eps) == true)
                cellNgc[ind] = -1;
              else if (K_FUNC::fEqualZero(lDonorCellN[noz][ind] - 2.,eps) == true)
                cellNgc[ind] = 0;
            }       
      }
      // Ouverture du fichier
      FILE* ptr_file = NULL;
      char* file = new char[K_ARRAY::VARSTRINGLENGTH];
      strcpy(file,PrefixFile);
      char* strId = new char[K_ARRAY::VARSTRINGLENGTH]; sprintf(strId,"%04d",BlockDonorId);
      strcat(file,strId);
      strcat(file,"_Blanking");
      ptr_file = fopen(file, "w");
      printf("Open file %s\n",file);fflush(stdout);
      // Ecriture du nombre de points du domaine d interpolation
      E_Int nptsInterp = cellNgc.getSize();
      fprintf(ptr_file,"%d\n",nptsInterp);
      for (E_Int np = 0; np < nptsInterp; np++)
      {
        fprintf(ptr_file, "%d ",
                cellNgc[np]);
      }
      fclose(ptr_file);
    }

    // 1- Fichier d'interpolations
    // Ouverture du fichier
    FILE* ptr_file = NULL;
    char* file = new char[K_ARRAY::VARSTRINGLENGTH];
    strcpy(file,PrefixFile);
    char* strId = new char[K_ARRAY::VARSTRINGLENGTH]; sprintf(strId,"%04d",BlockDonorId);
    strcat(file,strId);
    if (isEX) strcat(file,"_Int");
    ptr_file = fopen(file, "w");
    printf("Open file %s\n",file);fflush(stdout);
    // Ecriture du nombre de points d interpolations
    E_Int npts = nbInterpCells[noz];
    fprintf(ptr_file,"%d\n",npts);

    for (E_Int n = 0; n < nbRcvZones; n++)
    {
      E_Int nbOfDatas = donorCoefMap[noz][n].getSize();

      FldArrayF cellN = lDonorCellN[noz];
      
      E_Int blockRcvId = blockRcvIdMap[noz][n];
      E_Int RcvId; 
      // pointeurs sur le champ donorCoefMap de coefficients d'interpolation
      E_Float* donorCoefMap1 = donorCoefMap[noz][n].begin(1); E_Float* donorCoefMap2 = donorCoefMap[noz][n].begin(2); E_Float* donorCoefMap3 = donorCoefMap[noz][n].begin(3); 
      E_Float* donorCoefMap4 = donorCoefMap[noz][n].begin(4); E_Float* donorCoefMap5 = donorCoefMap[noz][n].begin(5); E_Float* donorCoefMap6 = donorCoefMap[noz][n].begin(6); 
      E_Float* donorCoefMap7 = donorCoefMap[noz][n].begin(7); E_Float* donorCoefMap8 = donorCoefMap[noz][n].begin(8);
      
      // Ecriture des donnees sous la forme :
      // indice_pt_interp indice_cell_donneuse id_bloc_rcv coef_interp_1 coef_interp_2 ... coef_interp_8
      for (E_Int np = 0; np < nbOfDatas; np++)
      {
        // indice pour la cellule donneuse exprimee dans le maillage en centres etendus
        E_Int indDonor = (E_Int) donorIndexMap[noz][n][np];
        // indice pour la cellule interpolee
        E_Int indRcv;
        E_Int ni = zoneDimMap[blockRcvId][0]; E_Int nj = zoneDimMap[blockRcvId][1]; E_Int nk = zoneDimMap[blockRcvId][2];
        E_Int ninj = ni*nj;
        E_Int kd  = rcvIndexMap[noz][n][np]/ninj; E_Int val = rcvIndexMap[noz][n][np] - kd*ninj;
        E_Int jd =  val/ni;
        E_Int id = val - jd*ni;
        E_Int nig = ni+2*nbOfGc; E_Int njg = nj+2*nbOfGc; E_Int nkg = nk+2*nbOfGc;
        E_Int nignjg = nig*njg; 
        E_Int kg = kd + nbOfGc; E_Int jg = jd + nbOfGc; E_Int ig = id + nbOfGc;
        typeInterp = (E_Int)interpTypesMap[noz][n][np];
        if (solver == 1) // elsA : prise en compte des cellules fictives pour l indicage
        {
          indRcv = kg*nignjg + jg*nig +ig;
        }
        else // Cassiopee : pas de cellules fictives
        {
          indRcv = kd*ni*nj + jd*ni +id;
        }
        if (solver == 1)  RcvId = blockRcvId; // elsA : les Id des blocs commencent a 0
        else RcvId = blockRcvId+1;// Cassiopee : les Id des blocs commencent a 1
        if (!isEX)
        {
          if (solver == 1) // elsA/Kernel
            fprintf(ptr_file, "%d %d %d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",
                    indDonor, indRcv, RcvId, typeInterp, 
                    donorCoefMap1[np], donorCoefMap2[np],donorCoefMap3[np], donorCoefMap4[np], 
                    donorCoefMap5[np], donorCoefMap6[np], donorCoefMap7[np], donorCoefMap8[np]);
          else //Cassiopee/Kernel
            fprintf(ptr_file, "%d %d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",
                    indDonor, indRcv, RcvId, donorCoefMap1[np], donorCoefMap2[np],donorCoefMap3[np], 
                    donorCoefMap4[np], donorCoefMap5[np], donorCoefMap6[np], donorCoefMap7[np], donorCoefMap8[np]);
        }
        else
        {
          // Calcul de l indice de l interface pour le point EX
          E_Int indIntEX=0;
          if (solver == 1) // elsA : prise en compte des cellules fictives pour l indicage
          {
            RcvId = blockRcvId; // elsA : les Id des blocs commencent a 0
            E_Int kg = kd + nbOfGc; E_Int jg = jd + nbOfGc; E_Int ig = id + nbOfGc;
            E_Int nbIntByDir = nignjg*nkg;
            if      (directionEXMap[noz][n][np] == 0) indIntEX = kg*nignjg+jg*nig+ig+1;                  // interface de frontiere Imax
            else if (directionEXMap[noz][n][np] == 1) indIntEX = kg*nignjg+jg*nig+ig;                    // interface de frontiere Imin
            else if (directionEXMap[noz][n][np] == 2) indIntEX = kg*nignjg+(jg+1)*nig+ig +   nbIntByDir; // interface de frontiere Jmax
            else if (directionEXMap[noz][n][np] == 3) indIntEX = kg*nignjg+jg*nig+ig     +   nbIntByDir; // interface de frontiere Jmin
            else if (directionEXMap[noz][n][np] == 4) indIntEX = (kg+1)*nignjg+jg*nig+ig + 2*nbIntByDir; // interface de frontiere Kmax
            else if (directionEXMap[noz][n][np] == 5) indIntEX = kg*nignjg+jg*nig+ig     + 2*nbIntByDir; // interface de frontiere Kmin
            fprintf(ptr_file, "%d %d %d %d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",
                    indDonor, indIntEX, directionEXMap[noz][n][np], RcvId, typeInterp, 
                    donorCoefMap1[np], donorCoefMap2[np],donorCoefMap3[np], donorCoefMap4[np], 
                    donorCoefMap5[np], donorCoefMap6[np], donorCoefMap7[np], donorCoefMap8[np]);
          }
          else // Cassiopee/Kernel : pas de cellules fictives
          {
            RcvId = blockRcvId+1;// Cassiopee : les Id des blocs commencent a 1
            E_Int  nbIntByDiri = (ni+1)*nj*nk;
            E_Int nbIntByDirj = ni*(nj+1)*nk;
            if      (directionEXMap[noz][n][np] == 0) indIntEX = kd*(ni+1)*nj+jd*(ni+1)+id+1;                               // interface de frontiere Imax
            else if (directionEXMap[noz][n][np] == 1) indIntEX = kd*(ni+1)*nj+jd*(ni+1)+id;                                 // interface de frontiere Imin
            else if (directionEXMap[noz][n][np] == 2) indIntEX = kd*ni*(nj+1)+(jd+1)*ni+id + nbIntByDiri;               // interface de frontiere Jmax
            else if (directionEXMap[noz][n][np] == 3) indIntEX = kd*ni*(nj+1)+jd*ni+id     + nbIntByDiri;               // interface de frontiere Jmin
            else if (directionEXMap[noz][n][np] == 4) indIntEX = (kd+1)*ninj+jd*ni+id + nbIntByDiri + nbIntByDirj; // interface de frontiere Kmax
            else if (directionEXMap[noz][n][np] == 5) indIntEX = kd*ninj+jd*ni+id     + nbIntByDiri + nbIntByDirj; // interface de frontiere Kmin      
            
            fprintf(ptr_file, "%d %d %d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",
                    indDonor, indIntEX, directionEXMap[noz][n][np], RcvId, donorCoefMap1[np], donorCoefMap2[np],donorCoefMap3[np], 
                    donorCoefMap4[np], donorCoefMap5[np], donorCoefMap6[np], donorCoefMap7[np], donorCoefMap8[np]);
          }          
        }
      }
    }
    fclose(ptr_file);
  }

  /*-------------------------------*/
  /* Destruction des objets crees  */
  /*-------------------------------*/
//   // destruction de rcvIndexMap
//   map<E_Int,vector<E_Int*> >::iterator itr;
//   for (itr=rcvIndexMap.begin(); itr != rcvIndexMap.end();itr++)
//   {
//     vector<E_Int*> vect = (*itr).second; 
//     for (E_Int w = 0; w < vect.size(); w++) delete [] vect[w]; 
//   }
//   // destruction de donorIndexMap
//   for (itr=donorIndexMap.begin(); itr != donorIndexMap.end();itr++)
//   {
//     vector<E_Int*> vect = (*itr).second; 
//     for (E_Int w = 0; w < vect.size(); w++) delete [] vect[w]; 
//   }

  Py_INCREF(Py_None);
  return Py_None;
}
