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

// Binary ADF CGNS file support

// Comportement de cette routine:
// A l'ecriture: on respecte les types d'entree, sauf pour les noeuds version
// et periodic (R4)
// A la lecture: on converti en numpy I4/R8.

# include "GenIO.h"
# include "kcore.h"
# include "ADF/ADF.h"
# define CGNSMAXLABEL 32
# define CGNSMAXDIM 20
using namespace std;

namespace K_IO
{
class GenIOAdf
{
  public:
    PyObject* loadOne(PyObject* tree, int depth);
    std::list<double>* getChildren(double parent);
    void getName(double node, char* name);
    void getLabel(double node, char* label);
    void getType(double node, char* type, int dim, int* dims);
    PyObject* createNode(double node);
    double writeNode(double node, PyObject* tree);
    int getSingleI4(double node);
    E_LONG getSingleI8(double node);
    float getSingleR4(double node);
    double getSingleR8(double node);
    PyObject* getArrayR8(double node);
    PyObject* getArrayR8Skel(double node);
    PyObject* getArrayR4(double node);
    PyObject* getArrayR4Skel(double node);
    PyObject* getArrayI4(double node);
    PyObject* getArrayI4Skel(double node);
    PyObject* getArrayI8(double node);
    PyObject* getArrayI8Skel(double node);
    char* getArrayC1(double node);
    double setSingleR4(double node, float data);
    double setSingleR8(double node, double data);
    double setSingleI4(double node, int data);
    double setSingleI8(double node, E_LONG data);
    double setArrayMT(double node);
    double setArrayR4(double node, float *data, int dim, int *dims);
    double setArrayR8(double node, double *data, int dim, int *dims);
    double setArrayI4(double node, int *data, int dim, int *dims);
    double setArrayI8(double node, E_LONG *data, int dim, int *dims);
    double setArrayC1(double node, char *data);
    PyObject* dumpOne(PyObject* tree);
    void blankAndCopy(char *s, const char *p, E_Int size);
    GenIOAdf() {_skeleton=0; _maxFloatSize=1e6; _maxDepth=1e6;};
  public:
    std::list<double> _fatherStack;
    int _errorFlag;
    int _skeleton;
    int _maxFloatSize;
    int _maxDepth;

    char _type[CGNSMAXLABEL+1];
    char _name[CGNSMAXLABEL+1];
    char _label[CGNSMAXLABEL+1];
    char _dtype[CGNSMAXLABEL+1];
    int _dims[CGNSMAXDIM];
};
}
//=============================================================================
/*
   adfcgnsread
   IN: file: nom du fichier a liste
   OUT: tree: le pyTree lu
   IN: skeleton: si 1, on ne lit que le squelette de l'arbre (les numpy
   float ne sont pas lus).
*/
//=============================================================================
E_Int K_IO::GenIO::adfcgnsread(char* file, PyObject*& tree, int skeleton,
                               int maxFloatSize, int maxDepth)
{
  /* Open file */
  tree = PyList_New(4);
  double rootId; int errorFlag;
  ADF_Database_Open(file, "READ_ONLY", "NATIVE", &rootId, &errorFlag);
  if (errorFlag > 0)
  {
    printf("Warning: adfcgnsread: cannot open file %s.\n", file);
    return 1;
  }

  /* Recursive load */
  GenIOAdf ADF;
  ADF._skeleton = skeleton;
  ADF._maxFloatSize = maxFloatSize;
  if (maxDepth == -1) ADF._maxDepth = 1e6;
  else ADF._maxDepth = maxDepth;
  ADF._fatherStack.push_front(rootId);
  ADF.loadOne(tree, 0);
  ADF._fatherStack.pop_front();

  /* Close */
  ADF._fatherStack.clear();
  ADF_Database_Close(rootId, &errorFlag);

  return 0;
}

//=============================================================================
PyObject* K_IO::GenIOAdf::loadOne(PyObject* tree, int depth)
{
  list<double>* sonList;
  PyObject* node;
  PyObject* l = PyList_New(0);
  PyList_SetItem(tree, 2, l);
  if (depth >= _maxDepth) return tree;
  sonList = getChildren(_fatherStack.front());
  for (list<double>::iterator it = sonList->begin();
       it != sonList->end(); it++)
  {
    node = createNode(*it);
    PyList_Append(l, node); Py_DECREF(node);
    _fatherStack.push_front(*it);
    loadOne(node, depth+1);
    _fatherStack.pop_front();
  }
  delete sonList;
  return tree;
}

//=============================================================================
list<double>* K_IO::GenIOAdf::getChildren(double parent)
{
  list<double> *r; // a desallouer a l'exterieur
  double p = parent;
  int lr;
  int nChildren;
  double id;
  ADF_Number_of_Children(p, &nChildren, &_errorFlag);
  r = new list<double>;
  for (int i = 1; i <= nChildren; i++)
  {
    ADF_Children_Names(p, i, 1, CGNSMAXLABEL+1,
                       &lr, (char*)_label, &_errorFlag);
    ADF_Get_Node_ID(p, _label, &id, &_errorFlag);
    r->push_back(id);
  }
  return r;
}

//=============================================================================
// Fills _name
//==============================================================================
void K_IO::GenIOAdf::getName(double node, char* name)
{
  ADF_Get_Name(node, name, &_errorFlag);
  // si errorFlag == 4: depassement
}

//=============================================================================
// Fill _label
//==============================================================================
void K_IO::GenIOAdf::getLabel(double node, char* label)
{
  ADF_Get_Label(node, label, &_errorFlag);
  // si errorFlag == 4: depassement
}

//=============================================================================
void K_IO::GenIOAdf::getType(double node, char* type, int dim, int* dims)
{
  ADF_Get_Data_Type(node, type, &_errorFlag);
  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, dims, &_errorFlag);
}

//=============================================================================
// Cree un noeud du pyTree
//=============================================================================
PyObject* K_IO::GenIOAdf::createNode(double node)
{
  // Nom et type du noeud de l'arbre
  getLabel(node, _type);
  getName(node, _name);

  // Valeur du noeud
  ADF_Get_Data_Type(node, _dtype, &_errorFlag);
  npy_intp npy_dim_vals2[1]; npy_dim_vals2[0] = 0;

  PyObject* v = NULL; // set everything in numpys
  if (strcmp(_dtype, "I4") == 0)
  {
    if (_skeleton == 1 && (strcmp(_type, "IndexArray_t") == 0 || strcmp(_type, "DataArray_t") == 0)) v = getArrayI4Skel(node);
    else v = getArrayI4(node);
  }
  else if (strcmp(_dtype, "I8") == 0)
  {
    if (_skeleton == 1 && (strcmp(_type, "IndexArray_t") == 0 || strcmp(_type, "DataArray_t") == 0)) v = getArrayI8Skel(node);
    else v = getArrayI8(node);
  }
  else if (strcmp(_dtype, "R4") == 0)
  {
    if (_skeleton == 1 && strcmp(_type, "DataArray_t") == 0) v = getArrayR4Skel(node);
    else v = getArrayR4(node);
  }
  else if (strcmp(_dtype, "R8") == 0)
  {
    if (_skeleton == 1 && strcmp(_type, "DataArray_t") == 0) v = getArrayR8Skel(node);
    else v = getArrayR8(node);
  }
  else if (strcmp(_dtype, "C1") == 0)
  {
    IMPORTNUMPY;
    char* s = getArrayC1(node);
    E_Int l = strlen(s); npy_dim_vals2[0] = l;
    v = (PyObject*)PyArray_EMPTY(1, npy_dim_vals2, NPY_CHAR, 1);
    memcpy(PyArray_DATA((PyArrayObject*)v), s, l*sizeof(char));
    free(s);
  }
  else if (strcmp(_dtype, "MT") == 0)
  {
    v = Py_None; Py_INCREF(Py_None);
  }
  else
  {
    v = Py_None; Py_INCREF(Py_None);
    printf("Warning: adfcgnsread: unknown type of node: %s (%s).\n", _dtype,
           _name);
  }

  PyObject* s = Py_BuildValue("[sOOs]", _name, v, Py_None, _type);
  Py_DECREF(v);
  return s;
}

//=============================================================================
int K_IO::GenIOAdf::getSingleI4(double node)
{
  int r;
  ADF_Read_All_Data(node, (char*)(&r), &_errorFlag);
  return r;
}

//=============================================================================
E_LONG K_IO::GenIOAdf::getSingleI8(double node)
{
  E_LONG r;
  ADF_Read_All_Data(node, (char*)(&r), &_errorFlag);
  return r;
}

//=============================================================================
float K_IO::GenIOAdf::getSingleR4(double node)
{
  float r;
  ADF_Read_All_Data(node, (char*)(&r), &_errorFlag);
  return r;
}

//=============================================================================
double K_IO::GenIOAdf::getSingleR8(double node)
{
  double r;
  ADF_Read_All_Data(node, (char*)(&r), &_errorFlag);
  return r;
}

//=============================================================================
PyObject* K_IO::GenIOAdf::getArrayR8Skel(double node)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s, dim, sizem;
  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*_dims[s];
  if (sizem < _maxFloatSize) return getArrayR8(node);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
PyObject* K_IO::GenIOAdf::getArrayR8(double node)
{
  IMPORTNUMPY;
  int  dim;
  PyObject* r = NULL;

  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);

  // Create numpy: toujours en double
  npy_intp npy_dim_vals[dim];
  for (E_Int nn = 0; nn < dim; nn++) npy_dim_vals[nn] = _dims[nn];
  r = (PyObject*)PyArray_EMPTY(dim, npy_dim_vals, NPY_DOUBLE, 1);
  ADF_Read_All_Data(node, (char*)PyArray_DATA((PyArrayObject*)r), &_errorFlag);

  return r;
}

//=============================================================================
PyObject* K_IO::GenIOAdf::getArrayR4Skel(double node)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s, dim, sizem;
  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*_dims[s];
  if (sizem < _maxFloatSize) return getArrayR4(node);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
PyObject* K_IO::GenIOAdf::getArrayR4(double node)
{
  IMPORTNUMPY;
  int  s, dim, sizem;
  PyArrayObject* r = NULL;

  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*_dims[s];

  float* ptr = (float*)::malloc(sizem*sizeof(float));
  if (!ptr) { Py_INCREF(Py_None); return Py_None; }

  ADF_Read_All_Data(node, (char*)ptr, &_errorFlag);

  double* ptr2 = (double*)::malloc(sizem*sizeof(double));
  for (int n = 0; n < sizem; n++)
  {
    ptr2[n] = (double)(ptr[n]);
  }
  free(ptr);

  // Create numpy : toujours en doubles
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++)
  {
    npy_dim_vals[nn] = _dims[nn];
  }
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_DOUBLE, 1);
  memcpy(PyArray_DATA(r), ptr2, sizem*sizeof(double));
  free(ptr2);
  return (PyObject*)r;
}

//=============================================================================
PyObject* K_IO::GenIOAdf::getArrayI4Skel(double node)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s, dim, sizem;
  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*_dims[s];
  if (sizem < _maxFloatSize) return getArrayI4(node);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
PyObject* K_IO::GenIOAdf::getArrayI4(double node)
{
  IMPORTNUMPY;
  int  dim;
  PyArrayObject* r = NULL;

  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);

  // Create numpy : toujours en INT
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++) npy_dim_vals[nn] = _dims[nn];
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_INT, 1);

  ADF_Read_All_Data(node, (char*)PyArray_DATA(r), &_errorFlag);
  return (PyObject*)r;
}

//=============================================================================
PyObject* K_IO::GenIOAdf::getArrayI8Skel(double node)
{
  if (_maxFloatSize == 0) { Py_INCREF(Py_None); return Py_None; }
  int  s, dim, sizem;
  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*_dims[s];
  if (sizem < _maxFloatSize) return getArrayI8(node);
  else { Py_INCREF(Py_None); return Py_None; }
}

//=============================================================================
PyObject* K_IO::GenIOAdf::getArrayI8(double node)
{
  IMPORTNUMPY;
  int  s, dim, sizem;
  PyArrayObject* r = NULL;

  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);
  sizem = 1;
  for (s = 0; s < dim; s++) sizem = sizem*_dims[s];

  E_LONG* ptr = (E_LONG*)::malloc(sizem*sizeof(E_LONG));
  if (!ptr) { Py_INCREF(Py_None); return Py_None; }

  ADF_Read_All_Data(node, (char*)ptr, &_errorFlag);

  int* ptr2 = (int*)::malloc(sizem*sizeof(int));
  for (int n = 0; n < sizem; n++)
  {
    ptr2[n] = (int)(ptr[n]);
  }
  free(ptr);

  // Create numpy : toujours en int32
  vector<npy_intp> npy_dim_vals(dim);
  for (E_Int nn = 0; nn < dim; nn++)
  {
    npy_dim_vals[nn] = _dims[nn];
  }
  r = (PyArrayObject*)PyArray_EMPTY(dim, &npy_dim_vals[0], NPY_INT, 1);
  memcpy(PyArray_DATA(r), ptr2, sizem*sizeof(int));
  free(ptr2);
  return (PyObject*)r;
}

//=============================================================================
char* K_IO::GenIOAdf::getArrayC1(double node)
{
  int  dim, sizem, st;
  char *s;

  ADF_Get_Number_of_Dimensions(node, &dim, &_errorFlag);
  ADF_Get_Dimension_Values(node, _dims, &_errorFlag);

  sizem = 1;
  for (st = 0; st < dim; st++) sizem = sizem*_dims[st];
  //stringsize = _dims[0];
  s = (char*)::malloc((sizem*sizeof(char))+1); // size of char !
  ADF_Read_All_Data(node, (char*)s, &_errorFlag);
  s[(sizem*sizeof(char))] = '\0'; // C arrays starts at zero
  return s;
}
//=============================================================================
/*
   adfcgnswrite
*/
//=============================================================================
E_Int K_IO::GenIO::adfcgnswrite(char* file, PyObject* tree)
{
  if (tree == Py_None)
  {
    // nothing to write
    return 1;
  }

  /* Ouverture du fichier pour l'ecriture */
  FILE* f = fopen(file, "r");
  if (f != NULL)
  {
    fclose(f);
    remove(file);
  }
  double rootId;
  int errorFlag;
  ADF_Database_Open(file, "NEW", "NATIVE", &rootId, &errorFlag);

  PyObject* o;
  GenIOAdf ADF;

  int listsize = PyList_Size(tree);
  for (int n = 0; n < listsize; n++) // pour chaque Base
  {
    o = PyList_GetItem(tree, n);
    ADF._fatherStack.push_front(rootId);
    ADF.dumpOne(o);
    ADF._fatherStack.pop_front();
  }

  ADF_Database_Close(rootId, &errorFlag);
  return 0;
}

//=============================================================================
PyObject* K_IO::GenIOAdf::dumpOne(PyObject* tree)
{
  // ecrit le noeud courant
  double node = _fatherStack.front();
  node = writeNode(node, tree);

  // Dump les enfants
  if (PyList_Check(tree) == true && PyList_Size(tree) > 3)
  {
    PyObject* children = PyList_GetItem(tree, 2);
    if (PyList_Check(children) == true)
    {
      int nChildren = PyList_Size(children);
      for (E_Int i = 0; i < nChildren; i++)
      {
        _fatherStack.push_front(node);
        dumpOne(PyList_GetItem(children, i));
        _fatherStack.pop_front();
      }
    }
  }
  return tree;
}

//=============================================================================
double K_IO::GenIOAdf::writeNode(double node, PyObject* tree)
{
  IMPORTNUMPY;
  char s1[CGNSMAXLABEL+1];
  char s2[CGNSMAXLABEL+1];
  PyObject* pname = PyList_GetItem(tree, 0);
  char* name = PyString_AsString(pname);
  PyObject* plabel = PyList_GetItem(tree, 3);
  char* label = PyString_AsString(plabel);
  blankAndCopy(s1, name, CGNSMAXLABEL);
  blankAndCopy(s2, label, CGNSMAXLABEL);

  double child;
  // Creation du noeud
  ADF_Create(node, s1, &child, &_errorFlag);
  ADF_Set_Label(child, s2, &_errorFlag);

  // Ecriture de la valeur
  PyObject* v = PyList_GetItem(tree, 1);

  if (v == Py_None)
  {
    setArrayMT(child);
  }
  else if (PyString_Check(v) == true)
  {
    setArrayC1(child, PyString_AsString(v));
  }
  else if (PyInt_Check(v) == true)
  {
    setSingleI4(child, PyInt_AsLong(v));
  }
  else if (PyFloat_Check(v) == true)
  {
    if (strcmp(name, "CGNSLibraryVersion") == 0)
      setSingleR4(child, PyFloat_AsDouble(v));
    else
      setSingleR8(child, PyFloat_AsDouble(v));
  }
  else if (PyArray_Check(v) == true)
  {
    PyArrayObject* ar = (PyArrayObject*)v;
    int dim = PyArray_NDIM(ar);
    int* dims = new int [dim];
    //int typeNum = ar->descr->type_num;
    //int elSize = ar->descr->elsize;
    int typeNum = PyArray_TYPE(ar);
    int elSize = PyArray_ITEMSIZE(ar);
    for (int n = 0; n < dim; n++)
    {
      dims[n] = PyArray_DIMS(ar)[n];
    }

    // En python CGNS, tout est en tableau
    if (dim == 1 && dims[0] == 1) // valeur simple
    {
      if (typeNum == NPY_DOUBLE)
      {
        if (strcmp(name, "CGNSLibraryVersion") == 0)
        {
          double* ptr = (double*)PyArray_DATA(ar);
          setSingleR4(child, (float)ptr[0]);
        }
        else
        {
          double* ptr = (double*)PyArray_DATA(ar);
          setSingleR8(child, ptr[0]);
        }
      }
      else if (typeNum == NPY_INT)
      {
	if (elSize == 4)
	{
	  int* ptr = (int*)PyArray_DATA(ar);
	  setSingleI4(child, ptr[0]);
	}
        else
        {
          E_LONG* ptr = (E_LONG*)PyArray_DATA(ar);
          setSingleI8(child, ptr[0]);
        }
      }
      else if (typeNum == NPY_CHAR ||
               typeNum == NPY_STRING ||
               typeNum == NPY_BYTE ||
               //typeNum == NPY_SBYTE ||
               typeNum == NPY_UBYTE )
      {
	E_Int diml = PyArray_DIMS(ar)[0];
        char* buf = new char [diml+1];
        strncpy(buf, (char*)PyArray_DATA(ar), diml);
        buf[diml] = '\0';
        setArrayC1(child, buf);
        delete [] buf;
      }
      else if (typeNum == NPY_LONG)
      {
        if (elSize == 4)
	{
	  int* ptr = (int*)PyArray_DATA(ar);
	  setSingleI4(child, ptr[0]);
	}
        else
        {
          E_LONG* ptr = (E_LONG*)PyArray_DATA(ar);
          setSingleI8(child, ptr[0]);
        }
      }
      else if (typeNum == NPY_FLOAT)
      {
        float* ptr = (float*)PyArray_DATA(ar);
        setSingleR4(child, ptr[0]);
      }
    }
    else
    {
      if (typeNum == NPY_DOUBLE)
      {
        // patch pour la norme ADF
        if (strcmp(name, "RotationCenter") == 0 ||
            strcmp(name, "RotationAngle") == 0 ||
            strcmp(name, "RotationRateVector") == 0 ||
            strcmp(name, "Translation") == 0)
        {
          E_Int s = PyArray_Size(v);
          float* buf = new float [s];
          double* ptr = (double*)PyArray_DATA(ar);
          for (int i = 0; i < s; i++) buf[i] = ptr[i];
          setArrayR4(child, buf, dim, dims);
          delete [] buf;
        }
        else
          setArrayR8(child, (double*)PyArray_DATA(ar), dim, dims);
      }
      else if (typeNum == NPY_INT)
      {
        if (elSize == 4)
	{
	  setArrayI4(child, (int*)PyArray_DATA(ar), dim, dims);
	}
        else
        {
          setArrayI8(child, (E_LONG*)PyArray_DATA(ar), dim, dims);
        }
      }
      else if (typeNum == NPY_CHAR ||
               typeNum == NPY_STRING ||
               typeNum == NPY_BYTE ||
               //typeNum == NPY_SBYTE ||
               typeNum == NPY_UBYTE )
      {
	E_Int diml = PyArray_DIMS(ar)[0];
        char* buf = new char [diml+1];
        strncpy(buf, (char*)PyArray_DATA(ar), diml);
        buf[diml] = '\0';
        setArrayC1(child, buf);
        delete [] buf;
      }
      else if (typeNum == NPY_LONG)
      {
	if (elSize == 4)
	{
	  setArrayI4(child, (int*)PyArray_DATA(ar), dim, dims);
	}
        else
        {
          setArrayI8(child, (E_LONG*)PyArray_DATA(ar), dim, dims);
        }
      }
      else if (typeNum == NPY_FLOAT)
      {
        setArrayR4(child, (float*)PyArray_DATA(ar), dim, dims);
      }
    }
    delete [] dims;
  }

  return child;
}

//=============================================================================
void K_IO::GenIOAdf::blankAndCopy(char *s, const char *p, E_Int size)
{
  int i = 0;
  while (p[i] != '\0' && i < size) { s[i] = p[i]; i++; }
  for (E_Int j = i; j < size+1; j++) s[j] = ' ';
  if (p[i] != '\0')
  {
    printf("Warning: adfcgnswrite: %s node name has been truncated.\n", p);
    i = size-1;
  }
  s[i] = '\0';
}
//=============================================================================
double K_IO::GenIOAdf::setSingleR4(double node, float data)
{
  int dim;
  int dims[1];
  dim = 1; dims[0] = 1;
  ADF_Put_Dimension_Information(node, "R4", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, (char*)&data, &_errorFlag);

  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setSingleR8(double node, double data)
{
  int dim;
  int dims[1];

  dim = 1; dims[0] = 1;
  ADF_Put_Dimension_Information(node, "R8", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, (char*)&data, &_errorFlag);

  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setSingleI4(double node, int data)
{
  int dim;
  int dims[1];

  dim = 1; dims[0] = 1;
  ADF_Put_Dimension_Information(node, "I4", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, (char*)&data, &_errorFlag);

  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setSingleI8(double node, E_LONG data)
{
  int dim;
  int dims[1];

  dim = 1; dims[0] = 1;
  ADF_Put_Dimension_Information(node, "I8", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, (char*)&data, &_errorFlag);

  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setArrayMT(double node)
{
  int dim;
  dim = 1; _dims[0] = 0;
  ADF_Put_Dimension_Information(node, "MT", dim, _dims, &_errorFlag);

  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setArrayR8(double node, double *data, int dim,
                                  int *dims)
{
  ADF_Put_Dimension_Information(node, "R8", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, (char*)data, &_errorFlag);
  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setArrayI8(double node, E_LONG *data, int dim,
                                  int *dims)
{
  ADF_Put_Dimension_Information(node, "I8", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, (char*)data, &_errorFlag);
  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setArrayR4(double node, float *data, int dim,
                                  int *dims)
{
  ADF_Put_Dimension_Information(node, "R4", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, (char*)data, &_errorFlag);
  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setArrayI4(double node, int *data, int dim, int *dims)
{
  ADF_Put_Dimension_Information(node, "I4", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, (char*)data, &_errorFlag);
  return node;
}

//=============================================================================
double K_IO::GenIOAdf::setArrayC1(double node, char *data)
{
  int dim;
  int dims[1];
  char *ptr;

  dim = 1;
  dims[0] = strlen(data); // no +1 for \0

  // use begin(i) to get i-th field
  ptr = (char*)data;

  ADF_Put_Dimension_Information(node, "C1", dim, dims, &_errorFlag);
  ADF_Write_All_Data(node, ptr, &_errorFlag);

  return node;
}
//=============================================================================
// soit un chemin /A/B retourne A et B
//=============================================================================
void K_IO::GenIO::getABFromPath(char* path, char*& A, char*& B)
{
  E_Int i, j;
  E_Int l = strlen(path);
  A = new char [l+1];
  B = new char [l+1];
  i = 0;
  if (path[0] == '/') i = 1;
  j = 0;
  while (path[i] != '/' && path[i] != '\0')
  {
    A[j] = path[i]; i++; j++;
  }
  A[j] = '\0';
  j = 0;
  if (i < l-1) i++;
  while (path[i] != '/' && path[i] != '\0')
  {
    B[j] = path[i]; i++; j++;
  }
  B[j] = '\0';
  //printf("%s : %s %s\n", path, A, B);
}

//=============================================================================
// Lit les paths specifies dans le fichier file.
// Retourne une liste d'objets pythons contenant les noeuds pointes par les
// chemins
//=============================================================================
PyObject* K_IO::GenIO::adfcgnsReadFromPaths(char* file, PyObject* paths)
{
  if (PyList_Check(paths) == false)
  {
    PyErr_SetString(PyExc_TypeError,
                    "adfcgnsread: paths must be a list of strings.");
    return NULL;
  }
  E_Int size = PyList_Size(paths);
  for (E_Int i = 0; i < size; i++)
  {
    if (PyString_Check(PyList_GetItem(paths, i)) == false)
    {
      PyErr_SetString(PyExc_TypeError,
                      "adfcgnsread: paths must be a list of strings.");
      return NULL;
    }
  }

  /* Open file */
  double rootId; int errorFlag;
  ADF_Database_Open(file, "read_only", " ", &rootId, &errorFlag);
  if (errorFlag != -1)
  {
    PyErr_SetString(PyExc_TypeError, "adfcgnsread: cannot open file.");
    return NULL;
  }

  PyObject* ret = PyList_New(0);
  PyObject* node;
  GenIOAdf ADF;

  double id, id2;
  int nChildren; int lr;
  char label[CGNSMAXLABEL+1];

  for (E_Int i = 0; i < size; i++)
  {
    char* path = PyString_AsString(PyList_GetItem(paths, i));
    char* base; char* zone;
    int found = 0;
    getABFromPath(path, base, zone);
    // Get base node
    ADF_Number_of_Children(rootId, &nChildren, &errorFlag);
    for (int j = 1; j <= nChildren; j++)
    {
      ADF_Children_Names(rootId, j, 1, CGNSMAXLABEL+1,
                         &lr, (char *)label, &errorFlag);
      if (strcmp(label, base) == 0)
      {
        ADF_Get_Node_ID(rootId, label, &id, &errorFlag); found = 1; break;
      }
    }
    if (found != 0)
    {
      // Get zone node
      ADF_Number_of_Children(id, &nChildren, &errorFlag);
      for (int j = 1; j <= nChildren; j++)
      {
        ADF_Children_Names(id, j, 1, CGNSMAXLABEL+1,
                           &lr, (char*)label, &errorFlag);
        if (strcmp(label, zone) == 0)
        {
          ADF_Get_Node_ID(id, label, &id2, &errorFlag); found = 1; break;
        }
      }
    }

    delete [] base; delete [] zone;
    if (found == 0)
    { printf("Warning: adfcgnsread: cannot find this path.\n"); }
    else
    {
      node = ADF.createNode(id2);
      ADF._fatherStack.push_front(id2);
      ADF.loadOne(node,0);
      ADF._fatherStack.pop_front();
      ADF._fatherStack.clear();
      PyList_Append(ret, node); Py_DECREF(node);
    }
  }
  ADF_Database_Close(rootId, &errorFlag);
  return ret;
}
//=============================================================================
// Ecrit seulement les chemins specifies de l'arbre
//=============================================================================
E_Int K_IO::GenIO::adfcgnsWritePaths(char* file, PyObject* treeList,
                                     PyObject* paths)
{
  if (PyList_Check(paths) == false)
  {
    PyErr_SetString(PyExc_TypeError,
                    "adfcgnswrite: paths must be a list of strings.");
    return 1;
  }
  E_Int size = PyList_Size(paths);
  for (E_Int i = 0; i < size; i++)
  {
    if (PyString_Check(PyList_GetItem(paths, i)) == false)
    {
      PyErr_SetString(PyExc_TypeError,
                      "adwrite: paths must be a list of strings.");
      return 1;
    }
  }

  /* Ouverture du fichier pour l'ecriture */
  double rootId; int errorFlag;
  ADF_Database_Open(file, "old", " ", &rootId, &errorFlag);
  if (errorFlag != -1)
  {
    PyErr_SetString(PyExc_TypeError, "adcgnsfwrite: cannot open file.");
    return 1;
  }

  GenIOAdf ADF;
  int nChildren; int lr; int found;
  double id;
  char label[CGNSMAXLABEL+1];

  for (E_Int i = 0; i < size; i++)
  {
    char* path = PyString_AsString(PyList_GetItem(paths, i));
    char* base; char* zone;
    getABFromPath(path, base, zone);

    PyObject* node = PyList_GetItem(treeList, i);
    // Find base node
    found = 0;
    ADF_Number_of_Children(rootId, &nChildren, &errorFlag);
    for (int j = 1; j <= nChildren; j++)
    {
      ADF_Children_Names(rootId, j, 1, CGNSMAXLABEL+1,
                         &lr, (char*)label, &errorFlag);
      if (strcmp(label, base) == 0)
      {
        ADF_Get_Node_ID(rootId, label, &id, &errorFlag); found = 1; break;
      }
    }
    delete [] base; delete [] zone;

    if (found == 0)
      printf("Warning: adfcgnswrite: cannot write this path %s.\n", path);
    else
    {
      ADF._fatherStack.push_front(id);
      ADF.dumpOne(node);
      ADF._fatherStack.pop_front();
    }
  }

  ADF_Database_Close(rootId, &errorFlag);
  return 0;
}
