/*    
    Copyright 2013-2019 Onera.

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

#define K_ARRAY_UNIQUE_SYMBOL
#include "distributor2.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pydistributor2 [] =
{
  {"distribute", K_DISTRIBUTOR2::distribute, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
struct module_state {
    PyObject *error;
};
static int myextension_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}
static int myextension_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "distributor2",
        NULL,
        sizeof(struct module_state),
        Pydistributor2,
        NULL,
        myextension_traverse,
        myextension_clear,
        NULL
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
  #if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_distributor2();
  PyMODINIT_FUNC PyInit_distributor2()
#else
  PyMODINIT_FUNC initdistributor2();
  PyMODINIT_FUNC initdistributor2()
#endif
  {
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("distributor2", Pydistributor2);
#endif
    import_array();
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}
