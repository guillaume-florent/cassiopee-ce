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
# include "cplot.h"
# include "Data.h"

//=============================================================================
/* Arrete l'export continu, finalize les MPEG, fait une barriere pour
   attendre la fin de l'ecriture */
//=============================================================================
PyObject* K_CPLOT::finalizeExport(PyObject* self, PyObject* args)
{
  int finalizeType = 0;
  if (!PyArg_ParseTuple(args, "i", &finalizeType)) return NULL;

  Data* d = Data::getInstance();
  // Bloc en attendant la fin de l'ecriture
  if (d->ptrState->continuousExport == 0)
  {
    pthread_mutex_lock(&d->ptrState->export_mutex);
    if (d->ptrState->shootScreen == 1)
      pthread_cond_wait (&d->ptrState->unlocked_export, &d->ptrState->export_mutex); 
    pthread_mutex_unlock(&d->ptrState->export_mutex);   
  }
  // Finalize mpeg
  if (finalizeType == 1 && strcmp(d->_pref.screenDump->extension, "mpeg") == 0)
    d->finalizeExport(); // force l'ecriture finale du fichier
  
  d->ptrState->continuousExport = 0;
  d->ptrState->shootScreen = 0;
  return Py_BuildValue("l", KSUCCESS);
}
