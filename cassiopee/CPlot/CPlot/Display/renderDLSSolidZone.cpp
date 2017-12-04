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
#include "../DataDL.h"
#include "../ZoneImplDL.h"

//=============================================================================
/*
  Display une zone en solid ou en material.
*/
//=============================================================================
void DataDL::renderGPUSSolidZone(StructZone* zonep, int zone)
{
  // Style
  float color1[3]; float color2[3];

  // Colormap
  float r, g, b;
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _plugins.colorMap->next->f;

  E_Float nz = 1./_numberOfStructZones;
#include "solidStyles.h"

  // Ecrasement si renderTag
  if (zonep->colorR > -0.5)
  {color1[0] = zonep->colorR; 
    color1[1] = zonep->colorG; 
    color1[2] = zonep->colorB;}

#include "selection.h"

  bool is1D = ((zonep->ni*zonep->nj == 1) | (zonep->ni*zonep->nk == 1) | (zonep->nj*zonep->nk == 1));
  if (is1D == true && ptrState->mode == RENDER) glLineWidth(1.+5*zonep->shaderParam1);
  else if (is1D == true) glLineWidth(3.);
  else glLineWidth(1.);

  // scale
  E_Float s = MAX(zonep->xmax-zonep->xmin, zonep->ymax-zonep->ymin);
  s = MAX(s, zonep->zmax-zonep->zmin);
  s = 100./(s+1.e-12);

#ifdef __SHADERS__
  if (ptrState->mode == RENDER)
  {
    if (zonep->selected == 1 && zonep->active == 1) 
      triggerShader(*zonep, zonep->material, s, color2);
    else triggerShader(*zonep, zonep->material, s, color1);
  }
  else
  {
    if (zonep->selected == 1 && zonep->active == 1) 
      triggerShader(*zonep, 0, s, color2);
    else triggerShader(*zonep, 0, s, color1);
  }
  //if (is1D == true) { light(2); glColor3f(1,0,0); printf("deactivating\n"); _shaders.deactivate(); }
#endif

  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  glCallList(zImpl->_DLsolid);
  glLineWidth(1.);
}
