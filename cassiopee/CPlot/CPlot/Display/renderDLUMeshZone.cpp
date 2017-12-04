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

#define PLOTNODE xi = x[i]; yi = y[i]; zi = z[i];               \
  dx = xi - xcam; dy = yi - ycam; dz = zi - zcam;               \
  dist = dx*dx + dy*dy + dz*dz;                                 \
  d = sqrt(dist)*dref;                                          \
  pru0 = d*(right[0] + up[0]);                                  \
  pru1 = d*(right[1] + up[1]);                                  \
  pru2 = d*(right[2] + up[2]);                                  \
  mru0 = d*(right[0] - up[0]);                                  \
  mru1 = d*(right[1] - up[1]);                                  \
  mru2 = d*(right[2] - up[2]);                                  \
  pt1[0] = xi - pru0;                                           \
  pt1[1] = yi - pru1;                                           \
  pt1[2] = zi - pru2;                                           \
  pt2[0] = xi + mru0;                                           \
  pt2[1] = yi + mru1;                                           \
  pt2[2] = zi + mru2;                                           \
  pt3[0] = xi + pru0;                                           \
  pt3[1] = yi + pru1;                                           \
  pt3[2] = zi + pru2;                                           \
  pt4[0] = xi - mru0;                                           \
  pt4[1] = yi - mru1;                                           \
  pt4[2] = zi - mru2;                                           \
  glVertex3dv(pt1); glVertex3dv(pt2);                           \
  glVertex3dv(pt3); glVertex3dv(pt4);

//=============================================================================
/*
  Display une zone en mesh avec une display list
  IN: zonep: pointeur sur la zone a afficher
  IN: zone: le no de zone non structure
  IN: zonet: le no de la zone dans liste globale des zones
  Le render des zones node est fait en direct.
*/
//=============================================================================
void DataDL::renderGPUUMeshZone(UnstructZone* zonep, int zone, int zonet)
{ 
  int i, ret;

  // Colormap
  float r, g, b;
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _plugins.colorMap->next->f;

  // Style
  float color1[3]; float color2[3];

  E_Float nz = 1./_numberOfUnstructZones;
#include "meshStyles.h"
      
#include "selection.h"

  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  //glCallList(zonep->_DLmesh);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glCallList(zImpl->_DLsolid);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // For node rendering (1D zones)
  double dref = 0.003;
  double xi, yi, zi;
  double viewMatrix[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, viewMatrix);
  double right[3];
  right[0] = viewMatrix[0];
  right[1] = viewMatrix[4];
  right[2] = viewMatrix[8];
  double up[3];
  up[0] = viewMatrix[1];
  up[1] = viewMatrix[5];
  up[2] = viewMatrix[9];
  double xcam = _view.xcam;
  double ycam = _view.ycam;
  double zcam = _view.zcam;
  double dx, dy, dz, dist, d;
  double pru0, pru1, pru2, mru0, mru1, mru2;
  double pt1[3]; double pt2[3]; double pt3[3]; double pt4[3];

  double* x = zonep->x;
  double* y = zonep->y;
  double* z = zonep->z;
  int eltType = zonep->eltType;

  // For BARS, NODE, 1D NGONS: display node
  if (eltType == 1 || eltType ==  0 || (eltType == 10 && zonep->nelts1D > 0)) 
  {
    glBegin(GL_QUADS);  
    if (zonep->blank == 0)
    {
      for (i = 0; i < zonep->np; i++)
      {
        PLOTNODE;
      }
    }
    else
    {
      for (i = 0; i < zonep->np; i++)
      {
        ret = _pref.blanking->f(this, i, zonep->blank, zone);
        if (ret != 0)
        {
          PLOTNODE;
        }
      }
    }
    glEnd();
  }
  glLineWidth(1.);
}