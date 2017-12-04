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
#define PLOT r = (f1[n1]-fmin1)*deltai;                                \
  g = (f2[n1]-fmin2)*deltai;                                           \
  b = (f3[n1]-fmin3)*deltai;                                           \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[n1s], surfy[n1s], surfz[n1s]);                       \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = (f1[n2]-fmin1)*deltai;                                           \
  g = (f2[n2]-fmin2)*deltai;                                           \
  b = (f3[n2]-fmin3)*deltai;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[n2s], surfy[n2s], surfz[n2s]);                       \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = (f1[n4]-fmin1)*deltai;                                           \
  g = (f2[n4]-fmin2)*deltai;                                           \
  b = (f3[n4]-fmin3)*deltai;                                           \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[n4s], surfy[n4s], surfz[n4s]);                       \
  glVertex3d(x[n4], y[n4], z[n4]);                                      \
  r = (f1[n3]-fmin1)*deltai;                                           \
  g = (f2[n3]-fmin2)*deltai;                                           \
  b = (f3[n3]-fmin3)*deltai;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[n3s], surfy[n3s], surfz[n3s]);                       \
  glVertex3d(x[n3], y[n3], z[n3]);                                      

#define PLOTB ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);   \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zone);               \
  ret4 = _pref.blanking->f(this, n4, zonep->blank, zone);               \
  if (ret1*ret2*ret3*ret4 != 0) { PLOT; }

  int i, j, k, n1, n2, n3, n4, n1s, n2s, n3s, n4s;
  float r, g, b;
  int ret1, ret2, ret3, ret4;

  // vector field
  double* f1 = zonep->f[nofield1];
  double* f2 = zonep->f[nofield2];
  double* f3 = zonep->f[nofield3];
  double fmin1, fmax1, fmin2, fmax2, fmin3, fmax3;
  fmax1 = maxf[nofield1]; fmin1 = minf[nofield1];
  fmax2 = maxf[nofield2]; fmin2 = minf[nofield2];
  fmax3 = maxf[nofield3]; fmin3 = minf[nofield3];
  double deltai1 = MAX(fmax1-fmin1, 1.e-6);
  double deltai2 = MAX(fmax2-fmin2, 1.e-6);
  double deltai3 = MAX(fmax3-fmin3, 1.e-6);
  double deltai = deltai1;
  deltai = MAX(deltai, deltai2);
  deltai = MAX(deltai, deltai3);
  deltai = 1./deltai;

  // Colormap
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _pref.colorMap->f;

  // Grid dimensions
  int ni = zonep->ni;
  int nj = zonep->nj;
  int nk = zonep->nk;
  if (ptrState->dim == 2) nk = 1;
  int nij = ni*nj;

  int nim = MIN(ni, 2);
  int njm = MIN(nj, 2);
  int nkm = MIN(nk, 2);
  int nbElti = nj*nk*nim;
  int nbEltj = ni*nk*njm;
  int nbEltk = nij*nkm;

  int nistepj = ni*stepj;
  int nijstepk = nij*stepk;
  int njstepk = nj*stepk;
  int nistepk = ni*stepk;
  int inci = (nim-1)*nj*nk;
  int incj = (njm-1)*ni*nk;
  int inck = (nkm-1)*nij;

  double* x = zonep->x; double* y = zonep->y; double* z = zonep->z;
  float* surfx = zonep->surf;
  float* surfy = surfx + nbElti;
  float* surfz = surfy + nbElti;

  glBegin(GL_QUADS);  

  // I Plane
  i = zonep->iPlane;
  if (i >= 0)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = i+j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = i+j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOTB;
        }
    }
  }
  else if (i == -1)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOT;
        }
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = inci+j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = ni-1+j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOTB;
        }
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = inci+j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = ni-1+j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOTB;
        }
    }
  }

  surfx = zonep->surf + 3*nbElti;
  surfy = surfx + nbEltj;
  surfz = surfy + nbEltj;

  // J Plane
  j = zonep->jPlane;
  if (j >= 0)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+j*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+j*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOTB;
        }
    }
  }
  else if (j == -1)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOT;
        }
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni+incj;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+(nj-1)*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOTB;
        }
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni+incj;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+(nj-1)*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOTB;
        }
    }
  }
  
  surfx = zonep->surf + 3*nbElti + 3*nbEltj;
  surfy = surfx + nbEltk;
  surfz = surfy + nbEltk;

  // K Plane
  k = zonep->kPlane;
  if (k >= 0)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOTB;
        }
    }
  }
  else if (k == -1)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOT;
        }
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni+inck;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni+(nk-1)*nij;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOTB;
        }
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni+inck;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni+(nk-1)*nij;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOTB;
        }
    }
  }
  glEnd();

  // Pour les lignes
  if (ni*nj == 1 || ni*nk == 1 || nj*nk == 1)
  {
    glLineWidth(3.);
    glPolygonOffset(-1.,-10.); // force offset
    glBegin(GL_LINES);
    int nie, nje, nke;
    nie = ni; nje = nj; nke = nk;
    if (ni*nj == 1) nke = nke-1;
    if (ni*nk == 1) nje = nje-1;
    if (nj*nk == 1) nie = nie-1;
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nke; k++)
        for (j = 0; j < nje; j++)
          for (i = 0; i < nie; i++)
          {
            n1 = i+j*ni+k*nij;
            n2 = n1+1;
            r = (f1[n1]-fmin1)*deltai;                                 
            g = (f2[n1]-fmin2)*deltai;                                 
            b = (f3[n1]-fmin3)*deltai;                                 
            glColor3f(r, g, b);             
            glVertex3d(x[n1], y[n1], z[n1]);
            r = (f1[n2]-fmin1)*deltai;                                 
            g = (f2[n2]-fmin2)*deltai;                                 
            b = (f3[n2]-fmin3)*deltai;
            glColor3f(r, g, b);
            glVertex3d(x[n2], y[n2], z[n2]);
          }
    }
    else
    {
      for (k = 0; k < nke; k++)
        for (j = 0; j < nje; j++)
          for (i = 0; i < nie; i++)
          {
            n1 = i+j*ni+k*nij;
            n2 = n1+1;
            ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);
            ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);
            if (ret1*ret2 != 0)
            { 
              r = (f1[n1]-fmin1)*deltai;                                 
              g = (f2[n1]-fmin2)*deltai;              
              b = (f3[n1]-fmin3)*deltai;     
              glVertex3d(x[n1], y[n1], z[n1]);
              r = (f1[n2]-fmin1)*deltai;                                 
              g = (f2[n2]-fmin2)*deltai;                                 
              b = (f3[n2]-fmin3)*deltai;
              glVertex3d(x[n2], y[n2], z[n2]);
            }
          }
    }
    glEnd();
    glLineWidth(1.);
  }
