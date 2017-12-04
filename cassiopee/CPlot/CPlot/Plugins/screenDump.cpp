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
#include "../Data.h"
#include <stdlib.h>
#include "gl2ps.h"

//=============================================================================
// Screen dump plugins
//=============================================================================

//=============================================================================
// Fait le rendu, dump et ecrit le fichier
//=============================================================================
void Data::exportFile()
{
  int stateHeader, stateInfo, stateMenu, stateBB;
  stateHeader = ptrState->header;
  stateInfo = ptrState->info;
  stateMenu = ptrState->menu;
  stateBB = ptrState->bb;
  ptrState->header = 0;
  ptrState->info = 0;
  ptrState->menu = 0;
  ptrState->bb = 0;
  dumpWindow();
  ptrState->header = stateHeader;
  ptrState->info = stateInfo;
  ptrState->menu = stateMenu;
  ptrState->bb = stateBB;
  if (ptrState->continuousExport == 0) { ptrState->shootScreen = 0; }
  pthread_cond_signal(&ptrState->unlocked_export); // signal end of export
}

//=============================================================================
// Display to an image using FBO
// Retourne un buffer contenant l'image RGB
// exportWidth et exportHeight doivent etre Pairs
// Si mode=1, ajoute depth au buffer
//=============================================================================
char* Data::export2Image(int exportWidth, int exportHeight, int mode) 
{
  printf("mode = %d\n", mode);

  // resolution
  GLuint fb, rb, db;
#ifdef __SHADERS__
  glGenFramebuffersEXT(1, &fb);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fb);
  // Create and attach a color buffer
  glGenRenderbuffersEXT(1, &rb);
  // We must bind color_rb before we call glRenderbufferStorageEXT
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rb);
  // The storage format is RGBA8
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA8, 
                           exportWidth, exportHeight);
  // Attach color buffer to FBO
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
                               GL_RENDERBUFFER_EXT, rb);
  
  glGenRenderbuffersEXT(1, &db);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, db);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, 
                           GL_DEPTH_COMPONENT24, exportWidth, exportHeight);
  
  // Attach depth buffer to FBO
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
                               GL_RENDERBUFFER_EXT, db);

  int status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  if (status == GL_FRAMEBUFFER_COMPLETE_EXT)
  {
    // SUCCESS
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fb);
  }
  else 
  { // FAIL
    printf("Warning: CPlot: unable to get the correct frame buffer.\n");
    exportWidth = _view.w; exportHeight = _view.h;
    exportWidth = (exportWidth/2)*2; exportHeight = (exportHeight/2)*2;
  }
#endif

  int s = exportWidth*exportHeight;
  char* buffer;
  if (mode == 0) buffer = (char*)malloc(s*3);
  else buffer = (char*)malloc(s*4); // depth at the end

  // Switch viewport
  int viewWSav = _view.w; int viewHSav = _view.h;
  _view.w = exportWidth; _view.h = exportHeight;
  glViewport(0, 0, (GLsizei) exportWidth, (GLsizei) exportHeight);
  _view.ratio = (double)_view.w/(double)_view.h;

  if (ptrState->stereo == 0) display();
  else displayAnaglyph();

  // Back viewport
  _view.w = viewWSav; _view.h = viewHSav;
  _view.ratio = (double)_view.w/(double)_view.h;
  glViewport(0, 0, (GLsizei) _view.w, (GLsizei) _view.h);

  if (ptrState->offscreen != 1)
  {
    // Read FBO contents from color buffer with glReadPixels
    glReadPixels(0, 0, exportWidth, exportHeight, 
                 GL_RGB, GL_UNSIGNED_BYTE, buffer);
    if (mode == 1) // get zbuffer
    {
      float* zbuf = new float [s];
      glReadPixels(0, 0, exportWidth, exportHeight, 
                   GL_DEPTH_COMPONENT, GL_FLOAT, zbuf);
      delete [] zbuf;
      for (int i = 0; i < s; i++) { printf("%f ", zbuf[i]); buffer[i] = int(zbuf[i]*255); }
      printf("\n");
    }
  }
  else 
  { // mesa offscreen rendering
    char* buffRGBA = (char*)ptrState->offscreenBuffer;
    for (int i = 0; i < s; i++)
    { 
      buffer[3*i] = buffRGBA[4*i];
      buffer[3*i+1] = buffRGBA[4*i+1];
      buffer[3*i+2] = buffRGBA[4*i+2];
    }
  }

#ifdef __SHADERS__
  // Delete FBO
  glDeleteRenderbuffersEXT(1, &rb);
  glDeleteRenderbuffersEXT(1, &db);
  //Bind 0, which means render to back buffer, as a result, fb is unbound
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  glDeleteFramebuffersEXT(1, &fb);
#endif

  return buffer;
}

//=============================================================================
/*
  Dump a window to a buffer suivant l'extension du plugins.
*/
//=============================================================================
void Data::dumpWindow() 
{
  char fileName[120];
  if (_pref.screenDump == NULL) return;

  // File name
  if (strcmp(ptrState->exportFile, "CPlot") == 0)
    sprintf(fileName, "%s/%s.%d.%s", ptrState->localPathName, 
            ptrState->exportFile, ptrState->exportNumber, 
            _pref.screenDump->extension);
  else strcpy(fileName, ptrState->exportFile);
  ptrState->exportNumber++;

  if (strcmp(_pref.screenDump->extension, "ps") == 0)
  {
    // Postscript
    int state = GL2PS_OVERFLOW, buffsize = 0;

    FILE* fp = fopen(fileName, "wb");
    while (state == GL2PS_OVERFLOW)
    {
      buffsize += 1024*1024;
      gl2psBeginPage("test", "CPlot", NULL, 
                     GL2PS_EPS, GL2PS_BSP_SORT, 
		     GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT, 
		     GL_RGBA, 0, NULL, 0, 0, 0,  buffsize, fp, fileName);
      display();
      state = gl2psEndPage();
    }
    printf("Wrote file %s.\n", fileName);
    fclose(fp);
  }
  else // Other formats
  {
    int mode = 0;
    if (strcmp(_pref.screenDump->extension, "dpng") == 0) mode = 1;

    int exportWidth = _view.w;
    int exportHeight = _view.h;
    double r = _view.h * 1. / _view.w;
    if (ptrState->exportWidth != -1 && ptrState->exportHeight != -1)
    {
      exportWidth = ptrState->exportWidth; 
      exportHeight = int(exportWidth * r);
    }
    else if (ptrState->exportWidth != -1 && ptrState->exportHeight == -1)
    {
      exportWidth = ptrState->exportWidth;
      exportHeight = int(exportWidth * r);
    }
    else if (ptrState->exportWidth == -1 && ptrState->exportHeight != -1)
    {
      exportHeight = ptrState->exportHeight;
      exportWidth = int(exportHeight * (1./r));
    }
    exportWidth = (exportWidth/2)*2;
    exportHeight = (exportHeight/2)*2; // doit etre pair
    
    char* buffer; char* buffer2;
    int antialiasing = 1;
    if (strcmp(_pref.screenDump->extension, "dpng") == 0) antialiasing = 0;

    if (antialiasing == 1)
    {
      // get image X2
      buffer = export2Image(2*exportWidth, 2*exportHeight, mode);
    
      // blur image X2
      buffer2 = (char*)malloc(3*2*exportWidth*2*exportHeight*sizeof(char));
      int nitBlur;
      if (ptrState->mode == MESH) nitBlur = int(exportWidth/500.)+1;
      else nitBlur = 1;
      gaussianBlur(2*exportWidth, 2*exportHeight, buffer, buffer2, nitBlur, 0.1);

      // supersample X2
      free(buffer);
      buffer = (char*)malloc(3*exportWidth*exportHeight*sizeof(char));
      superSample(exportWidth, exportHeight, buffer2, buffer, 2);
      free(buffer2);
    }
    else
    {
      buffer = export2Image(exportWidth, exportHeight, mode);
    }
    
    if (strcmp(_pref.screenDump->extension, "mpeg") == 0)
    {
      buffer2 = (char*)malloc(3*exportWidth*exportHeight*sizeof(char));
      gaussianBlur(exportWidth, exportHeight, buffer, buffer2, 1, 0.1);
      //sharpenImage(exportWidth, exportHeight, buffer, buffer2, 5.,
      //             2, 2);
      for (int i = 0; i < 3*exportWidth*exportHeight; i++) buffer[i] = buffer2[i];
      free(buffer2);
    }

    // Dump the buffer to a file
    _pref.screenDump->f(this, fileName, buffer, exportWidth, exportHeight, 0);    
    free(buffer);
  }
}

//=============================================================================
void Data::finalizeExport()
{
  int exportWidth = _view.w;
  int exportHeight = _view.h;
  if (ptrState->exportWidth != -1) exportWidth = ptrState->exportWidth;
  if (ptrState->exportHeight != -1) exportHeight = ptrState->exportHeight;
  exportWidth = (exportWidth/2)*2;
  exportHeight = (exportHeight/2)*2; // doit etre pair
  char* buffer = NULL;
  // Dump the buffer to a file
  _pref.screenDump->f(this, ptrState->exportFile, buffer, 
                      exportWidth, exportHeight, 1);
}
