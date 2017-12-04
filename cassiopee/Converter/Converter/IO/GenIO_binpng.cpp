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

// Binary PNG (Portable Network Graphics) file support

# include <png.h>
# include "GenIO.h"
# include "Array/Array.h"
# include "Connect/connect.h"
# include <vector>
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   pngread 
   Un seul domaine peut-etre stocke dans ce type de fichier.
   Limite a un bit_depth de 8.
*/
//=============================================================================
E_Int K_IO::GenIO::pngread( 
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, 
  vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: pngread: cannot open file %s.\n", file);
    return 1;
  }

  int width, height;
  png_structp png_ptr;
  png_infop info_ptr, end_info_ptr;
  //int number_of_passes;
  png_bytep* row_pointers;
  png_byte color_type;
  png_byte bit_depth;

  // Header
  png_byte header[8];    // 8 is the maximum size that can be checked
  fread(header, 1, 8, ptrFile);
  if (png_sig_cmp(header, 0, 8))
  { printf("File %s is an invalid png file.\n", file); fclose(ptrFile); return 1; }
  
  /* initialize stuff */
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr)
  { printf("Can not create png structures.\n"); fclose(ptrFile); return 1; }
  
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  { printf("Can not create png structures.\n"); fclose(ptrFile); return 1; }
               
  end_info_ptr = png_create_info_struct(png_ptr);
  if (!end_info_ptr)
  { printf("Can not create png structures.\n"); fclose(ptrFile); return 1; }
 
  if (setjmp(png_jmpbuf(png_ptr)))
  { printf("Can not create png structures.\n"); fclose(ptrFile); return 1; }
              
  png_init_io(png_ptr, ptrFile);
  png_set_sig_bytes(png_ptr, 8);

  png_read_info(png_ptr, info_ptr);

  width = png_get_image_width(png_ptr, info_ptr);
  height = png_get_image_height(png_ptr, info_ptr);
  color_type = png_get_color_type(png_ptr, info_ptr);
  //printf("color type %d\n", color_type);
  if (color_type == 3)
  {
    printf("Palette indexed colors is not supported.\n"); 
    fclose(ptrFile); return 1;
  }

  bit_depth = png_get_bit_depth(png_ptr, info_ptr);
  //printf("bit depth %d\n", bit_depth);
  if (bit_depth != 8)
  {
    printf("Color depth must be 8.\n"); 
    fclose(ptrFile); return 1;
  }

  //number_of_passes = png_set_interlace_handling(png_ptr);
  png_read_update_info(png_ptr, info_ptr);

  varString = new char [128];

  /* read file */
  if (setjmp(png_jmpbuf(png_ptr)))
  { printf("Error during read of png image.\n"); fclose(ptrFile); return 1; }

  row_pointers = (png_bytep*)malloc(sizeof(png_bytep)*height);
  for (int y = 0; y < height; y++)
    row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png_ptr, info_ptr));

  png_read_image(png_ptr, row_pointers);

  // Stockage du champ
  E_Int nil = width;
  E_Int njl = height;
  FldArrayF* f;
  if (color_type == 0) // greyscale
  {
    strcpy(varString, "x,y,z,r");
    f = new FldArrayF(nil*njl, 4);
  }
  else if (color_type == 2) // RGB
  {
    strcpy(varString, "x,y,z,r,g,b");
    f = new FldArrayF(nil*njl, 6);
  }
  else if (color_type == 4) // greyscale + alpha
  {
    strcpy(varString, "x,y,z,r,a");
    f = new FldArrayF(nil*njl, 5);
  }
  else if (color_type == 6) // RGB + alpha
  {
    strcpy(varString, "x,y,z,r,g,b,a");
    f = new FldArrayF(nil*njl, 7);
  }
 
  f->setAllValuesAtNull();
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);

  structField.push_back(f);
  ni.push_back(nil); nj.push_back(njl); nk.push_back(1);
  
  // Cree les noms de zones
  char* zoneName = new char [128];
  sprintf(zoneName, "Zone0");
  zoneNames.push_back(zoneName);

  for (E_Int j = 0; j < njl; j++)
    for (E_Int i = 0; i < nil; i++)
    {
      fx[i+j*nil] = i*1.;
      fy[i+j*nil] = (njl-j)*1.;
      fz[i+j*nil] = 0.;
    }

  if (color_type == 0) // greyscale
  {
    E_Float* r = f->begin(4);
    for (E_Int j = 0; j < njl; j++)
    {
      png_byte* p = row_pointers[j];
      for (E_Int i = 0; i < nil; i++)
      {
        r[i+j*nil] = *p; p++;
      }
    }
  }
  else if (color_type == 2)
  {
    E_Float* r = f->begin(4);
    E_Float* g = f->begin(5);
    E_Float* b = f->begin(6);
    for (E_Int j = 0; j < njl; j++)
    {
      png_byte* p = row_pointers[j];
      for (E_Int i = 0; i < nil; i++)
      {
        r[i+j*nil] = *p; p++;
        g[i+j*nil] = *p; p++;
        b[i+j*nil] = *p; p++;
      }
    }
  }
  else if (color_type == 4)
  {
    E_Float* r = f->begin(4);
    E_Float* alpha = f->begin(5);
    for (E_Int j = 0; j < njl; j++)
    {
      png_byte* p = row_pointers[j];
      for (E_Int i = 0; i < nil; i++)
      {
        r[i+j*nil] = *p; p++;
        alpha[i+j*nil] = *p; p++;
      }
    }
  }
  else if (color_type == 6)
  {
    E_Float* r = f->begin(4);
    E_Float* g = f->begin(5);
    E_Float* b = f->begin(6);
    E_Float* alpha = f->begin(7);
    for (E_Int j = 0; j < njl; j++)
    {
      png_byte* p = row_pointers[j];
      for (E_Int i = 0; i < nil; i++)
      {
        r[i+j*nil] = *p; p++;
        g[i+j*nil] = *p; p++;
        b[i+j*nil] = *p; p++;
        alpha[i+j*nil] = *p; p++;
      }
    }
  }

  for (int y = 0; y < height; y++) free(row_pointers[y]);
  free(row_pointers);
  
  png_read_end(png_ptr, end_info_ptr);
  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info_ptr);

  fclose(ptrFile);
  return 0;
}

//=============================================================================
// This is a non-sense function
//=============================================================================
E_Int K_IO::GenIO::pngwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  printf("Warning: pngwrite: Can not be used for output.");
  return 1;
}
