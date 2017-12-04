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

//=============================================================================
// Colormaps plugins
//=============================================================================
// Attention: si on utilise les SHADERS, les colormaps utilisees sont 
// dans createColormapTextures

// Parameters:
// IN: f -> scalar field value
// OUT: r, g, b -> red, green, blue color values.

//============================================================================
/*
  Blue to red colormap (par defaut pour les isos)
*/
//============================================================================
void colBlueToRed(Data* d, double f, float* r, float* g, float* b)
{
  if (f < 0.25) 
  {
    // bleu, vert augmente
    *r = 0.;
    *g = 4.*f;
    *b = 1.;
  }
  else if (f < 0.5)
  {
    // bleu diminue
    *r = 0.;
    *g = 1.;
    *b = 4.*(0.5-f);
  }
  else if (f < 0.75)
  {
    // rouge augmente
    *r = 4*(f-0.5);
    *g = 1.;
    *b = 0.;
  }
  else
  {
    // vert diminue, on finit rouge
    *r = 1.;
    *g = 4.*(1.-f);
    *b = 0.;
  }
  //printf("%f %f %f %f\n", f, *r, *g, *b);
}

//============================================================================
/*
  Green to red colormap (par defaut pour les maillages)
  Dans cette colormap, on evite le (0.1,0.1,1.) utilise pour la selection.
*/
//============================================================================
void colGreenToRed(Data* d, double f, float* r, float* g, float* b)
{
  if (f < 0.25)
  {
    // vert, bleu augmente
    *r = 0.;
    *g = 0.9;
    *b = 4.*f;
  }
  else if (f < 0.5)
  {
    // vert diminue (evite le bleu de selection)
    *r = 0.4;
    *g = 4.*(0.5-f);
    *b = 0.9;
  }
  else if (f < 0.75)
  {
    // rouge augmente
    *r = 4.*(f-0.5);
    *g = 0.9;
    *b = 4.*(0.75-f);
  }
  else
  {
    *r = 0.9;
    *g = 4.*(f-0.75);
    *b = 0.;
  }

  //printf("%f %f %f %f\n", f, *r, *g, *b);
}

//============================================================================
/*
  Grey colormap
*/
//============================================================================
void colGrey(Data* d, double f, float* r, float* g, float* b)
{
  *r = f;
  *g = f;
  *b = f;
}

//============================================================================
/*
  Grey colormap 2
*/
//============================================================================
void colGrey2(Data* d, double f, float* r, float* g, float* b)
{
  *r = 1.-f;
  *g = 1.-f;
  *b = 1.-f;
}


//=============================================================================
/*
  Diverging colormap
*/
//=============================================================================
void diverging(Data* d, double f, float* r, float* g, float* b)
{
  if (f < 0.5)
  {
    if (f < 0.25) 
    {
      if (f < 0.125)
      {
        if (f < 0.03125) 
        {*r = 0.23137254902+1.12941176471*f;
          *g = 0.298039215686+1.7568627451*f;
          *b = 0.752941176471+239.937254902*f; }
        else if (f < 0.0625) 
        {*r = 0.23137254902+1.12941176471*f;
          *g = 0.298039215686+1.7568627451*f;
          *b = 15.6588235294-237.050980392*f;  }
        else if (f < 0.09375) 
        {*r = 0.223529411765+1.25490196078*f;
          *g = 0.305882352941+1.63137254902*f;
          *b = 0.764705882353+1.25490196078*f; }
        else 
        {*r = 0.211764705882+1.38039215686*f;
          *g = 0.305882352941+1.63137254902*f;
          *b = 0.776470588235+1.12941176471*f;  }
      }
      else
      {   
        if (f < 0.15625) 
        { *r = 0.227450980392+1.25490196078*f;
          *g = 0.321568627451+1.50588235294*f;
          *b = 0.807843137255+0.878431372549*f; }
        else if (f < 0.1875) 
        { *r = 0.207843137255+1.38039215686*f;
          *g = 0.321568627451+1.50588235294*f;
          *b = 0.827450980392+0.752941176471*f; }
        else if (f < 0.21875) 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.345098039216+1.38039215686*f;
          *b = 0.874509803922+0.501960784314*f;  }
        else 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.345098039216+1.38039215686*f;
          *b = 0.901960784314+0.376470588235*f;  } 
      }
    }
    else
    { 
      if (f < 0.375)
      {
        if (f < 0.28125) 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.407843137255+1.12941176471*f;
          *b = 0.964705882353+0.125490196078*f;  }
        else if (f < 0.3125) 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.407843137255+1.12941176471*f;
          *b = 1.0+0.0*f;  }
        else if (f < 0.34375) 
        {*r = 0.207843137255+1.38039215686*f;
          *g = 0.486274509804+0.878431372549*f;
          *b = 1.07843137255-0.250980392157*f;  }
        else 
        {*r = 0.250980392157+1.25490196078*f;
          *g = 0.486274509804+0.878431372549*f;
          *b = 1.16470588235-0.501960784314*f;  }
      }
      else
      {
        if (f < 0.40625) 
        {*r = 0.250980392157+1.25490196078*f;
          *g = 0.580392156863+0.627450980392*f;
          *b = 1.21176470588-0.627450980392*f;  }
        else if (f < 0.4375) 
        {*r = 0.250980392157+1.25490196078*f;
          *g = 0.63137254902+0.501960784314*f;
          *b = 1.26274509804-0.752941176471*f;  }
        else if (f < 0.46875) 
        {*r = 0.305882352941+1.12941176471*f;
          *g = 0.741176470588+0.250980392157*f;
          *b = 1.37254901961-1.00392156863*f;  }
        else 
        {*r = 0.364705882353+1.00392156863*f;
          *g = 0.858823529412+0.0*f;
          *b = 1.43137254902-1.12941176471*f;  }
      }     
    }
  }
  else 
  {
    if (f < 0.75) 
    {
      if (f < 0.625)
      {
        if (f < 0.53125) 
        { *r = 0.364705882353+1.00392156863*f;
          *g = 0.733333333333+0.250980392157*f;
          *b = 1.61960784314-1.50588235294*f; }
        else if (f < 0.5625) 
        {*r = 0.43137254902+0.878431372549*f;
          *g = 1.2-0.627450980392*f;
          *b = 1.61960784314-1.50588235294*f;  }
        else if (f < 0.59375) 
        { *r = 0.572549019608+0.627450980392*f;
          *g = 1.2-0.627450980392*f;
          *b = 1.61960784314-1.50588235294*f; }
        else 
        {*r = 0.647058823529+0.501960784314*f;
          *g = 1.34901960784-0.878431372549*f;
          *b = 1.61960784314-1.50588235294*f;  }
      }
      else
      {
        if (f < 0.65625) 
        {*r = 0.803921568627+0.250980392157*f;
          *g = 1.42745098039-1.00392156863*f;
          *b = 1.69803921569-1.63137254902*f;  }
        else if (f < 0.6875) 
        { *r = 0.96862745098+0.0*f;
          *g = 1.50980392157-1.12941176471*f;
          *b = 1.61568627451-1.50588235294*f; }
        else if (f < 0.71875) 
        { *r = 0.96862745098+0.0*f;
          *g = 1.59607843137-1.25490196078*f;
          *b = 1.70196078431-1.63137254902*f; }
        else 
        {*r = 1.23921568627-0.376470588235*f;
          *g = 1.6862745098-1.38039215686*f;
          *b = 1.61176470588-1.50588235294*f;  }
      }
    }
    else
    {
      if (f < 0.875)
      {
        if (f < 0.78125) 
        { *r = 1.23921568627-0.376470588235*f;
          *g = 1.78039215686-1.50588235294*f;
          *b = 1.61176470588-1.50588235294*f; }
        else if (f < 0.8125) 
        {*r = 1.43529411765-0.627450980392*f;
          *g = 1.87843137255-1.63137254902*f;
          *b = 1.61176470588-1.50588235294*f; }
        else if (f < 0.84375) 
        { *r = 1.63921568627-0.878431372549*f;
          *g = 1.98039215686-1.7568627451*f;
          *b = 1.50980392157-1.38039215686*f; }
        else 
        {*r = 1.63921568627-0.878431372549*f;
          *g = 2.0862745098-1.88235294118*f;
          *b = 1.50980392157-1.38039215686*f;  }
      }
      else
      {
        if (f < 0.900625) 
        { *r = 1.85882352941-1.12941176471*f;
          *g = 2.19607843137-2.00784313725*f;
          *b = 1.50980392157-1.38039215686*f; }
        else if (f < 0.9375) 
        {*r = 1.97254901961-1.25490196078*f;
          *g = 4.2431372549-4.26666666667*f;
          *b = 1.39607843137-1.25490196078*f;  }
        else if (f < 0.96875) 
        {*r = 2.09019607843-1.38039215686*f;
          *g = 2.83137254902-2.76078431373*f;
          *b = 1.27843137255-1.12941176471*f;  }
        else 
        {*r = 2.21176470588-1.50588235294*f;
          *g = 4.53333333333-4.51764705882*f;
          *b = 1.27843137255-1.12941176471*f;  }
      }
    }
  }
}
