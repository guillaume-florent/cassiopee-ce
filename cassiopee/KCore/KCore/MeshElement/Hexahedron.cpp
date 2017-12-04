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

#include "Hexahedron.h"

namespace K_MESH
{
  
const E_Int K_MESH::Hexahedron::NB_NODES=8;
const E_Int K_MESH::Hexahedron::NB_TRIS=12;

void Hexahedron::triangulate(E_Int* target)
{
  // WARNING: connectT3 is Apended (not cleared upon entry)
  
  K_MESH::Quadrangle q4;
  
  for (size_t f = 0; f < 6; ++f)
  {
    getBoundary(f, q4);
    K_MESH::Quadrangle::triangulate(q4.nodes(), target + f*6);
  }
}

}
