C  
C    Copyright 2013-2017 Onera.
C
C    This file is part of Cassiopee.
C
C    Cassiopee is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    Cassiopee is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.

C =============================================================================
C Compute the volume of a 3D structured cell
C =============================================================================
      SUBROUTINE k6compvolofstructcell(ni, nj, nk, indcell,
     &                                 xt, yt, zt, vol)
C
      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk
      INTEGER_E indcell
      REAL_E xt(0:ni*nj*nk-1)
      REAL_E yt(0:ni*nj*nk-1)
      REAL_E zt(0:ni*nj*nk-1)
C_OUT
      REAL_E vol
C_LOCAL
      INTEGER_E ninj, i, j ,k
      INTEGER_E indA, indB, indC, indD, indE, indF, indG, indH
      REAL_E xA, xB, xC, xD, xE, xF, xG, xH
      REAL_E yA, yB, yC, yD, yE, yF, yG, yH
      REAL_E zA, zB, zC, zD, zE, zF, zG, zH
      REAL_E xcc, ycc, zcc
      REAL_E nx, ny, nz
      REAL_E v1, v2, v3, v4, v5, v6
      REAL_E const
C==============================================================================
C     1./24.
      const = 0.041666666667D0
C
      ninj = ni*nj
        
C     Warning : here i,j,k start from 0
      k = indcell/ninj
      j = indcell/ni - k*nj
      i = indcell - k*ninj - j*ni
      
C     Increments
      if (i .eq. ni-1) then 
         i = i-1
      endif
      if (j .eq. nj-1) then
         j = j-1
      endif
      if (k .eq. nk-1) then 
         k = k-1
      endif

      indA = i + j * ni + k * ninj
      indB = indA + 1           !(i+1,j,k)
      indC = indB + ni          !(i+1,j+1,k)
      indD = indA + ni          !(i,j+1,k)
      
      indE = indA + ninj        ! (i,j,k+1) 
      indF = indB + ninj        ! (i+1,j,k+1)
      indG = indC + ninj        ! (i+1,j+1,k+1)
      indH = indD + ninj        ! (i,j+1,k+1) 

      xA = xt(indA)
      yA = yt(indA)
      zA = zt(indA)
      
      xB = xt(indB)
      yB = yt(indB)
      zB = zt(indB)

      xC = xt(indC)
      yC = yt(indC)
      zC = zt(indC)

      xD = xt(indD)
      yD = yt(indD)
      zD = zt(indD)

      xE = xt(indE)
      yE = yt(indE)
      zE = zt(indE)

      xF = xt(indF)
      yF = yt(indF)
      zF = zt(indF)
      
      xG = xt(indG)
      yG = yt(indG)
      zG = zt(indG)

      xH = xt(indH)
      yH = yt(indH)
      zH = zt(indH)

C     Face ABCD
      xcc = xA + xB + xC + xD
      ycc = yA + yB + yC + yD
      zcc = zA + zB + zC + zD
      call k6crossproduct(xA, yA, zA, xB, yB, zB, xC, yC, zC, 
     &                    xD, yD, zD, nx, ny, nz)
      v1 = xcc * nx + ycc * ny + zcc * nz

C     Face EFGH
      xcc = xE + xF + xG + xH
      ycc = yE + yF + yG + yH
      zcc = zE + zF + zG + zH
      
      call k6crossproduct(xE, yE, zE, xF, yF, zF, xG, yG, zG, 
     &                    xH, yH, zH, nx, ny, nz) 
      v2 = xcc * nx + ycc * ny + zcc * nz

C     Face EADH
      xcc = xA + xD + xH + xE
      ycc = yA + yD + yH + yE
      zcc = zA + zD + zH + zE
      call k6crossproduct(xE, yE, zE, xA, yA, zA, xD, yD, zD, 
     &                    xH, yH, zH, nx, ny, nz) 
      v3 = xcc * nx + ycc * ny + zcc * nz
  
C     Face FBCG
      xcc = xB + xF + xG + xC
      ycc = yB + yF + yG + yC
      zcc = zB + zF + zG + zC
      call k6crossproduct(xF, yF, zF, xB, yB, zB, xC, yC, zC,
     &                    xG, yG, zG, nx, ny, nz)
      v4 = xcc * nx + ycc * ny + zcc * nz

C     Face EFBA
      xcc = xA + xB + xF + xE
      ycc = yA + yB + yF + yE
      zcc = zA + zB + zF + zE
      call k6crossproduct(xE, yE, zE, xF, yF, zF, xB, yB, zB,  
     &                    xA, yA, zA, nx, ny, nz)
      v5 = xcc * nx + ycc * ny + zcc * nz

C     Face HGCD
      xcc = xD + xC + xG + xH
      ycc = yD + yC + yG + yH
      zcc = zD + zC + zG + zH
      call k6crossproduct(xH, yH, zH, xG, yG, zG, xC, yC, zC, 
     &                    xD, yD, zD, nx, ny, nz)
      v6 = xcc * nx + ycc * ny + zcc * nz

      vol = v2 - v1 + v4 - v3 + v6 - v5
      vol = ABS(vol) * const 
      END

C=========================================================================
C computes the cross product AC ^ BD 
C=========================================================================
      subroutine k6crossproduct(xA, yA, zA, xB, yB, zB, xC, yC, zC,
     &                          xD, yD, zD, nx, ny, nz )

C
      IMPLICIT NONE
C
C==========================================================================
C_IN
      REAL_E xA, yA, zA, xB, yB, zB, xC, yC, zC,xD, yD, zD
C_OUT
      REAL_E nx, ny, nz 
C_LOCAL
      REAL_E xAC, yAC, zAC, xBD, yBD, zBD

      xAC = xC - xA
      yAC = yC - yA
      zAC = zC - zA

      xBD = xD - xB
      yBD = yD - yB
      zBD = zD - zB

      nx = yAC * zBD - zAC * yBD
      ny = zAC * xBD - xAC * zBD
      nz = xAC * yBD - yAC * xBD

      END
C ======= KCore/Metric/CompVolOfCellF.for ================================
