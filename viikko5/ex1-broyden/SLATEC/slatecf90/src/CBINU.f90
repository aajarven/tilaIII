subroutine CBINU (Z, FNU, KODE, N, CY, NZ, RL, FNUL, TOL, ELIM, ALIM)
!
!! CBINU is subsidiary to CAIRY, CBESH, CBESI, CBESJ, CBESK and CBIRY.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBINU-A, ZBINU-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
!
!***SEE ALSO  CAIRY, CBESH, CBESI, CBESJ, CBESK, CBIRY
!***ROUTINES CALLED  CASYI, CBUNI, CMLRI, CSERI, CUOIK, CWRSK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CBINU
  COMPLEX CW, CY, CZERO, Z
  REAL ALIM, AZ, DFNU, ELIM, FNU, FNUL, RL, TOL
  INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
  DIMENSION CY(N), CW(2)
  DATA CZERO / (0.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CBINU
  NZ = 0
  AZ = ABS(Z)
  NN = N
  DFNU = FNU + (N-1)
  if (AZ <= 2.0E0) go to 10
  if (AZ*AZ*0.25E0 > DFNU+1.0E0) go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     POWER SERIES
!-----------------------------------------------------------------------
  call CSERI(Z, FNU, KODE, NN, CY, NW, TOL, ELIM, ALIM)
  INW = ABS(NW)
  NZ = NZ + INW
  NN = NN - INW
  if (NN == 0) RETURN
  if (NW >= 0) go to 120
  DFNU = FNU + (NN-1)
   20 CONTINUE
  if (AZ < RL) go to 40
  if (DFNU <= 1.0E0) go to 30
  if (AZ+AZ < DFNU*DFNU) go to 50
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z
!-----------------------------------------------------------------------
   30 CONTINUE
  call CASYI(Z, FNU, KODE, NN, CY, NW, RL, TOL, ELIM, ALIM)
  if (NW < 0) go to 130
  go to 120
   40 CONTINUE
  if (DFNU <= 1.0E0) go to 70
   50 CONTINUE
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
  call CUOIK(Z, FNU, KODE, 1, NN, CY, NW, TOL, ELIM, ALIM)
  if (NW < 0) go to 130
  NZ = NZ + NW
  NN = NN - NW
  if (NN == 0) RETURN
  DFNU = FNU+(NN-1)
  if (DFNU > FNUL) go to 110
  if (AZ > FNUL) go to 110
   60 CONTINUE
  if (AZ > RL) go to 80
   70 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES
!-----------------------------------------------------------------------
  call CMLRI(Z, FNU, KODE, NN, CY, NW, TOL)
  if ( NW < 0) go to 130
  go to 120
   80 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
!-----------------------------------------------------------------------
  call CUOIK(Z, FNU, KODE, 2, 2, CW, NW, TOL, ELIM, ALIM)
  if (NW >= 0) go to 100
  NZ = NN
  DO 90 I=1,NN
    CY(I) = CZERO
   90 CONTINUE
  return
  100 CONTINUE
  if (NW > 0) go to 130
  call CWRSK(Z, FNU, KODE, NN, CY, NW, CW, TOL, ELIM, ALIM)
  if (NW < 0) go to 130
  go to 120
  110 CONTINUE
!-----------------------------------------------------------------------
!     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
!-----------------------------------------------------------------------
  NUI = FNUL-DFNU + 1
  NUI = MAX(NUI,0)
  call CBUNI(Z, FNU, KODE, NN, CY, NW, NUI, NLAST, FNUL, TOL, ELIM, &
   ALIM)
  if (NW < 0) go to 130
  NZ = NZ + NW
  if (NLAST == 0) go to 120
  NN = NLAST
  go to 60
  120 CONTINUE
  return
  130 CONTINUE
  NZ = -1
  if ( NW == (-2)) NZ=-2
  return
end
