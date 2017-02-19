subroutine ZBINU (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, &
     TOL, ELIM, ALIM)
!
!! ZBINU is subsidiary to ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK and ZBIRY.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBINU-A, ZBINU-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
!
!***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBIRY
!***ROUTINES CALLED  ZABS, ZASYI, ZBUNI, ZMLRI, ZSERI, ZUOIK, ZWRSK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZBINU
  DOUBLE PRECISION ALIM, AZ, CWI, CWR, CYI, CYR, DFNU, ELIM, FNU, &
   FNUL, RL, TOL, ZEROI, ZEROR, ZI, ZR, ZABS
  INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
  DIMENSION CYR(N), CYI(N), CWR(2), CWI(2)
  EXTERNAL ZABS
  DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
!***FIRST EXECUTABLE STATEMENT  ZBINU
  NZ = 0
  AZ = ZABS(ZR,ZI)
  NN = N
  DFNU = FNU + (N-1)
  if (AZ <= 2.0D0) go to 10
  if (AZ*AZ*0.25D0 > DFNU+1.0D0) go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     POWER SERIES
!-----------------------------------------------------------------------
  call ZSERI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
  INW = ABS(NW)
  NZ = NZ + INW
  NN = NN - INW
  if (NN == 0) RETURN
  if (NW >= 0) go to 120
  DFNU = FNU + (NN-1)
   20 CONTINUE
  if (AZ < RL) go to 40
  if (DFNU <= 1.0D0) go to 30
  if (AZ+AZ < DFNU*DFNU) go to 50
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z
!-----------------------------------------------------------------------
   30 CONTINUE
  call ZASYI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM, &
   ALIM)
  if (NW < 0) go to 130
  go to 120
   40 CONTINUE
  if (DFNU <= 1.0D0) go to 70
   50 CONTINUE
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
  call ZUOIK(ZR, ZI, FNU, KODE, 1, NN, CYR, CYI, NW, TOL, ELIM, &
   ALIM)
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
  call ZMLRI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL)
  if ( NW < 0) go to 130
  go to 120
   80 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
!-----------------------------------------------------------------------
  call ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM, &
   ALIM)
  if (NW >= 0) go to 100
  NZ = NN
  DO 90 I=1,NN
    CYR(I) = ZEROR
    CYI(I) = ZEROI
   90 CONTINUE
  return
  100 CONTINUE
  if (NW > 0) go to 130
  call ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL, &
   ELIM, ALIM)
  if (NW < 0) go to 130
  go to 120
  110 CONTINUE
!-----------------------------------------------------------------------
!     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
!-----------------------------------------------------------------------
  NUI = FNUL-DFNU + 1
  NUI = MAX(NUI,0)
  call ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL, &
   TOL, ELIM, ALIM)
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
