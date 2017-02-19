subroutine ZBUNK (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, &
     ALIM)
!
!! ZBUNK is subsidiary to ZBESH and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBUNI-A, ZBUNI-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU > FNUL.
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
!     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
!
!***SEE ALSO  ZBESH, ZBESK
!***ROUTINES CALLED  ZUNK1, ZUNK2
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZBUNK
!     COMPLEX Y,Z
  DOUBLE PRECISION ALIM, AX, AY, ELIM, FNU, TOL, YI, YR, ZI, ZR
  INTEGER KODE, MR, N, NZ
  DIMENSION YR(N), YI(N)
!***FIRST EXECUTABLE STATEMENT  ZBUNK
  NZ = 0
  AX = ABS(ZR)*1.7321D0
  AY = ABS(ZI)
  if (AY > AX) go to 10
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  call ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
  go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
  return
end
