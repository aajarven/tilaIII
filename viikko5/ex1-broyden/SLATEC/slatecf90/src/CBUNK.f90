subroutine CBUNK (Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
!
!! CBUNK is subsidiary to CBESH and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CBUNK-A, ZBUNK-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU > FNUL.
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
!     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2
!
!***SEE ALSO  CBESH, CBESK
!***ROUTINES CALLED  CUNK1, CUNK2
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CBUNK
  COMPLEX Y, Z
  REAL ALIM, AX, AY, ELIM, FNU, TOL, XX, YY
  INTEGER KODE, MR, N, NZ
  DIMENSION Y(N)
!***FIRST EXECUTABLE STATEMENT  CBUNK
  NZ = 0
  XX = REAL(Z)
  YY = AIMAG(Z)
  AX = ABS(XX)*1.7321E0
  AY = ABS(YY)
  if (AY > AX) go to 10
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  call CUNK1(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
  go to 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
  call CUNK2(Z, FNU, KODE, MR, N, Y, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
  return
end
