subroutine CS1S2 (ZR, S1, S2, NZ, ASCLE, ALIM, IUF)
!
!! CS1S2 is subsidiary to CAIRY and CBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CS1S2-A, ZS1S2-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     CS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
!     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
!     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
!     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
!     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
!     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
!     PRECISION ABOVE THE UNDERFLOW LIMIT.
!
!***SEE ALSO  CAIRY, CBESK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CS1S2
  COMPLEX CZERO, C1, S1, S1D, S2, ZR
  REAL AA, ALIM, ALN, ASCLE, AS1, AS2, XX
  INTEGER IUF, NZ
  DATA CZERO / (0.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CS1S2
  NZ = 0
  AS1 = ABS(S1)
  AS2 = ABS(S2)
  AA = REAL(S1)
  ALN = AIMAG(S1)
  if (AA == 0.0E0 .AND. ALN == 0.0E0) go to 10
  if (AS1 == 0.0E0) go to 10
  XX = REAL(ZR)
  ALN = -XX - XX + ALOG(AS1)
  S1D = S1
  S1 = CZERO
  AS1 = 0.0E0
  if (ALN < (-ALIM)) go to 10
  C1 = CLOG(S1D) - ZR - ZR
  S1 = CEXP(C1)
  AS1 = ABS(S1)
  IUF = IUF + 1
   10 CONTINUE
  AA = MAX(AS1,AS2)
  if (AA > ASCLE) RETURN
  S1 = CZERO
  S2 = CZERO
  NZ = 1
  IUF = 0
  return
end