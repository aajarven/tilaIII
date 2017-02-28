subroutine ZS1S2 (ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM, IUF)
!
!! ZS1S2 is subsidiary to ZAIRY and ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CS1S2-A, ZS1S2-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
!     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
!     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
!     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
!     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
!     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
!     PRECISION ABOVE THE UNDERFLOW LIMIT.
!
!***SEE ALSO  ZAIRY, ZBESK
!***ROUTINES CALLED  ZABS, ZEXP, ZLOG
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)
!***END PROLOGUE  ZS1S2
!     COMPLEX CZERO,C1,S1,S1D,S2,ZR
  DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI, &
   S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR, ZABS
  INTEGER IUF, IDUM, NZ
  EXTERNAL ZABS, ZEXP, ZLOG
  DATA ZEROR,ZEROI  / 0.0D0 , 0.0D0 /
!***FIRST EXECUTABLE STATEMENT  ZS1S2
  NZ = 0
  AS1 = ZABS(S1R,S1I)
  AS2 = ZABS(S2R,S2I)
  if (S1R == 0.0D0 .AND. S1I == 0.0D0) go to 10
  if (AS1 == 0.0D0) go to 10
  ALN = -ZRR - ZRR + LOG(AS1)
  S1DR = S1R
  S1DI = S1I
  S1R = ZEROR
  S1I = ZEROI
  AS1 = ZEROR
  if (ALN < (-ALIM)) go to 10
  call ZLOG(S1DR, S1DI, C1R, C1I, IDUM)
  C1R = C1R - ZRR - ZRR
  C1I = C1I - ZRI - ZRI
  call ZEXP(C1R, C1I, S1R, S1I)
  AS1 = ZABS(S1R,S1I)
  IUF = IUF + 1
   10 CONTINUE
  AA = MAX(AS1,AS2)
  if (AA > ASCLE) RETURN
  S1R = ZEROR
  S1I = ZEROI
  S2R = ZEROR
  S2I = ZEROI
  NZ = 1
  IUF = 0
  return
end