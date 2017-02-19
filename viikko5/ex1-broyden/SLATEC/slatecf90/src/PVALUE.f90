subroutine PVALUE (L, NDER, X, YFIT, YP, A)
!
!! PVALUE uses the coefficients generated by POLFIT to evaluate the ...
!            polynomial fit of degree L, along with the first NDER of ...
!            its derivatives, at a specified point.
!
!***LIBRARY   SLATEC
!***CATEGORY  K6
!***TYPE      SINGLE PRECISION (PVALUE-S, DP1VLU-D)
!***KEYWORDS  CURVE FITTING, LEAST SQUARES, POLYNOMIAL APPROXIMATION
!***AUTHOR  Shampine, L. F., (SNLA)
!           Davenport, S. M., (SNLA)
!***DESCRIPTION
!
!     Written by L. F. Shampine and S. M. Davenport.
!
!     Abstract
!
!     The subroutine  PVALUE  uses the coefficients generated by  POLFIT
!     to evaluate the polynomial fit of degree  L , along with the first
!     NDER  of its derivatives, at a specified point.  Computationally
!     stable recurrence relations are used to perform this task.
!
!     The parameters for  PVALUE  are
!
!     Input --
!         L -      the degree of polynomial to be evaluated.  L  may be
!                  any non-negative integer which is less than or equal
!                  to  NDEG , the highest degree polynomial provided
!                  by  POLFIT .
!         NDER -   the number of derivatives to be evaluated.  NDER
!                  may be 0 or any positive value.  If NDER is less
!                  than 0, it will be treated as 0.
!         X -      the argument at which the polynomial and its
!                  derivatives are to be evaluated.
!         A -      work and output array containing values from last
!                  call to  POLFIT .
!
!     Output --
!         YFIT -   value of the fitting polynomial of degree  L  at  X
!         YP -     array containing the first through  NDER  derivatives
!                  of the polynomial of degree  L .  YP  must be
!                  dimensioned at least  NDER  in the calling program.
!
!***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!                 Curve fitting by polynomials in one variable, Report
!                 SLA-74-0270, Sandia Laboratories, June 1974.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   740601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  PVALUE
  DIMENSION YP(*),A(*)
  CHARACTER*8 XERN1, XERN2
!***FIRST EXECUTABLE STATEMENT  PVALUE
  if (L  <  0) go to 12
  NDO = MAX(NDER,0)
  NDO = MIN(NDO,L)
  MAXORD = A(1) + 0.5
  K1 = MAXORD + 1
  K2 = K1 + MAXORD
  K3 = K2 + MAXORD + 2
  NORD = A(K3) + 0.5
  if (L  >  NORD) go to 11
  K4 = K3 + L + 1
  if (NDER  <  1) go to 2
  DO 1 I = 1,NDER
 1      YP(I) = 0.0
 2    if (L  >=  2) go to 4
  if (L  ==  1) go to 3
!
! L IS 0
!
  VAL = A(K2+1)
  go to 10
!
! L IS 1
!
 3    CC = A(K2+2)
  VAL = A(K2+1) + (X-A(2))*CC
  if (NDER  >=  1) YP(1) = CC
  go to 10
!
! L IS GREATER THAN 1
!
 4    NDP1 = NDO + 1
  K3P1 = K3 + 1
  K4P1 = K4 + 1
  LP1 = L + 1
  LM1 = L - 1
  ILO = K3 + 3
  IUP = K4 + NDP1
  DO 5 I = ILO,IUP
 5      A(I) = 0.0
  DIF = X - A(LP1)
  KC = K2 + LP1
  A(K4P1) = A(KC)
  A(K3P1) = A(KC-1) + DIF*A(K4P1)
  A(K3+2) = A(K4P1)
!
! EVALUATE RECURRENCE RELATIONS FOR FUNCTION VALUE AND DERIVATIVES
!
  DO 9 I = 1,LM1
    IN = L - I
    INP1 = IN + 1
    K1I = K1 + INP1
    IC = K2 + IN
    DIF = X - A(INP1)
    VAL = A(IC) + DIF*A(K3P1) - A(K1I)*A(K4P1)
    if (NDO  <=  0) go to 8
    DO 6 N = 1,NDO
      K3PN = K3P1 + N
      K4PN = K4P1 + N
 6        YP(N) = DIF*A(K3PN) + N*A(K3PN-1) - A(K1I)*A(K4PN)
!
! SAVE VALUES NEEDED FOR NEXT EVALUATION OF RECURRENCE RELATIONS
!
    DO 7 N = 1,NDO
      K3PN = K3P1 + N
      K4PN = K4P1 + N
      A(K4PN) = A(K3PN)
 7        A(K3PN) = YP(N)
 8      A(K4P1) = A(K3P1)
 9      A(K3P1) = VAL
!
! NORMAL RETURN OR ABORT DUE TO ERROR
!
 10   YFIT = VAL
  return
!
   11 WRITE (XERN1, '(I8)') L
  WRITE (XERN2, '(I8)') NORD
  call XERMSG ('SLATEC', 'PVALUE', &
     'THE ORDER OF POLYNOMIAL EVALUATION, L = ' // XERN1 // &
     ' REQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD = ' // XERN2 // &
     ', COMPUTED BY POLFIT -- EXECUTION TERMINATED.', 8, 2)
  return
!
   12 call XERMSG ('SLATEC', 'PVALUE', &
     'INVALID INPUT PARAMETER.  ORDER OF POLYNOMIAL EVALUATION ' // &
     'REQUESTED IS NEGATIVE -- EXECUTION TERMINATED.', 2, 2)
  return
end
