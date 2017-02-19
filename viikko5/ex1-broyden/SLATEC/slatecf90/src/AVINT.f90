subroutine AVINT (X, Y, N, XLO, XUP, ANS, IERR)
!
!! AVINT integrates a function tabulated at arbitrarily spaced abscissas...
!  using overlapping parabolas.
!
!***LIBRARY   SLATEC
!***CATEGORY  H2A1B2
!***TYPE      SINGLE PRECISION (AVINT-S, DAVINT-D)
!***KEYWORDS  INTEGRATION, QUADRATURE, TABULATED DATA
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!         AVINT integrates a function tabulated at arbitrarily spaced
!         abscissas.  The limits of integration need not coincide
!         with the tabulated abscissas.
!
!         A method of overlapping parabolas fitted to the data is used
!         provided that there are at least 3 abscissas between the
!         limits of integration.  AVINT also handles two special cases.
!         If the limits of integration are equal, AVINT returns a result
!         of zero regardless of the number of tabulated values.
!         If there are only two function values, AVINT uses the
!         trapezoid rule.
!
!     Description of Parameters
!         The user must dimension all arrays appearing in the call list
!              X(N), Y(N).
!
!         Input--
!         X    - real array of abscissas, which must be in increasing
!                order.
!         Y    - real array of functional values. i.e., Y(I)=FUNC(X(I)).
!         N    - the integer number of function values supplied.
!                N  >=  2 unless XLO = XUP.
!         XLO  - real lower limit of integration.
!         XUP  - real upper limit of integration.
!                Must have XLO  <=  XUP.
!
!         Output--
!         ANS  - computed approximate value of integral
!         IERR - a status code
!              --normal code
!                =1 means the requested integration was performed.
!              --abnormal codes
!                =2 means XUP was less than XLO.
!                =3 means the number of X(I) between XLO and XUP
!                   (inclusive) was less than 3 and neither of the two
!                   special cases described in the Abstract occurred.
!                   No integration was performed.
!                =4 means the restriction X(I+1)  >  X(I) was violated.
!                =5 means the number N of function values was  <  2.
!                ANS is set to zero if IERR=2,3,4,or 5.
!
!     AVINT is documented completely in SC-M-69-335
!     Original program from "Numerical Integration" by Davis &
!     Rabinowitz.
!     Adaptation and modifications for Sandia Mathematical Program
!     Library by Rondall E. Jones.
!
!***REFERENCES  R. E. Jones, Approximate integrator of functions
!                 tabulated at arbitrarily spaced abscissas,
!                 Report SC-M-69-335, Sandia Laboratories, 1969.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   690901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  AVINT
!
  DOUBLE PRECISION R3,RP5,SUM,SYL,SYL2,SYL3,SYU,SYU2,SYU3,X1,X2,X3 &
  ,X12,X13,X23,TERM1,TERM2,TERM3,A,B,C,CA,CB,CC
  DIMENSION X(*),Y(*)
!***FIRST EXECUTABLE STATEMENT  AVINT
  IERR=1
  ANS =0.0
  if (XLO-XUP) 3,100,200
    3 if (N < 2) go to 215
  DO 5 I=2,N
  if (X(I) <= X(I-1)) go to 210
  if (X(I) > XUP) go to 6
    5 CONTINUE
    6 CONTINUE
  if (N >= 3) go to 9
!
!     SPECIAL N=2 CASE
  SLOPE = (Y(2)-Y(1))/(X(2)-X(1))
  FL = Y(1) + SLOPE*(XLO-X(1))
  FR = Y(2) + SLOPE*(XUP-X(2))
  ANS = 0.5*(FL+FR)*(XUP-XLO)
  return
    9 CONTINUE
  if (X(N-2) < XLO)  go to 205
  if (X(3) > XUP)    go to 205
  I = 1
   10 if (X(I) >= XLO) go to 15
  I = I+1
  go to 10
   15 INLFT = I
  I = N
   20 if (X(I) <= XUP) go to 25
  I = I-1
  go to 20
   25 INRT = I
  if ((INRT-INLFT) < 2) go to 205
  ISTART = INLFT
  if (INLFT == 1) ISTART = 2
  ISTOP  = INRT
  if (INRT == N)  ISTOP  = N-1
!
  R3 = 3.0D0
  RP5= 0.5D0
  SUM = 0.0
  SYL = XLO
  SYL2= SYL*SYL
  SYL3= SYL2*SYL
!
  DO 50 I=ISTART,ISTOP
  X1 = X(I-1)
  X2 = X(I)
  X3 = X(I+1)
  X12 = X1-X2
  X13 = X1-X3
  X23 = X2-X3
  TERM1 = DBLE(Y(I-1))/(X12*X13)
  TERM2 =-DBLE(Y(I)) /(X12*X23)
  TERM3 = DBLE(Y(I+1))/(X13*X23)
  A = TERM1+TERM2+TERM3
  B = -(X2+X3)*TERM1 - (X1+X3)*TERM2 - (X1+X2)*TERM3
  C = X2*X3*TERM1 + X1*X3*TERM2 + X1*X2*TERM3
  if (I-ISTART) 30,30,35
   30 CA = A
  CB = B
  CC = C
  go to 40
   35 CA = 0.5*(A+CA)
  CB = 0.5*(B+CB)
  CC = 0.5*(C+CC)
   40 SYU = X2
  SYU2= SYU*SYU
  SYU3= SYU2*SYU
  SUM = SUM + CA*(SYU3-SYL3)/R3  + CB*RP5*(SYU2-SYL2) + CC*(SYU-SYL)
  CA  = A
  CB  = B
  CC  = C
  SYL = SYU
  SYL2= SYU2
  SYL3= SYU3
   50 CONTINUE
  SYU = XUP
  ANS = SUM + CA*(SYU**3-SYL3)/R3 + CB*RP5*(SYU**2-SYL2) &
    + CC*(SYU-SYL)
  100 RETURN
  200 IERR=2
  call XERMSG ('SLATEC', 'AVINT', &
     'THE UPPER LIMIT OF INTEGRATION WAS NOT GREATER THAN THE ' // &
     'LOWER LIMIT.', 4, 1)
  return
  205 IERR=3
  call XERMSG ('SLATEC', 'AVINT', &
     'THERE WERE LESS THAN THREE FUNCTION VALUES BETWEEN THE ' // &
     'LIMITS OF INTEGRATION.', 4, 1)
  return
  210 IERR=4
  call XERMSG ('SLATEC', 'AVINT', &
     'THE ABSCISSAS WERE NOT STRICTLY INCREASING.  MUST HAVE ' // &
     'X(I-1)  <  X(I) FOR ALL I.', 4, 1)
  return
  215 IERR=5
  call XERMSG ('SLATEC', 'AVINT', &
     'LESS THAN TWO FUNCTION VALUES WERE SUPPLIED.', 4, 1)
  return
end
