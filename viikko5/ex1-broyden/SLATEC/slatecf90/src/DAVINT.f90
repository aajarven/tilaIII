subroutine DAVINT (X, Y, N, XLO, XUP, ANS, IERR)
!
!! DAVINT integrates a function tabulated at arbitrarily spaced abscissas...
!  using overlapping parabolas.
!
!***LIBRARY   SLATEC
!***CATEGORY  H2A1B2
!***TYPE      DOUBLE PRECISION (AVINT-S, DAVINT-D)
!***KEYWORDS  INTEGRATION, QUADRATURE, TABULATED DATA
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!         DAVINT integrates a function tabulated at arbitrarily spaced
!         abscissas.  The limits of integration need not coincide
!         with the tabulated abscissas.
!
!         A method of overlapping parabolas fitted to the data is used
!         provided that there are at least 3 abscissas between the
!         limits of integration.  DAVINT also handles two special cases.
!         If the limits of integration are equal, DAVINT returns a
!         result of zero regardless of the number of tabulated values.
!         If there are only two function values, DAVINT uses the
!         trapezoid rule.
!
!     Description of Parameters
!         The user must dimension all arrays appearing in the call list
!              X(N), Y(N)
!
!         Input--
!      X    - DOUBLE PRECISION array of abscissas, which must be in
!             increasing order.
!      Y    - DOUBLE PRECISION array of function values. i.e.,
!                Y(I)=FUNC(X(I))
!      N    - The integer number of function values supplied.
!                N  >=  2 unless XLO = XUP.
!      XLO  - DOUBLE PRECISION lower limit of integration
!      XUP  - DOUBLE PRECISION upper limit of integration.  Must have
!              XLO <= XUP
!
!         Output--
!      ANS  - Double Precision computed approximate value of integral
!      IERR - A status code
!           --Normal Code
!                =1 Means the requested integration was performed.
!           --Abnormal Codes
!                =2 Means XUP was less than XLO.
!                =3 Means the number of X(I) between XLO and XUP
!                   (inclusive) was less than 3 and neither of the two
!                   special cases described in the abstract occurred.
!                   No integration was performed.
!                =4 Means the restriction X(I+1) > X(I) was violated.
!                =5 Means the number N of function values was .lt. 2.
!                   ANS is set to zero if IERR=2,3,4,or 5.
!
!    DAVINT is documented completely in SC-M-69-335
!    Original program from *Numerical Integration* by Davis & Rabinowitz
!    Adaptation and modifications by Rondall E Jones.
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
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DAVINT
!
  INTEGER I, IERR, INLFT, INRT, ISTART, ISTOP, N
  DOUBLE PRECISION A, ANS, B, C, CA, CB, CC, FL, FR, R3, RP5, &
       SLOPE, SUM, SYL, SYL2, SYL3, SYU, SYU2, SYU3, TERM1, TERM2, &
       TERM3, X, X1, X12, X13, X2, X23, X3, XLO, XUP, Y
  DIMENSION X(*),Y(*)
!     BEGIN BLOCK PERMITTING ...EXITS TO 190
!        BEGIN BLOCK PERMITTING ...EXITS TO 180
!***FIRST EXECUTABLE STATEMENT  DAVINT
        IERR = 1
        ANS = 0.0D0
        if (XLO  >  XUP) go to 160
           if (XLO  ==  XUP) go to 150
              if (N  >=  2) go to 10
                 IERR = 5
                 call XERMSG ('SLATEC', 'DAVINT', &
                    'LESS THAN TWO FUNCTION VALUES WERE SUPPLIED.', &
                    4, 1)
!     ...............EXIT
                 go to 190
   10             CONTINUE
              DO 20 I = 2, N
!        ............EXIT
                 if (X(I)  <=  X(I-1)) go to 180
!                 ...EXIT
                 if (X(I)  >  XUP) go to 30
   20             CONTINUE
   30             CONTINUE
              if (N  >=  3) go to 40
!
!                    SPECIAL N=2 CASE
                 SLOPE = (Y(2) - Y(1))/(X(2) - X(1))
                 FL = Y(1) + SLOPE*(XLO - X(1))
                 FR = Y(2) + SLOPE*(XUP - X(2))
                 ANS = 0.5D0*(FL + FR)*(XUP - XLO)
!     ...............EXIT
                 go to 190
   40             CONTINUE
              if (X(N-2)  >=  XLO) go to 50
                 IERR = 3
                 call XERMSG ('SLATEC', 'DAVINT', &
                    'THERE WERE LESS THAN THREE FUNCTION VALUES ' // &
                    'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
!     ...............EXIT
                 go to 190
   50             CONTINUE
              if (X(3)  <=  XUP) go to 60
                 IERR = 3
                 call XERMSG ('SLATEC', 'DAVINT', &
                    'THERE WERE LESS THAN THREE FUNCTION VALUES ' // &
                    'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
!     ...............EXIT
                 go to 190
   60             CONTINUE
              I = 1
   70             if (X(I)  >=  XLO) go to 80
                 I = I + 1
              go to 70
   80             CONTINUE
              INLFT = I
              I = N
   90             if (X(I)  <=  XUP) go to 100
                 I = I - 1
              go to 90
  100             CONTINUE
              INRT = I
              if ((INRT - INLFT)  >=  2) go to 110
                 IERR = 3
                 call XERMSG ('SLATEC', 'DAVINT', &
                    'THERE WERE LESS THAN THREE FUNCTION VALUES ' // &
                    'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
!     ...............EXIT
                 go to 190
  110             CONTINUE
              ISTART = INLFT
              if (INLFT  ==  1) ISTART = 2
              ISTOP = INRT
              if (INRT  ==  N) ISTOP = N - 1
!
              R3 = 3.0D0
              RP5 = 0.5D0
              SUM = 0.0D0
              SYL = XLO
              SYL2 = SYL*SYL
              SYL3 = SYL2*SYL
!
              DO 140 I = ISTART, ISTOP
                 X1 = X(I-1)
                 X2 = X(I)
                 X3 = X(I+1)
                 X12 = X1 - X2
                 X13 = X1 - X3
                 X23 = X2 - X3
                 TERM1 = Y(I-1)/(X12*X13)
                 TERM2 = -Y(I)/(X12*X23)
                 TERM3 = Y(I+1)/(X13*X23)
                 A = TERM1 + TERM2 + TERM3
                 B = -(X2 + X3)*TERM1 - (X1 + X3)*TERM2 &
                     - (X1 + X2)*TERM3
                 C = X2*X3*TERM1 + X1*X3*TERM2 + X1*X2*TERM3
                 if (I  >  ISTART) go to 120
                    CA = A
                    CB = B
                    CC = C
                 go to 130
  120                CONTINUE
                    CA = 0.5D0*(A + CA)
                    CB = 0.5D0*(B + CB)
                    CC = 0.5D0*(C + CC)
  130                CONTINUE
                 SYU = X2
                 SYU2 = SYU*SYU
                 SYU3 = SYU2*SYU
                 SUM = SUM + CA*(SYU3 - SYL3)/R3 &
                       + CB*RP5*(SYU2 - SYL2) + CC*(SYU - SYL)
                 CA = A
                 CB = B
                 CC = C
                 SYL = SYU
                 SYL2 = SYU2
                 SYL3 = SYU3
  140             CONTINUE
              SYU = XUP
              ANS = SUM + CA*(SYU**3 - SYL3)/R3 &
                    + CB*RP5*(SYU**2 - SYL2) + CC*(SYU - SYL)
  150          CONTINUE
        go to 170
  160       CONTINUE
           IERR = 2
           call XERMSG ('SLATEC', 'DAVINT', &
              'THE UPPER LIMIT OF INTEGRATION WAS NOT GREATER ' // &
              'THAN THE LOWER LIMIT.', 4, 1)
  170       CONTINUE
!     ......EXIT
        go to 190
  180    CONTINUE
     IERR = 4
     call XERMSG ('SLATEC', 'DAVINT', &
        'THE ABSCISSAS WERE NOT STRICTLY INCREASING.  MUST HAVE ' // &
        'X(I-1)  <  X(I) FOR ALL I.', 4, 1)
  190 CONTINUE
  return
end
