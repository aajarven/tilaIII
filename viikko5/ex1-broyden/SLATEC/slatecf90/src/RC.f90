FUNCTION RC (X, Y, IER)
!
!! RC approximates the elliptic integral RC(X,Y).
!
!***PURPOSE  Calculate an approximation to
!             RC(X,Y) = Integral from zero to infinity of
!                              -1/2     -1
!                    (1/2)(t+X)    (t+Y)  dt,
!            where X is nonnegative and Y is positive.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      SINGLE PRECISION (RC-S, DRC-D)
!***KEYWORDS  DUPLICATION THEOREM, ELEMENTARY FUNCTIONS,
!             ELLIPTIC INTEGRAL, TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     RC
!          Standard FORTRAN function routine
!          Single precision version
!          The routine calculates an approximation to
!           RC(X,Y) = Integral from zero to infinity of
!
!                              -1/2     -1
!                    (1/2)(t+X)    (t+Y)  dt,
!
!          where X is nonnegative and Y is positive.  The duplication
!          theorem is iterated until the variables are nearly equal,
!          and the function is then expanded in Taylor series to fifth
!          order.  Logarithmic, inverse circular, and inverse hyper-
!          bolic functions can be expressed in terms of RC.
!
!
!   2.     Calling Sequence
!          RC( X, Y, IER )
!
!          Parameters on Entry
!          Values assigned by the calling routine
!
!          X      - Single precision, nonnegative variable
!
!          Y      - Single precision, positive variable
!
!
!
!          On Return  (values assigned by the RC routine)
!
!          RC     - Single precision approximation to the integral
!
!          IER    - Integer to indicate normal or abnormal termination.
!
!                     IER = 0 Normal and reliable termination of the
!                             routine.  It is assumed that the requested
!                             accuracy has been achieved.
!
!                     IER > 0 Abnormal termination of the routine
!
!          X and Y are unaltered.
!
!
!   3.    Error Messages
!
!         Value of IER assigned by the RC routine
!
!                  Value Assigned         Error Message Printed
!                  IER = 1                X < 0.0E0.OR.Y <= 0.0E0
!                      = 2                X+Y < LOLIM
!                      = 3                MAX(X,Y)  >  UPLIM
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X and Y
!
!          LOLIM  - Lower limit of valid arguments
!
!                   Not less  than 5 * (machine minimum)  .
!
!          UPLIM  - Upper limit of valid arguments
!
!                   Not greater than (machine maximum) / 5 .
!
!
!                     Acceptable values for:   LOLIM       UPLIM
!                     IBM 360/370 SERIES   :   3.0E-78     1.0E+75
!                     CDC 6000/7000 SERIES :   1.0E-292    1.0E+321
!                     UNIVAC 1100 SERIES   :   1.0E-37     1.0E+37
!                     CRAY                 :   2.3E-2466   1.09E+2465
!                     VAX 11 SERIES        :   1.5E-38     3.0E+37
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!
!          ERRTOL  - Relative error due to truncation is less than
!                    16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).
!
!
!              The accuracy of the computed approximation to the inte-
!              gral can be controlled by choosing the value of ERRTOL.
!              Truncation of a Taylor series after terms of fifth order
!              introduces an error less than the amount shown in the
!              second column of the following table for each value of
!              ERRTOL in the first column.  In addition to the trunca-
!              tion error there will be round-off error, but in prac-
!              tice the total error from both sources is usually less
!              than the amount given in the table.
!
!
!
!          Sample Choices:  ERRTOL   Relative Truncation
!                                    error less than
!                           1.0E-3    2.0E-17
!                           3.0E-3    2.0E-14
!                           1.0E-2    2.0E-11
!                           3.0E-2    2.0E-8
!                           1.0E-1    2.0E-5
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   RC Special Comments
!
!
!
!
!                  Check: RC(X,X+Z) + RC(Y,Y+Z) = RC(0,Z)
!
!                  where X, Y, and Z are positive and X * Y = Z * Z
!
!
!          On Input:
!
!          X and Y are the variables in the integral RC(X,Y).
!
!          On Output:
!
!          X and Y are unaltered.
!
!
!
!                    RC(0,1/4)=RC(1/16,1/8)=PI=3.14159...
!
!                    RC(9/4,2)=LN(2)
!
!
!
!          ********************************************************
!
!          Warning: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!   --------------------------------------------------------------------
!
!   Special Functions via RC
!
!
!
!                  LN X                X  >  0
!
!                                            2
!                  LN(X) = (X-1) RC(((1+X)/2)  , X )
!
!
!   --------------------------------------------------------------------
!
!                  ARCSIN X            -1  <=  X  <=  1
!
!                                      2
!                  ARCSIN X = X RC (1-X  ,1 )
!
!   --------------------------------------------------------------------
!
!                  ARCCOS X            0  <=  X  <=  1
!
!
!                                     2      2
!                  ARCCOS X = SQRT(1-X ) RC(X  ,1 )
!
!   --------------------------------------------------------------------
!
!                  ARCTAN X            -INF  <  X  <  +INF
!
!                                       2
!                  ARCTAN X = X RC(1,1+X  )
!
!   --------------------------------------------------------------------
!
!                  ARCCOT X            0  <=  X  <  INF
!
!                                 2   2
!                  ARCCOT X = RC(X  ,X +1 )
!
!   --------------------------------------------------------------------
!
!                  ARCSINH X           -INF  <  X  <  +INF
!
!                                      2
!                  ARCSINH X = X RC(1+X  ,1 )
!
!   --------------------------------------------------------------------
!
!                  ARCCOSH X           X  >=  1
!
!                                    2        2
!                  ARCCOSH X = SQRT(X -1) RC(X  ,1 )
!
!   --------------------------------------------------------------------
!
!                  ARCTANH X           -1  <  X  <  1
!
!                                        2
!                  ARCTANH X = X RC(1,1-X  )
!
!   --------------------------------------------------------------------
!
!                  ARCCOTH X           X  >  1
!
!                                  2   2
!                  ARCCOTH X = RC(X  ,X -1 )
!
!   --------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Changed calls to XERMSG to standard form, and some
!           editorial changes.  (RWC))
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RC
  REAL RC
  CHARACTER*16 XERN3, XERN4, XERN5
  INTEGER IER
  REAL C1, C2, ERRTOL, LAMDA, LOLIM
  REAL MU, S, SN, UPLIM, X, XN, Y, YN
  LOGICAL FIRST
  SAVE ERRTOL,LOLIM,UPLIM,C1,C2,FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  RC
  if (FIRST) THEN
     ERRTOL = (R1MACH(3)/16.0E0)**(1.0E0/6.0E0)
     LOLIM  = 5.0E0 * R1MACH(1)
     UPLIM  = R1MACH(2) / 5.0E0
!
     C1 = 1.0E0/7.0E0
     C2 = 9.0E0/22.0E0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  RC = 0.0E0
  if (X < 0.0E0.OR.Y <= 0.0E0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     call XERMSG ('SLATEC', 'RC', &
        'X < 0 .OR. Y <= 0 WHERE X = ' // XERN3 // ' AND Y = ' // &
        XERN4, 1, 1)
     return
  end if
!
  if (MAX(X,Y) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'RC', &
        'MAX(X,Y) > UPLIM WHERE X = '  // XERN3 // ' Y = ' // &
        XERN4 // ' AND UPLIM = ' // XERN5, 3, 1)
     return
  end if
!
  if (X+Y < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'RC', &
        'X+Y < LOLIM WHERE X = ' // XERN3 // ' Y = ' // XERN4 // &
        ' AND LOLIM = ' // XERN5, 2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
!
   30 MU = (XN+YN+YN)/3.0E0
  SN = (YN+MU)/MU - 2.0E0
  if (ABS(SN) < ERRTOL) go to 40
  LAMDA = 2.0E0*SQRT(XN)*SQRT(YN) + YN
  XN = (XN+LAMDA)*0.250E0
  YN = (YN+LAMDA)*0.250E0
  go to 30
!
   40 S = SN*SN*(0.30E0+SN*(C1+SN*(0.3750E0+SN*C2)))
  RC = (1.0E0+S)/SQRT(MU)
  return
end
