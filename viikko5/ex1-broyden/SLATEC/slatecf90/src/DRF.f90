  DOUBLE PRECISION FUNCTION DRF (X, Y, Z, IER)
!
!! DRF computes the incomplete or complete elliptic integral of 1st kind.
!
!***PURPOSE  Compute the incomplete or complete elliptic integral of the
!            1st kind.  For X, Y, and Z non-negative and at most one of
!            them zero, RF(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -1/2
!                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X, Y or Z is zero, the integral is complete.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      DOUBLE PRECISION (RF-S, DRF-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
!             TAYLOR SERIES
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
!   1.     DRF
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the first kind
!          Standard FORTRAN function routine
!          Double precision version
!          The routine calculates an approximation result to
!          DRF(X,Y,Z) = Integral from zero to infinity of
!
!                               -1/2     -1/2     -1/2
!                     (1/2)(t+X)    (t+Y)    (t+Z)    dt,
!
!          where X, Y, and Z are nonnegative and at most one of them
!          is zero.  If one of them  is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling sequence
!          DRF( X, Y, Z, IER )
!
!          Parameters On entry
!          Values assigned by the calling routine
!
!          X      - Double precision, nonnegative variable
!
!          Y      - Double precision, nonnegative variable
!
!          Z      - Double precision, nonnegative variable
!
!
!
!          On Return    (values assigned by the DRF routine)
!
!          DRF     - Double precision approximation to the integral
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine. It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!          X, Y, Z are unaltered.
!
!
!   3.    Error Messages
!
!
!         Value of IER assigned by the DRF routine
!
!                  Value assigned         Error Message Printed
!                  IER = 1                MIN(X,Y,Z)  <  0.0D0
!                      = 2                MIN(X+Y,X+Z,Y+Z)  <  LOLIM
!                      = 3                MAX(X,Y,Z)  >  UPLIM
!
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                   Not less than 5 * (machine minimum).
!
!          UPLIM  - Upper limit of valid arguments
!
!                   Not greater than (machine maximum) / 5.
!
!
!                     Acceptable values for:   LOLIM      UPLIM
!                     IBM 360/370 SERIES   :   3.0D-78     1.0D+75
!                     CDC 6000/7000 SERIES :   1.0D-292    1.0D+321
!                     UNIVAC 1100 SERIES   :   1.0D-307    1.0D+307
!                     CRAY                 :   2.3D-2466   1.09D+2465
!                     VAX 11 SERIES        :   1.5D-38     3.0D+37
!
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!
!
!          ERRTOL - Relative error due to truncation is less than
!                   ERRTOL ** 6 / (4 * (1-ERRTOL)  .
!
!
!
!        The accuracy of the computed approximation to the integral
!        can be controlled by choosing the value of ERRTOL.
!        Truncation of a Taylor series after terms of fifth order
!        introduces an error less than the amount shown in the
!        second column of the following table for each value of
!        ERRTOL in the first column.  In addition to the truncation
!        error there will be round-off error, but in practice the
!        total error from both sources is usually less than the
!        amount given in the table.
!
!
!
!
!
!          Sample choices:  ERRTOL   Relative Truncation
!                                    error less than
!                           1.0D-3    3.0D-19
!                           3.0D-3    2.0D-16
!                           1.0D-2    3.0D-13
!                           3.0D-2    2.0D-10
!                           1.0D-1    3.0D-7
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   DRF Special Comments
!
!
!
!          Check by addition theorem: DRF(X,X+Z,X+W) + DRF(Y,Y+Z,Y+W)
!          = DRF(0,Z,W), where X,Y,Z,W are positive and X * Y = Z * W.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral DRF(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, Z are unaltered.
!
!
!
!          ********************************************************
!
!          WARNING: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!   Special double precision functions via DRF
!
!
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 1st kind
!
!                  -----------------------------------------
!
!
!
!                                             2         2   2
!                  F(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1)
!
!
!                                  2
!                  K(K) = DRF(0,1-K ,1)
!
!
!                         PI/2     2   2      -1/2
!                       = INT  (1-K SIN (PHI) )   D PHI
!                          0
!
!
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 1st kind
!
!                  -----------------------------------------
!
!
!                                          22    2
!                  EL1(X,KC) = X DRF(1,1+KC X ,1+X )
!
!
!                  Lemniscate constant A
!
!                  -----------------------------------------
!
!
!                       1      4 -1/2
!                  A = INT (1-S )    DS = DRF(0,1,2) = DRF(0,2,1)
!                       0
!
!
!
!    -------------------------------------------------------------------
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
!***ROUTINES CALLED  D1MACH, XERMSG
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
!***END PROLOGUE  DRF
  CHARACTER*16 XERN3, XERN4, XERN5, XERN6
  INTEGER IER
  DOUBLE PRECISION LOLIM, UPLIM, EPSLON, ERRTOL, D1MACH
  DOUBLE PRECISION C1, C2, C3, E2, E3, LAMDA
  DOUBLE PRECISION MU, S, X, XN, XNDEV
  DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
   ZNROOT
  LOGICAL FIRST
  SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DRF
!
  if (FIRST) THEN
     ERRTOL = (4.0D0*D1MACH(3))**(1.0D0/6.0D0)
     LOLIM  = 5.0D0 * D1MACH(1)
     UPLIM  = D1MACH(2)/5.0D0
!
     C1 = 1.0D0/24.0D0
     C2 = 3.0D0/44.0D0
     C3 = 1.0D0/14.0D0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  DRF = 0.0D0
  if (MIN(X,Y,Z) < 0.0D0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     call XERMSG ('SLATEC', 'DRF', &
        'MIN(X,Y,Z) < 0 WHERE X = ' // XERN3 // ' Y = ' // XERN4 // &
        ' AND Z = ' // XERN5, 1, 1)
     return
  end if
!
  if (MAX(X,Y,Z) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'DRF', &
        'MAX(X,Y,Z) > UPLIM WHERE X = '  // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6, 3, 1)
     return
  end if
!
  if (MIN(X+Y,X+Z,Y+Z) < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'DRF', &
        'MIN(X+Y,X+Z,Y+Z) < LOLIM WHERE X = ' // XERN3 // &
        ' Y = ' // XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' // &
        XERN6, 2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
  ZN = Z
!
   30 MU = (XN+YN+ZN)/3.0D0
  XNDEV = 2.0D0 - (MU+XN)/MU
  YNDEV = 2.0D0 - (MU+YN)/MU
  ZNDEV = 2.0D0 - (MU+ZN)/MU
  EPSLON = MAX(ABS(XNDEV),ABS(YNDEV),ABS(ZNDEV))
  if (EPSLON < ERRTOL) go to 40
  XNROOT = SQRT(XN)
  YNROOT = SQRT(YN)
  ZNROOT = SQRT(ZN)
  LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
  XN = (XN+LAMDA)*0.250D0
  YN = (YN+LAMDA)*0.250D0
  ZN = (ZN+LAMDA)*0.250D0
  go to 30
!
   40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
  E3 = XNDEV*YNDEV*ZNDEV
  S  = 1.0D0 + (C1*E2-0.10D0-C2*E3)*E2 + C3*E3
  DRF = S/SQRT(MU)
!
  return
end
