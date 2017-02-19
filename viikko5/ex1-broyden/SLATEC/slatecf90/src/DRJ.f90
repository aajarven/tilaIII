  DOUBLE PRECISION FUNCTION DRJ (X, Y, Z, P, IER)
!
!! DRJ computes the incomplete or complete elliptic integral of 3rd kind.
!
!***PURPOSE  Compute the incomplete or complete (X or Y or Z is zero)
!            elliptic integral of the 3rd kind.  For X, Y, and Z non-
!            negative, at most one of them zero, and P positive,
!             RJ(X,Y,Z,P) = Integral from zero to infinity of
!                              -1/2     -1/2     -1/2     -1
!                    (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt.
!
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      DOUBLE PRECISION (RJ-S, DRJ-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE THIRD KIND,
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
!   1.     DRJ
!          Standard FORTRAN function routine
!          Double precision version
!          The routine calculates an approximation result to
!          DRJ(X,Y,Z,P) = Integral from zero to infinity of
!
!                                -1/2     -1/2     -1/2     -1
!                      (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt,
!
!          where X, Y, and Z are nonnegative, at most one of them is
!          zero, and P is positive.  If X or Y or Z is zero, the
!          integral is COMPLETE.  The duplication theorem is iterated
!          until the variables are nearly equal, and the function is
!          then expanded in Taylor series to fifth order.
!
!
!   2.     Calling Sequence
!          DRJ( X, Y, Z, P, IER )
!
!          Parameters on Entry
!          Values assigned by the calling routine
!
!          X      - Double precision, nonnegative variable
!
!          Y      - Double precision, nonnegative variable
!
!          Z      - Double precision, nonnegative variable
!
!          P      - Double precision, positive variable
!
!
!          On  Return    (values assigned by the DRJ routine)
!
!          DRJ     - Double precision approximation to the integral
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine. It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!
!          X, Y, Z, P are unaltered.
!
!
!   3.    Error Messages
!
!         Value of IER assigned by the DRJ routine
!
!              Value assigned         Error Message printed
!              IER = 1                MIN(X,Y,Z)  <  0.0D0
!                  = 2                MIN(X+Y,X+Z,Y+Z,P)  <  LOLIM
!                  = 3                MAX(X,Y,Z,P)  >  UPLIM
!
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!
!          LOLIM and UPLIM determine the valid range of X, Y, Z, and P
!
!          LOLIM is not less than the cube root of the value
!          of LOLIM used in the routine for DRC.
!
!          UPLIM is not greater than 0.3 times the cube root of
!          the value of UPLIM used in the routine for DRC.
!
!
!                     Acceptable values for:   LOLIM      UPLIM
!                     IBM 360/370 SERIES   :   2.0D-26     3.0D+24
!                     CDC 6000/7000 SERIES :   5.0D-98     3.0D+106
!                     UNIVAC 1100 SERIES   :   5.0D-103    6.0D+101
!                     CRAY                 :   1.32D-822   1.4D+821
!                     VAX 11 SERIES        :   2.5D-13     9.0D+11
!
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 the value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!
!
!
!          Relative error due to truncation of the series for DRJ
!          is less than 3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2.
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
!          Sample choices:  ERRTOL   Relative truncation
!                                    error less than
!                           1.0D-3    4.0D-18
!                           3.0D-3    3.0D-15
!                           1.0D-2    4.0D-12
!                           3.0D-2    3.0D-9
!                           1.0D-1    4.0D-6
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   DRJ Special Comments
!
!
!     Check by addition theorem: DRJ(X,X+Z,X+W,X+P)
!     + DRJ(Y,Y+Z,Y+W,Y+P) + (A-B) * DRJ(A,B,B,A) + 3.0D0 / SQRT(A)
!     = DRJ(0,Z,W,P), where X,Y,Z,W,P are positive and X * Y
!     = Z * W,  A = P * P * (X+Y+Z+W),  B = P * (P+X) * (P+Y),
!     and B - A = P * (P-Z) * (P-W).  The sum of the third and
!     fourth terms on the left side is 3.0D0 * DRC(A,B).
!
!
!          On Input:
!
!     X, Y, Z, and P are the variables in the integral DRJ(X,Y,Z,P).
!
!
!          On Output:
!
!
!          X, Y, Z, P are unaltered.
!
!          ********************************************************
!
!          WARNING: Changes in the program may improve speed at the
!                   expense of robustness.
!
!    -------------------------------------------------------------------
!
!
!   Special double precision functions via DRJ and DRF
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 3rd kind
!                  -----------------------------------------
!
!
!                          PHI         2         -1
!             P(PHI,K,N) = INT (1+N SIN (THETA) )   *
!                           0
!
!
!                                  2    2         -1/2
!                             *(1-K  SIN (THETA) )     D THETA
!
!
!                                           2          2   2
!                        = SIN (PHI) DRF(COS (PHI), 1-K SIN (PHI),1)
!
!                                   3             2         2   2
!                         -(N/3) SIN (PHI) DRJ(COS (PHI),1-K SIN (PHI),
!
!                                  2
!                         1,1+N SIN (PHI))
!
!
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 3rd kind
!                  -----------------------------------------
!
!
!                                            22    2
!                  EL3(X,KC,P) = X DRF(1,1+KC X ,1+X ) +
!
!                                            3           22    2     2
!                               +(1/3)(1-P) X  DRJ(1,1+KC X ,1+X ,1+PX )
!
!
!                                           2
!                  CEL(KC,P,A,B) = A RF(0,KC ,1) +
!
!
!                                                      2
!                                 +(1/3)(B-PA) DRJ(0,KC ,1,P)
!
!
!                  Heuman's LAMBDA function
!                  -----------------------------------------
!
!
!                                2                      2      2    1/2
!                  L(A,B,P) =(COS (A)SIN(B)COS(B)/(1-COS (A)SIN (B))   )
!
!                                            2         2       2
!                            *(SIN(P) DRF(COS (P),1-SIN (A) SIN (P),1)
!
!                                 2       3            2       2
!                            +(SIN (A) SIN (P)/(3(1-COS (A) SIN (B))))
!
!                                    2         2       2
!                            *DRJ(COS (P),1-SIN (A) SIN (P),1,1-
!
!                                2       2          2       2
!                            -SIN (A) SIN (P)/(1-COS (A) SIN (B))))
!
!
!
!                  (PI/2) LAMBDA0(A,B) =L(A,B,PI/2) =
!
!                   2                         2       2    -1/2
!              = COS (A)  SIN(B) COS(B) (1-COS (A) SIN (B))
!
!                           2                  2       2
!                 *DRF(0,COS (A),1) + (1/3) SIN (A) COS (A)
!
!                                      2       2    -3/2
!                 *SIN(B) COS(B) (1-COS (A) SIN (B))
!
!                           2         2       2          2       2
!                 *DRJ(0,COS (A),1,COS (A) COS (B)/(1-COS (A) SIN (B)))
!
!
!                  Jacobi ZETA function
!                  -----------------------------------------
!
!                        2                     2   2    1/2
!             Z(B,K) = (K/3) SIN(B) COS(B) (1-K SIN (B))
!
!
!                                  2      2   2                 2
!                        *DRJ(0,1-K ,1,1-K SIN (B)) / DRF (0,1-K ,1)
!
!
!  ---------------------------------------------------------------------
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
!***ROUTINES CALLED  D1MACH, DRC, XERMSG
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
!           editorial changes.  (RWC)).
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DRJ
  INTEGER IER
  CHARACTER*16 XERN3, XERN4, XERN5, XERN6, XERN7
  DOUBLE PRECISION ALFA, BETA, C1, C2, C3, C4, EA, EB, EC, E2, E3
  DOUBLE PRECISION LOLIM, UPLIM, EPSLON, ERRTOL, D1MACH
  DOUBLE PRECISION LAMDA, MU, P, PN, PNDEV
  DOUBLE PRECISION POWER4, DRC, SIGMA, S1, S2, S3, X, XN, XNDEV
  DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
   ZNROOT
  LOGICAL FIRST
  SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,C4,FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DRJ
  if (FIRST) THEN
     ERRTOL = (D1MACH(3)/3.0D0)**(1.0D0/6.0D0)
     LOLIM  = (5.0D0 * D1MACH(1))**(1.0D0/3.0D0)
     UPLIM  = 0.30D0*( D1MACH(2) / 5.0D0)**(1.0D0/3.0D0)
!
     C1 = 3.0D0/14.0D0
     C2 = 1.0D0/3.0D0
     C3 = 3.0D0/22.0D0
     C4 = 3.0D0/26.0D0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  DRJ = 0.0D0
  if (MIN(X,Y,Z) < 0.0D0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     call XERMSG ('SLATEC', 'DRJ', &
        'MIN(X,Y,Z) < 0 WHERE X = ' // XERN3 // ' Y = ' // XERN4 // &
        ' AND Z = ' // XERN5, 1, 1)
     return
  end if
!
  if (MAX(X,Y,Z,P) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') P
     WRITE (XERN7, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'DRJ', &
        'MAX(X,Y,Z,P) > UPLIM WHERE X = ' // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' P = ' // XERN6 // &
        ' AND UPLIM = ' // XERN7, 3, 1)
     return
  end if
!
  if (MIN(X+Y,X+Z,Y+Z,P) < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') P
     WRITE (XERN7, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'RJ', &
        'MIN(X+Y,X+Z,Y+Z,P) < LOLIM WHERE X = ' // XERN3 // &
        ' Y = ' // XERN4 // ' Z = '  // XERN5 // ' P = ' // XERN6 // &
        ' AND LOLIM = ', 2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
  ZN = Z
  PN = P
  SIGMA = 0.0D0
  POWER4 = 1.0D0
!
   30 MU = (XN+YN+ZN+PN+PN)*0.20D0
  XNDEV = (MU-XN)/MU
  YNDEV = (MU-YN)/MU
  ZNDEV = (MU-ZN)/MU
  PNDEV = (MU-PN)/MU
  EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV), ABS(PNDEV))
  if (EPSLON < ERRTOL) go to 40
  XNROOT =  SQRT(XN)
  YNROOT =  SQRT(YN)
  ZNROOT =  SQRT(ZN)
  LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
  ALFA = PN*(XNROOT+YNROOT+ZNROOT) + XNROOT*YNROOT*ZNROOT
  ALFA = ALFA*ALFA
  BETA = PN*(PN+LAMDA)*(PN+LAMDA)
  SIGMA = SIGMA + POWER4*DRC(ALFA,BETA,IER)
  POWER4 = POWER4*0.250D0
  XN = (XN+LAMDA)*0.250D0
  YN = (YN+LAMDA)*0.250D0
  ZN = (ZN+LAMDA)*0.250D0
  PN = (PN+LAMDA)*0.250D0
  go to 30
!
   40 EA = XNDEV*(YNDEV+ZNDEV) + YNDEV*ZNDEV
  EB = XNDEV*YNDEV*ZNDEV
  EC = PNDEV*PNDEV
  E2 = EA - 3.0D0*EC
  E3 = EB + 2.0D0*PNDEV*(EA-EC)
  S1 = 1.0D0 + E2*(-C1+0.750D0*C3*E2-1.50D0*C4*E3)
  S2 = EB*(0.50D0*C2+PNDEV*(-C3-C3+PNDEV*C4))
  S3 = PNDEV*EA*(C2-PNDEV*C3) - C2*PNDEV*EC
  DRJ = 3.0D0*SIGMA + POWER4*(S1+S2+S3)/(MU* SQRT(MU))
  return
end
