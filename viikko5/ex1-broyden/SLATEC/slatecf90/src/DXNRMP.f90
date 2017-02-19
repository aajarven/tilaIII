subroutine DXNRMP (NU, MU1, MU2, DARG, MODE, DPN, IPN, ISIG, &
     IERROR)
!
!! DXNRMP computes normalized Legendre polynomials.
!
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      DOUBLE PRECISION (XNRMP-S, DXNRMP-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
!           Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!
!        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS
!        (XNRMP is single-precision version)
!        DXNRMP calculates normalized Legendre polynomials of varying
!        order and fixed argument and degree. The order MU and degree
!        NU are non-negative integers and the argument is real. Because
!        the algorithm requires the use of numbers outside the normal
!        machine range, this subroutine employs a special arithmetic
!        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver,
!        and D.W. Lozier, Extended-Range Arithmetic and Normalized
!        Legendre Polynomials, ACM Transactions on Mathematical Soft-
!        ware, 93-105, March 1981, for a complete description of the
!        algorithm and special arithmetic. Also see program comments
!        in DXSET.
!
!        The normalized Legendre polynomials are multiples of the
!        associated Legendre polynomials of the first kind where the
!        normalizing coefficients are chosen so as to make the integral
!        from -1 to 1 of the square of each function equal to 1. See
!        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions,
!        McGraw-Hill, New York, 1960, p. 121.
!
!        The input values to DXNRMP are NU, MU1, MU2, DARG, and MODE.
!        These must satisfy
!          1. NU  >=  0 specifies the degree of the normalized Legendre
!             polynomial that is wanted.
!          2. MU1  >=  0 specifies the lowest-order normalized Legendre
!             polynomial that is wanted.
!          3. MU2  >=  MU1 specifies the highest-order normalized Leg-
!             endre polynomial that is wanted.
!         4a. MODE = 1 and -1.0D0  <=  DARG  <=  1.0D0 specifies that
!             Normalized Legendre(NU, MU, DARG) is wanted for MU = MU1,
!             MU1 + 1, ..., MU2.
!         4b. MODE = 2 and -3.14159...  <  DARG  <  3.14159... spec-
!             ifies that Normalized Legendre(NU, MU, COS(DARG)) is
!             wanted for MU = MU1, MU1 + 1, ..., MU2.
!
!        The output of DXNRMP consists of the two vectors DPN and IPN
!        and the error estimate ISIG. The computed values are stored as
!        extended-range numbers such that
!             (DPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,DX)
!             (DPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,DX)
!                .
!                .
!             (DPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,DX)
!        where K = MU2 - MU1 + 1 and DX = DARG or COS(DARG) according
!        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the
!        number of decimal digits lost through rounding errors in the
!        computation. For example if DARG is accurate to 12 significant
!        decimals, then the computed function values are accurate to
!        12 - ISIG significant decimals (except in neighborhoods of
!        zeros).
!
!        The interpretation of (DPN(I),IPN(I)) is DPN(I)*(IR**IPN(I))
!        where IR is the internal radix of the computer arithmetic. When
!        IPN(I) = 0 the value of the normalized Legendre polynomial is
!        contained entirely in DPN(I) and subsequent double-precision
!        computations can be performed without further consideration of
!        extended-range arithmetic. However, if IPN(I)  /=  0 the corre-
!        sponding value of the normalized Legendre polynomial cannot be
!        represented in double-precision because of overflow or under-
!        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the case
!        that IPN(I) is nonzero, the user could rewrite his/her program
!        to use extended range arithmetic.
!
!
!
!        The interpretation of (DPN(I),IPN(I)) can be changed to
!        DPN(I)*(10**IPN(I)) by calling the extended-range subroutine
!        DXCON. This should be done before printing the computed values.
!        As an example of usage, the Fortran coding
!              J = K
!              DO 20 I = 1, K
!              call DXCON(DPN(I), IPN(I),IERROR)
!              if (IERROR /= 0) RETURN
!              PRINT 10, DPN(I), IPN(I)
!           10 FORMAT(1X, D30.18 , I15)
!              if ((IPN(I)  ==  0) .OR. (J  <  K)) go to 20
!              J = I - 1
!           20 CONTINUE
!        will print all computed values and determine the largest J
!        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the
!        change of representation caused by calling DXCON, (DPN(I),
!        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent
!        extended-range computations.
!
!        IERROR is an error indicator. If no errors are detected,
!        IERROR=0 when control returns to the calling routine. If
!        an error is detected, IERROR is returned as nonzero. The
!        calling routine must check the value of IERROR.
!
!        If IERROR=212 or 213, invalid input was provided to DXNRMP.
!        If IERROR=201,202,203, or 204, invalid input was provided
!        to DXSET.
!        If IERROR=205 or 206, an internal consistency error occurred
!        in DXSET (probably due to a software malfunction in the
!        library routine I1MACH).
!        If IERROR=207, an overflow or underflow of an extended-range
!        number was detected in DXADJ.
!        If IERROR=208, an overflow or underflow of an extended-range
!        number was detected in DXC210.
!
!***SEE ALSO  DXSET
!***REFERENCES  Smith, Olver and Lozier, Extended-Range Arithmetic and
!                 Normalized Legendre Polynomials, ACM Trans on Math
!                 Softw, v 7, n 1, March 1981, pp 93--105.
!***ROUTINES CALLED  DXADD, DXADJ, DXRED, DXSET, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820712  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXNRMP
  INTEGER NU, MU1, MU2, MODE, IPN, ISIG
  DOUBLE PRECISION DARG, DPN
  DIMENSION DPN(*), IPN(*)
  DOUBLE PRECISION C1,C2,P,P1,P2,P3,S,SX,T,TX,X,DK
! call DXSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE DXSET
! LISTING FOR DETAILS)
!***FIRST EXECUTABLE STATEMENT  DXNRMP
  IERROR=0
  call DXSET (0, 0, 0.0D0, 0,IERROR)
  if (IERROR /= 0) RETURN
!
!        TEST FOR PROPER INPUT VALUES.
!
  if (NU < 0) go to 110
  if (MU1 < 0) go to 110
  if (MU1 > MU2) go to 110
  if (NU == 0) go to 90
  if (MODE < 1 .OR. MODE > 2) go to 110
  go to (10, 20), MODE
   10 if (ABS(DARG) > 1.0D0) go to 120
  if (ABS(DARG) == 1.0D0) go to 90
  X = DARG
  SX = SQRT((1.0D0+ABS(X))*((0.5D0-ABS(X))+0.5D0))
  TX = X/SX
  ISIG = LOG10(2.0D0*NU*(5.0D0+TX**2))
  go to 30
   20 if (ABS(DARG) > 4.0D0*ATAN(1.0D0)) go to 120
  if (DARG == 0.0D0) go to 90
  X = COS(DARG)
  SX = ABS(SIN(DARG))
  TX = X/SX
  ISIG = LOG10(2.0D0*NU*(5.0D0+ABS(DARG*TX)))
!
!        BEGIN CALCULATION
!
   30 MU = MU2
  I = MU2 - MU1 + 1
!
!        if MU > NU, NORMALIZED LEGENDRE(NU,MU,X)=0.
!
   40 if (MU <= NU) go to 50
  DPN(I) = 0.0D0
  IPN(I) = 0
  I = I - 1
  MU = MU - 1
  if (I  >  0) go to 40
  ISIG = 0
  go to 160
   50 MU = NU
!
!        P1 = 0. = NORMALIZED LEGENDRE(NU,NU+1,X)
!
  P1 = 0.0D0
  IP1 = 0
!
!        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X)
!
  P2 = 1.0D0
  IP2 = 0
  P3 = 0.5D0
  DK = 2.0D0
  DO 60 J=1,NU
    P3 = ((DK+1.0D0)/DK)*P3
    P2 = P2*SX
    call DXADJ(P2, IP2,IERROR)
    if (IERROR /= 0) RETURN
    DK = DK + 2.0D0
   60 CONTINUE
  P2 = P2*SQRT(P3)
  call DXADJ(P2, IP2,IERROR)
  if (IERROR /= 0) RETURN
  S = 2.0D0*TX
  T = 1.0D0/NU
  if (MU2 < NU) go to 70
  DPN(I) = P2
  IPN(I) = IP2
  I = I - 1
  if (I  ==  0) go to 140
!
!        RECURRENCE PROCESS
!
   70 P = MU*T
  C1 = 1.0D0/SQRT((1.0D0-P+T)*(1.0D0+P))
  C2 = S*P*C1*P2
  C1 = -SQRT((1.0D0+P+T)*(1.0D0-P))*C1*P1
  call DXADD(C2, IP2, C1, IP1, P, IP,IERROR)
  if (IERROR /= 0) RETURN
  MU = MU - 1
  if (MU > MU2) go to 80
!
!        STORE IN ARRAY DPN FOR RETURN TO CALLING ROUTINE.
!
  DPN(I) = P
  IPN(I) = IP
  I = I - 1
  if (I  ==  0) go to 140
   80 P1 = P2
  IP1 = IP2
  P2 = P
  IP2 = IP
  if (MU <= MU1) go to 140
  go to 70
!
!        SPECIAL CASE WHEN X=-1 OR +1, OR NU=0.
!
   90 K = MU2 - MU1 + 1
  DO 100 I=1,K
    DPN(I) = 0.0D0
    IPN(I) = 0
  100 CONTINUE
  ISIG = 0
  if (MU1 > 0) go to 160
  ISIG = 1
  DPN(1) = SQRT(NU+0.5D0)
  IPN(1) = 0
  if (MOD(NU,2) == 0) go to 160
  if (MODE == 1 .AND. DARG == 1.0D0) go to 160
  if (MODE == 2) go to 160
  DPN(1) = -DPN(1)
  go to 160
!
!          ERROR PRINTOUTS AND TERMINATION.
!
  110 call XERMSG ('SLATEC', 'DXNRMP', 'NU, MU1, MU2 or MODE not valid', &
               212, 1)
  IERROR=212
  return
  120 call XERMSG ('SLATEC', 'DXNRMP', 'DARG out of range', 213, 1)
  IERROR=213
  return
!
!        return TO CALLING PROGRAM
!
  140 K = MU2 - MU1 + 1
  DO 150 I=1,K
    call DXRED(DPN(I),IPN(I),IERROR)
    if (IERROR /= 0) RETURN
  150 CONTINUE
  160 RETURN
end
