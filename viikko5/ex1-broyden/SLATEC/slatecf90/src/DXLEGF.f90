subroutine DXLEGF (DNU1, NUDIFF, MU1, MU2, THETA, ID, PQA, IPQA, &
     IERROR)
!
!! DXLEGF computes normalized Legendre polynomials and associated Legendre ...
!  functions.
!
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      DOUBLE PRECISION (XLEGF-S, DXLEGF-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!
!   DXLEGF: Extended-range Double-precision Legendre Functions
!
!   A feature of the DXLEGF subroutine for Legendre functions is
! the use of extended-range arithmetic, a software extension of
! ordinary floating-point arithmetic that greatly increases the
! exponent range of the representable numbers. This avoids the
! need for scaling the solutions to lie within the exponent range
! of the most restrictive manufacturer's hardware. The increased
! exponent range is achieved by allocating an integer storage
! location together with each floating-point storage location.
!
!   The interpretation of the pair (X,I) where X is floating-point
! and I is integer is X*(IR**I) where IR is the internal radix of
! the computer arithmetic.
!
!   This subroutine computes one of the following vectors:
!
! 1. Legendre function of the first kind of negative order, either
!    a. P(-MU1,NU,X), P(-MU1-1,NU,X), ..., P(-MU2,NU,X) or
!    b. P(-MU,NU1,X), P(-MU,NU1+1,X), ..., P(-MU,NU2,X)
! 2. Legendre function of the second kind, either
!    a. Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X) or
!    b. Q(MU,NU1,X), Q(MU,NU1+1,X), ..., Q(MU,NU2,X)
! 3. Legendre function of the first kind of positive order, either
!    a. P(MU1,NU,X), P(MU1+1,NU,X), ..., P(MU2,NU,X) or
!    b. P(MU,NU1,X), P(MU,NU1+1,X), ..., P(MU,NU2,X)
! 4. Normalized Legendre polynomials, either
!    a. PN(MU1,NU,X), PN(MU1+1,NU,X), ..., PN(MU2,NU,X) or
!    b. PN(MU,NU1,X), PN(MU,NU1+1,X), ..., PN(MU,NU2,X)
!
! where X = COS(THETA).
!
!   The input values to DXLEGF are DNU1, NUDIFF, MU1, MU2, THETA,
! and ID. These must satisfy
!
!    DNU1 is DOUBLE PRECISION and greater than or equal to -0.5;
!    NUDIFF is INTEGER and non-negative;
!    MU1 is INTEGER and non-negative;
!    MU2 is INTEGER and greater than or equal to MU1;
!    THETA is DOUBLE PRECISION and in the half-open interval (0,PI/2];
!    ID is INTEGER and equal to 1, 2, 3 or 4;
!
! and  additionally either NUDIFF = 0 or MU2 = MU1.
!
!   If ID=1 and NUDIFF=0, a vector of type 1a above is computed
! with NU=DNU1.
!
!   If ID=1 and MU1=MU2, a vector of type 1b above is computed
! with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
!
!   If ID=2 and NUDIFF=0, a vector of type 2a above is computed
! with NU=DNU1.
!
!   If ID=2 and MU1=MU2, a vector of type 2b above is computed
! with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
!
!   If ID=3 and NUDIFF=0, a vector of type 3a above is computed
! with NU=DNU1.
!
!   If ID=3 and MU1=MU2, a vector of type 3b above is computed
! with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
!
!   If ID=4 and NUDIFF=0, a vector of type 4a above is computed
! with NU=DNU1.
!
!   If ID=4 and MU1=MU2, a vector of type 4b above is computed
! with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
!
!   In each case the vector of computed Legendre function values
! is returned in the extended-range vector (PQA(I),IPQA(I)). The
! length of this vector is either MU2-MU1+1 or NUDIFF+1.
!
!   Where possible, DXLEGF returns IPQA(I) as zero. In this case the
! value of the Legendre function is contained entirely in PQA(I),
! so it can be used in subsequent computations without further
! consideration of extended-range arithmetic. If IPQA(I) is nonzero,
! then the value of the Legendre function is not representable in
! floating-point because of underflow or overflow. The program that
! calls DXLEGF must test IPQA(I) to ensure correct usage.
!
!   IERROR is an error indicator. If no errors are detected, IERROR=0
! when control returns to the calling routine. If an error is detected,
! IERROR is returned as nonzero. The calling routine must check the
! value of IERROR.
!
!   If IERROR=210 or 211, invalid input was provided to DXLEGF.
!   If IERROR=201,202,203, or 204, invalid input was provided to DXSET.
!   If IERROR=205 or 206, an internal consistency error occurred in
! DXSET (probably due to a software malfunction in the library routine
! I1MACH).
!   If IERROR=207, an overflow or underflow of an extended-range number
! was detected in DXADJ.
!   If IERROR=208, an overflow or underflow of an extended-range number
! was detected in DXC210.
!
!***SEE ALSO  DXSET
!***REFERENCES  Olver and Smith, Associated Legendre Functions on the
!                 Cut, J Comp Phys, v 51, n 3, Sept 1983, pp 502--518.
!               Smith, Olver and Lozier, Extended-Range Arithmetic and
!                 Normalized Legendre Polynomials, ACM Trans on Math
!                 Softw, v 7, n 1, March 1981, pp 93--105.
!***ROUTINES CALLED  DXPMU, DXPMUP, DXPNRM, DXPQNU, DXQMU, DXQNU, DXRED,
!                    DXSET, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXLEGF
  DOUBLE PRECISION PQA,DNU1,DNU2,SX,THETA,X,PI2
  DIMENSION PQA(*),IPQA(*)
!
!***FIRST EXECUTABLE STATEMENT  DXLEGF
  IERROR=0
  call DXSET (0, 0, 0.0D0, 0,IERROR)
  if (IERROR /= 0) RETURN
  PI2=2.D0*ATAN(1.D0)
!
!        ZERO OUTPUT ARRAYS
!
  L=(MU2-MU1)+NUDIFF+1
  DO 290 I=1,L
  PQA(I)=0.D0
  290 IPQA(I)=0
!
!        CHECK FOR VALID INPUT VALUES
!
  if ( NUDIFF < 0) go to 400
  if ( DNU1 < -.5D0) go to 400
  if ( MU2 < MU1) go to 400
  if ( MU1 < 0) go to 400
  if ( THETA <= 0.D0.OR.THETA > PI2) go to 420
  if ( ID < 1.OR.ID > 4) go to 400
  if ( (MU1 /= MU2).AND.(NUDIFF > 0)) go to 400
!
!        if DNU1 IS NOT AN INTEGER, NORMALIZED P(MU,DNU,X)
!        CANNOT BE CALCULATED.  if DNU1 IS AN INTEGER AND
!        MU1 > DNU2 THEN ALL VALUES OF P(+MU,DNU,X) AND
!        NORMALIZED P(MU,NU,X) WILL BE ZERO.
!
  DNU2=DNU1+NUDIFF
  if ( (ID == 3).AND.(MOD(DNU1,1.D0) /= 0.D0)) go to 295
  if ( (ID == 4).AND.(MOD(DNU1,1.D0) /= 0.D0)) go to 400
  if ( (ID == 3.OR.ID == 4).AND.MU1 > DNU2) RETURN
  295 CONTINUE
!
  X=COS(THETA)
  SX=1.D0/SIN(THETA)
  if ( ID == 2) go to 300
  if ( MU2-MU1 <= 0) go to 360
!
!        FIXED NU, VARIABLE MU
!        call DXPMU TO CALCULATE P(-MU1,NU,X),....,P(-MU2,NU,X)
!
  call DXPMU(DNU1,DNU2,MU1,MU2,THETA,X,SX,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  go to 380
!
  300 if ( MU2 == MU1) go to 320
!
!        FIXED NU, VARIABLE MU
!        call DXQMU TO CALCULATE Q(MU1,NU,X),....,Q(MU2,NU,X)
!
  call DXQMU(DNU1,DNU2,MU1,MU2,THETA,X,SX,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  go to 390
!
!        FIXED MU, VARIABLE NU
!        call DXQNU TO CALCULATE Q(MU,DNU1,X),....,Q(MU,DNU2,X)
!
  320 call DXQNU(DNU1,DNU2,MU1,THETA,X,SX,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  go to 390
!
!        FIXED MU, VARIABLE NU
!        call DXPQNU TO CALCULATE P(-MU,DNU1,X),....,P(-MU,DNU2,X)
!
  360 call DXPQNU(DNU1,DNU2,MU1,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
!
!        if ID = 3, TRANSFORM P(-MU,NU,X) VECTOR INTO
!        P(MU,NU,X) VECTOR.
!
  380 if ( ID == 3) call DXPMUP(DNU1,DNU2,MU1,MU2,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
!
!        if ID = 4, TRANSFORM P(-MU,NU,X) VECTOR INTO
!        NORMALIZED P(MU,NU,X) VECTOR.
!
  if ( ID == 4) call DXPNRM(DNU1,DNU2,MU1,MU2,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
!
!        PLACE RESULTS IN REDUCED FORM if POSSIBLE
!        AND RETURN TO MAIN PROGRAM.
!
  390 DO 395 I=1,L
  call DXRED(PQA(I),IPQA(I),IERROR)
  if (IERROR /= 0) RETURN
  395 CONTINUE
  return
!
!        *****     ERROR TERMINATION     *****
!
  400 call XERMSG ('SLATEC', 'DXLEGF', &
               'DNU1, NUDIFF, MU1, MU2, or ID not valid', 210, 1)
  IERROR=210
  return
  420 call XERMSG ('SLATEC', 'DXLEGF', 'THETA out of range', 211, 1)
  IERROR=211
  return
end
