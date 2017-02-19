subroutine DXPNRM (NU1, NU2, MU1, MU2, PQA, IPQA, IERROR)
!
!! DXPNRM computes the values of Legendre functions for DXLEGF.
!
!            This subroutine transforms an array of Legendre functions
!            of the first kind of negative order stored in array PQA
!            into normalized Legendre polynomials stored in array PQA.
!            The original array is destroyed.
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      DOUBLE PRECISION (XPNRM-S, DXPNRM-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  Smith, John M., (NBS and George Mason University)
!***ROUTINES CALLED  DXADJ
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXPNRM
  DOUBLE PRECISION C1,DMU,NU,NU1,NU2,PQA,PROD
  DIMENSION PQA(*),IPQA(*)
!***FIRST EXECUTABLE STATEMENT  DXPNRM
  IERROR=0
  L=(MU2-MU1)+(NU2-NU1+1.5D0)
  MU=MU1
  DMU=MU1
  NU=NU1
!
!         if MU  > NU, NORM P =0.
!
  J=1
  500 if ( DMU <= NU) go to 505
  PQA(J)=0.D0
  IPQA(J)=0
  J=J+1
  if ( J > L) RETURN
!
!        INCREMENT EITHER MU OR NU AS APPROPRIATE.
!
  if ( MU2 > MU1) DMU=DMU+1.D0
  if ( NU2-NU1 > .5D0) NU=NU+1.D0
  go to 500
!
!         TRANSFORM P(-MU,NU,X) INTO NORMALIZED P(MU,NU,X) USING
!              NORM P(MU,NU,X)=
!                 SQRT((NU+.5)*FACTORIAL(NU+MU)/FACTORIAL(NU-MU))
!                              *P(-MU,NU,X)
!
  505 PROD=1.D0
  IPROD=0
  K=2*MU
  if ( K <= 0) go to 520
  DO 510 I=1,K
  PROD=PROD*SQRT(NU+DMU+1.D0-I)
  510 call DXADJ(PROD,IPROD,IERROR)
  if (IERROR /= 0) RETURN
  520 DO 540 I=J,L
  C1=PROD*SQRT(NU+.5D0)
  PQA(I)=PQA(I)*C1
  IPQA(I)=IPQA(I)+IPROD
  call DXADJ(PQA(I),IPQA(I),IERROR)
  if (IERROR /= 0) RETURN
  if ( NU2-NU1 > .5D0) go to 530
  if ( DMU >= NU) go to 525
  PROD=SQRT(NU+DMU+1.D0)*PROD
  if ( NU > DMU) PROD=PROD*SQRT(NU-DMU)
  call DXADJ(PROD,IPROD,IERROR)
  if (IERROR /= 0) RETURN
  MU=MU+1
  DMU=DMU+1.D0
  go to 540
  525 PROD=0.D0
  IPROD=0
  MU=MU+1
  DMU=DMU+1.D0
  go to 540
  530 PROD=SQRT(NU+DMU+1.D0)*PROD
  if ( NU /= DMU-1.D0) PROD=PROD/SQRT(NU-DMU+1.D0)
  call DXADJ(PROD,IPROD,IERROR)
  if (IERROR /= 0) RETURN
  NU=NU+1.D0
  540 CONTINUE
  return
end
