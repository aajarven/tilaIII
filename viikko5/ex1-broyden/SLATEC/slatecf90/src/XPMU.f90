subroutine XPMU (NU1, NU2, MU1, MU2, THETA, X, SX, ID, PQA, IPQA, &
     IERROR)
!
!! XPMU computes the values of Legendre functions for XLEGF.
!
!            Method: backward mu-wise recurrence for P(-MU,NU,X) for
!            fixed nu to obtain P(-MU2,NU1,X), P(-(MU2-1),NU1,X), ...,
!            P(-MU1,NU1,X) and store in ascending mu order.
!
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      SINGLE PRECISION (XPMU-S, DXPMU-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  Smith, John M., (NBS and George Mason University)
!***ROUTINES CALLED  XADD, XADJ, XPQNU
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  XPMU
  REAL PQA,NU1,NU2,P0,X,SX,THETA,X1,X2
  DIMENSION PQA(*),IPQA(*)
!
!        call XPQNU TO OBTAIN P(-MU2,NU,X)
!
!***FIRST EXECUTABLE STATEMENT  XPMU
  IERROR=0
  call XPQNU(NU1,NU2,MU2,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  P0=PQA(1)
  IP0=IPQA(1)
  MU=MU2-1
!
!        call XPQNU TO OBTAIN P(-MU2-1,NU,X)
!
  call XPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  N=MU2-MU1+1
  PQA(N)=P0
  IPQA(N)=IP0
  if ( N == 1) go to 300
  PQA(N-1)=PQA(1)
  IPQA(N-1)=IPQA(1)
  if ( N == 2) go to 300
  J=N-2
  290 CONTINUE
!
!        BACKWARD RECURRENCE IN MU TO OBTAIN
!              P(-MU2,NU1,X),P(-(MU2-1),NU1,X),....P(-MU1,NU1,X)
!              USING
!              (NU-MU)*(NU+MU+1.)*P(-(MU+1),NU,X)=
!                2.*MU*X*SQRT((1./(1.-X**2))*P(-MU,NU,X)-P(-(MU-1),NU,X)
!
  X1=2.*MU*X*SX*PQA(J+1)
  X2=-(NU1-MU)*(NU1+MU+1.)*PQA(J+2)
  call XADD(X1,IPQA(J+1),X2,IPQA(J+2),PQA(J),IPQA(J),IERROR)
  if (IERROR /= 0) RETURN
  call XADJ(PQA(J),IPQA(J),IERROR)
  if (IERROR /= 0) RETURN
  if ( J == 1) go to 300
  J=J-1
  MU=MU-1
  go to 290
  300 RETURN
end
