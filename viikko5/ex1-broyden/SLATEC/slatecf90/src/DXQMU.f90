subroutine DXQMU (NU1, NU2, MU1, MU2, THETA, X, SX, ID, PQA, IPQA, &
     IERROR)
!
!! DXQMU computes the values of Legendre functions for DXLEGF.
!
!            Method: forward mu-wise recurrence for Q(MU,NU,X) for fixed
!            nu to obtain Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X).
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      DOUBLE PRECISION (XQMU-S, DXQMU-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  Smith, John M., (NBS and George Mason University)
!***ROUTINES CALLED  DXADD, DXADJ, DXPQNU
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXQMU
  DIMENSION PQA(*),IPQA(*)
  DOUBLE PRECISION DMU,NU,NU1,NU2,PQ,PQA,PQ1,PQ2,SX,X,X1,X2
  DOUBLE PRECISION THETA
!***FIRST EXECUTABLE STATEMENT  DXQMU
  IERROR=0
  MU=0
!
!        call DXPQNU TO OBTAIN Q(0.,NU1,X)
!
  call DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  PQ2=PQA(1)
  IPQ2=IPQA(1)
  MU=1
!
!        call DXPQNU TO OBTAIN Q(1.,NU1,X)
!
  call DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  NU=NU1
  K=0
  MU=1
  DMU=1.D0
  PQ1=PQA(1)
  IPQ1=IPQA(1)
  if ( MU1 > 0) go to 310
  K=K+1
  PQA(K)=PQ2
  IPQA(K)=IPQ2
  if ( MU2 < 1) go to 330
  310 if ( MU1 > 1) go to 320
  K=K+1
  PQA(K)=PQ1
  IPQA(K)=IPQ1
  if ( MU2 <= 1) go to 330
  320 CONTINUE
!
!        FORWARD RECURRENCE IN MU TO OBTAIN
!                  Q(MU1,NU,X),Q(MU1+1,NU,X),....,Q(MU2,NU,X) USING
!             Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X)
!                               -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X)
!
  X1=-2.D0*DMU*X*SX*PQ1
  X2=(NU+DMU)*(NU-DMU+1.D0)*PQ2
  call DXADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  call DXADJ(PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  PQ2=PQ1
  IPQ2=IPQ1
  PQ1=PQ
  IPQ1=IPQ
  MU=MU+1
  DMU=DMU+1.D0
  if ( MU < MU1) go to 320
  K=K+1
  PQA(K)=PQ
  IPQA(K)=IPQ
  if ( MU2 > MU) go to 320
  330 RETURN
end
