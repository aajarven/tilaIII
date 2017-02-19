subroutine XQMU (NU1, NU2, MU1, MU2, THETA, X, SX, ID, PQA, IPQA, &
     IERROR)
!
!! XQMU computes the values of Legendre functions for XLEGF.
!
!            Method: forward mu-wise recurrence for Q(MU,NU,X) for fixed
!            nu to obtain Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X).
!
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      SINGLE PRECISION (XQMU-S, DXQMU-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  Smith, John M., (NBS and George Mason University)
!***ROUTINES CALLED  XADD, XADJ, XPQNU
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  XQMU
  DIMENSION PQA(*),IPQA(*)
  REAL DMU,NU,NU1,NU2,PQ,PQA,PQ1,PQ2,SX,X,X1,X2
  REAL THETA
!***FIRST EXECUTABLE STATEMENT  XQMU
  IERROR=0
  MU=0
!
!        call XPQNU TO OBTAIN Q(0.,NU1,X)
!
  call XPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  PQ2=PQA(1)
  IPQ2=IPQA(1)
  MU=1
!
!        call XPQNU TO OBTAIN Q(1.,NU1,X)
!
  call XPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  NU=NU1
  K=0
  MU=1
  DMU=1.
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
  X1=-2.*DMU*X*SX*PQ1
  X2=(NU+DMU)*(NU-DMU+1.)*PQ2
  call XADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  call XADJ(PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  PQ2=PQ1
  IPQ2=IPQ1
  PQ1=PQ
  IPQ1=IPQ
  MU=MU+1
  DMU=DMU+1.
  if ( MU < MU1) go to 320
  K=K+1
  PQA(K)=PQ
  IPQA(K)=IPQ
  if ( MU2 > MU) go to 320
  330 RETURN
end
