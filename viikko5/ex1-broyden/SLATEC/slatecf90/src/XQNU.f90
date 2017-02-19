subroutine XQNU (NU1, NU2, MU1, THETA, X, SX, ID, PQA, IPQA, &
     IERROR)
!
!! XQNU computes the values of Legendre functions for XLEGF.
!
!            Method: backward nu-wise recurrence for Q(MU,NU,X) for
!            fixed mu to obtain Q(MU1,NU1,X), Q(MU1,NU1+1,X), ...,
!            Q(MU1,NU2,X).
!
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      SINGLE PRECISION (XQNU-S, DXQNU-D)
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
!***END PROLOGUE  XQNU
  DIMENSION PQA(*),IPQA(*)
  REAL DMU,NU,NU1,NU2,PQ,PQA,PQ1,PQ2,SX,X,X1,X2
  REAL THETA,PQL1,PQL2
!***FIRST EXECUTABLE STATEMENT  XQNU
  IERROR=0
  K=0
  PQ2=0.0
  IPQ2=0
  PQL2=0.0
  IPQL2=0
  if ( MU1 == 1) go to 290
  MU=0
!
!        call XPQNU TO OBTAIN Q(0.,NU2,X) AND Q(0.,NU2-1,X)
!
  call XPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  if ( MU1 == 0) RETURN
  K=(NU2-NU1+1.5)
  PQ2=PQA(K)
  IPQ2=IPQA(K)
  PQL2=PQA(K-1)
  IPQL2=IPQA(K-1)
  290 MU=1
!
!        call XPQNU TO OBTAIN Q(1.,NU2,X) AND Q(1.,NU2-1,X)
!
  call XPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
  if (IERROR /= 0) RETURN
  if ( MU1 == 1) RETURN
  NU=NU2
  PQ1=PQA(K)
  IPQ1=IPQA(K)
  PQL1=PQA(K-1)
  IPQL1=IPQA(K-1)
  300 MU=1
  DMU=1.
  320 CONTINUE
!
!        FORWARD RECURRENCE IN MU TO OBTAIN Q(MU1,NU2,X) AND
!              Q(MU1,NU2-1,X) USING
!              Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X)
!                   -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X)
!
!              FIRST FOR NU=NU2
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
  PQA(K)=PQ
  IPQA(K)=IPQ
  if ( K == 1) RETURN
  if ( NU < NU2) go to 340
!
!              THEN FOR NU=NU2-1
!
  NU=NU-1.
  PQ2=PQL2
  IPQ2=IPQL2
  PQ1=PQL1
  IPQ1=IPQL1
  K=K-1
  go to 300
!
!         BACKWARD RECURRENCE IN NU TO OBTAIN
!              Q(MU1,NU1,X),Q(MU1,NU1+1,X),....,Q(MU1,NU2,X)
!              USING
!              (NU-MU+1.)*Q(MU,NU+1,X)=
!                       (2.*NU+1.)*X*Q(MU,NU,X)-(NU+MU)*Q(MU,NU-1,X)
!
  340 PQ1=PQA(K)
  IPQ1=IPQA(K)
  PQ2=PQA(K+1)
  IPQ2=IPQA(K+1)
  350 if ( NU <= NU1) RETURN
  K=K-1
  X1=(2.*NU+1.)*X*PQ1/(NU+DMU)
  X2=-(NU-DMU+1.)*PQ2/(NU+DMU)
  call XADD(X1,IPQ1,X2,IPQ2,PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  call XADJ(PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  PQ2=PQ1
  IPQ2=IPQ1
  PQ1=PQ
  IPQ1=IPQ
  PQA(K)=PQ
  IPQA(K)=IPQ
  NU=NU-1.
  go to 350
end
