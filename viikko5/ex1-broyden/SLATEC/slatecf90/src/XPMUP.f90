subroutine XPMUP (NU1, NU2, MU1, MU2, PQA, IPQA, IERROR)
!
!! XPMUP computes the values of Legendre functions for XLEGF.
!
!            This subroutine transforms an array of Legendre functions
!            of the first kind of negative order stored in array PQA
!            into Legendre functions of the first kind of positive
!            order stored in array PQA. The original array is destroyed.
!
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      SINGLE PRECISION (XPMUP-S, DXPMUP-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  Smith, John M., (NBS and George Mason University)
!***ROUTINES CALLED  XADJ
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  XPMUP
  REAL DMU,NU,NU1,NU2,PQA,PROD
  DIMENSION PQA(*),IPQA(*)
!***FIRST EXECUTABLE STATEMENT  XPMUP
  IERROR=0
  NU=NU1
  MU=MU1
  DMU=MU
  N=INT(NU2-NU1+.1)+(MU2-MU1)+1
  J=1
  if ( MOD(NU,1.) /= 0.) go to 210
  200 if ( DMU < NU+1.) go to 210
  PQA(J)=0.
  IPQA(J)=0
  J=J+1
  if ( J > N) RETURN
!        INCREMENT EITHER MU OR NU AS APPROPRIATE.
  if ( NU2-NU1 > .5) NU=NU+1.
  if ( MU2 > MU1) MU=MU+1
  go to 200
!
!        TRANSFORM P(-MU,NU,X) TO P(MU,NU,X) USING
!        P(MU,NU,X)=(NU-MU+1)*(NU-MU+2)*...*(NU+MU)*P(-MU,NU,X)*(-1)**MU
!
  210 PROD=1.
  IPROD=0
  K=2*MU
  if ( K == 0) go to 222
  DO 220 L=1,K
  PROD=PROD*(DMU-NU-L)
  220 call XADJ(PROD,IPROD,IERROR)
  if (IERROR /= 0) RETURN
  222 CONTINUE
  DO 240 I=J,N
  if ( MU == 0) go to 225
  PQA(I)=PQA(I)*PROD*(-1)**MU
  IPQA(I)=IPQA(I)+IPROD
  call XADJ(PQA(I),IPQA(I),IERROR)
  if (IERROR /= 0) RETURN
  225 if ( NU2-NU1 > .5) go to 230
  PROD=(DMU-NU)*PROD*(-DMU-NU-1.)
  call XADJ(PROD,IPROD,IERROR)
  if (IERROR /= 0) RETURN
  MU=MU+1
  DMU=DMU+1.
  go to 240
  230 PROD=PROD*(-DMU-NU-1.)/(DMU-NU-1.)
  call XADJ(PROD,IPROD,IERROR)
  if (IERROR /= 0) RETURN
  NU=NU+1.
  240 CONTINUE
  return
end
