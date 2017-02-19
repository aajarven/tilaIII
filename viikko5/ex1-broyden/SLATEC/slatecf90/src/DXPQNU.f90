subroutine DXPQNU (NU1, NU2, MU, THETA, ID, PQA, IPQA, IERROR)
!
!! DXPQNU computes the values of Legendre functions for DXLEGF.
!
!            This subroutine calculates initial values of P or Q using
!            power series, then performs forward nu-wise recurrence to
!            obtain P(-MU,NU,X), Q(0,NU,X), or Q(1,NU,X). The nu-wise
!            recurrence is stable for P for all mu and for Q for mu=0,1.
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      DOUBLE PRECISION (XPQNU-S, DXPQNU-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  Smith, John M., (NBS and George Mason University)
!***ROUTINES CALLED  DXADD, DXADJ, DXPSI
!***COMMON BLOCKS    DXBLK1
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXPQNU
  DOUBLE PRECISION A,NU,NU1,NU2,PQ,PQA,DXPSI,R,THETA,W,X,X1,X2,XS, &
   Y,Z
  DOUBLE PRECISION DI,DMU,PQ1,PQ2,FACTMU,FLOK
  DIMENSION PQA(*),IPQA(*)
  COMMON /DXBLK1/ NBITSF
  SAVE /DXBLK1/
!
!        J0, IPSIK, AND IPSIX ARE INITIALIZED IN THIS SUBROUTINE.
!        J0 IS THE NUMBER OF TERMS USED IN SERIES EXPANSION
!        IN SUBROUTINE DXPQNU.
!        IPSIK, IPSIX ARE VALUES OF K AND X RESPECTIVELY
!        USED IN THE CALCULATION OF THE DXPSI FUNCTION.
!
!***FIRST EXECUTABLE STATEMENT  DXPQNU
  IERROR=0
  J0=NBITSF
  IPSIK=1+(NBITSF/10)
  IPSIX=5*IPSIK
  IPQ=0
!        FIND NU IN INTERVAL [-.5,.5) if ID=2  ( CALCULATION OF Q )
  NU=MOD(NU1,1.D0)
  if ( NU >= .5D0) NU=NU-1.D0
!        FIND NU IN INTERVAL (-1.5,-.5] if ID=1,3, OR 4  ( CALC. OF P )
  if ( ID /= 2.AND.NU > -.5D0) NU=NU-1.D0
!        CALCULATE MU FACTORIAL
  K=MU
  DMU=MU
  if ( MU <= 0) go to 60
  FACTMU=1.D0
  IF=0
  DO 50 I=1,K
  FACTMU=FACTMU*I
   50 call DXADJ(FACTMU,IF,IERROR)
  if (IERROR /= 0) RETURN
   60 if ( K == 0) FACTMU=1.D0
  if ( K == 0) IF=0
!
!        X=COS(THETA)
!        Y=SIN(THETA/2)**2=(1-X)/2=.5-.5*X
!        R=TAN(THETA/2)=SQRT((1-X)/(1+X)
!
  X=COS(THETA)
  Y=SIN(THETA/2.D0)**2
  R=TAN(THETA/2.D0)
!
!        USE ASCENDING SERIES TO CALCULATE TWO VALUES OF P OR Q
!        FOR USE AS STARTING VALUES IN RECURRENCE RELATION.
!
  PQ2=0.0D0
  DO 100 J=1,2
  IPQ1=0
  if ( ID == 2) go to 80
!
!        SERIES FOR P ( ID = 1, 3, OR 4 )
!        P(-MU,NU,X)=1./FACTORIAL(MU)*SQRT(((1.-X)/(1.+X))**MU)
!                *SUM(FROM 0 TO J0-1)A(J)*(.5-.5*X)**J
!
  IPQ=0
  PQ=1.D0
  A=1.D0
  IA=0
  DO 65 I=2,J0
  DI=I
  A=A*Y*(DI-2.D0-NU)*(DI-1.D0+NU)/((DI-1.D0+DMU)*(DI-1.D0))
  call DXADJ(A,IA,IERROR)
  if (IERROR /= 0) RETURN
  if ( A == 0.D0) go to 66
  call DXADD(PQ,IPQ,A,IA,PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
   65 CONTINUE
   66 CONTINUE
  if ( MU <= 0) go to 90
  X2=R
  X1=PQ
  K=MU
  DO 77 I=1,K
  X1=X1*X2
   77 call DXADJ(X1,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  PQ=X1/FACTMU
  IPQ=IPQ-IF
  call DXADJ(PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  go to 90
!
!        Z=-LN(R)=.5*LN((1+X)/(1-X))
!
   80 Z=-LOG(R)
  W=DXPSI(NU+1.D0,IPSIK,IPSIX)
  XS=1.D0/SIN(THETA)
!
!        SERIES SUMMATION FOR Q ( ID = 2 )
!        Q(0,NU,X)=SUM(FROM 0 TO J0-1)((.5*LN((1+X)/(1-X))
!    +DXPSI(J+1,IPSIK,IPSIX)-DXPSI(NU+1,IPSIK,IPSIX)))*A(J)*(.5-.5*X)**J
!
!        Q(1,NU,X)=-SQRT(1./(1.-X**2))+SQRT((1-X)/(1+X))
!             *SUM(FROM 0 T0 J0-1)(-NU*(NU+1)/2*LN((1+X)/(1-X))
!                 +(J-NU)*(J+NU+1)/(2*(J+1))+NU*(NU+1)*
!     (DXPSI(NU+1,IPSIK,IPSIX)-DXPSI(J+1,IPSIK,IPSIX))*A(J)*(.5-.5*X)**J
!
!        NOTE, IN THIS LOOP K=J+1
!
  PQ=0.D0
  IPQ=0
  IA=0
  A=1.D0
  DO 85 K=1,J0
  FLOK=K
  if ( K == 1) go to 81
  A=A*Y*(FLOK-2.D0-NU)*(FLOK-1.D0+NU)/((FLOK-1.D0+DMU)*(FLOK-1.D0))
  call DXADJ(A,IA,IERROR)
  if (IERROR /= 0) RETURN
   81 CONTINUE
  if ( MU >= 1) go to 83
  X1=(DXPSI(FLOK,IPSIK,IPSIX)-W+Z)*A
  IX1=IA
  call DXADD(PQ,IPQ,X1,IX1,PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  go to 85
   83 X1=(NU*(NU+1.D0)*(Z-W+DXPSI(FLOK,IPSIK,IPSIX))+(NU-FLOK+1.D0) &
    *(NU+FLOK)/(2.D0*FLOK))*A
  IX1=IA
  call DXADD(PQ,IPQ,X1,IX1,PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
   85 CONTINUE
  if ( MU >= 1) PQ=-R*PQ
  IXS=0
  if ( MU >= 1) call DXADD(PQ,IPQ,-XS,IXS,PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  if ( J == 2) MU=-MU
  if ( J == 2) DMU=-DMU
   90 if ( J == 1) PQ2=PQ
  if ( J == 1) IPQ2=IPQ
  NU=NU+1.D0
  100 CONTINUE
  K=0
  if ( NU-1.5D0 < NU1) go to 120
  K=K+1
  PQA(K)=PQ2
  IPQA(K)=IPQ2
  if ( NU > NU2+.5D0) RETURN
  120 PQ1=PQ
  IPQ1=IPQ
  if ( NU < NU1+.5D0) go to 130
  K=K+1
  PQA(K)=PQ
  IPQA(K)=IPQ
  if ( NU > NU2+.5D0) RETURN
!
!        FORWARD NU-WISE RECURRENCE FOR F(MU,NU,X) FOR FIXED MU
!        USING
!        (NU+MU+1)*F(MU,NU,X)=(2.*NU+1)*F(MU,NU,X)-(NU-MU)*F(MU,NU-1,X)
!        WHERE F(MU,NU,X) MAY BE P(-MU,NU,X) OR if MU IS REPLACED
!        BY -MU THEN F(MU,NU,X) MAY BE Q(MU,NU,X).
!        NOTE, IN THIS LOOP, NU=NU+1
!
  130 X1=(2.D0*NU-1.D0)/(NU+DMU)*X*PQ1
  X2=(NU-1.D0-DMU)/(NU+DMU)*PQ2
  call DXADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  call DXADJ(PQ,IPQ,IERROR)
  if (IERROR /= 0) RETURN
  NU=NU+1.D0
  PQ2=PQ1
  IPQ2=IPQ1
  go to 120
!
end
