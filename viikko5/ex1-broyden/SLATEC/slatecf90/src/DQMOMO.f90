subroutine DQMOMO (ALFA, BETA, RI, RJ, RG, RH, INTEGR)
!
!! DQMOMO computes modified Chebyshev moments.
!  The K-th
!            modified Chebyshev moment is defined as the integral over
!            (-1,1) of W(X)*T(K,X), where T(K,X) is the Chebyshev
!            polynomial of degree K.
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A1, C3A2
!***TYPE      DOUBLE PRECISION (QMOMO-S, DQMOMO-D)
!***KEYWORDS  MODIFIED CHEBYSHEV MOMENTS, QUADPACK, QUADRATURE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        MODIFIED CHEBYSHEV MOMENTS
!        STANDARD FORTRAN SUBROUTINE
!        DOUBLE PRECISION VERSION
!
!        PARAMETERS
!           ALFA   - Double precision
!                    Parameter in the weight function W(X), ALFA > (-1)
!
!           BETA   - Double precision
!                    Parameter in the weight function W(X), BETA > (-1)
!
!           RI     - Double precision
!                    Vector of dimension 25
!                    RI(K) is the integral over (-1,1) of
!                    (1+X)**ALFA*T(K-1,X), K = 1, ..., 25.
!
!           RJ     - Double precision
!                    Vector of dimension 25
!                    RJ(K) is the integral over (-1,1) of
!                    (1-X)**BETA*T(K-1,X), K = 1, ..., 25.
!
!           RG     - Double precision
!                    Vector of dimension 25
!                    RG(K) is the integral over (-1,1) of
!                    (1+X)**ALFA*LOG((1+X)/2)*T(K-1,X), K = 1, ..., 25.
!
!           RH     - Double precision
!                    Vector of dimension 25
!                    RH(K) is the integral over (-1,1) of
!                    (1-X)**BETA*LOG((1-X)/2)*T(K-1,X), K = 1, ..., 25.
!
!           INTEGR - Integer
!                    Input parameter indicating the modified
!                    Moments to be computed
!                    INTEGR = 1 compute RI, RJ
!                           = 2 compute RI, RJ, RG
!                           = 3 compute RI, RJ, RH
!                           = 4 compute RI, RJ, RG, RH
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820101  DATE WRITTEN
!   891009  Removed unreferenced statement label.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQMOMO
!
  DOUBLE PRECISION ALFA,ALFP1,ALFP2,AN,ANM1,BETA,BETP1,BETP2,RALF, &
    RBET,RG,RH,RI,RJ
  INTEGER I,IM1,INTEGR
!
  DIMENSION RG(25),RH(25),RI(25),RJ(25)
!
!
!***FIRST EXECUTABLE STATEMENT  DQMOMO
  ALFP1 = ALFA+0.1D+01
  BETP1 = BETA+0.1D+01
  ALFP2 = ALFA+0.2D+01
  BETP2 = BETA+0.2D+01
  RALF = 0.2D+01**ALFP1
  RBET = 0.2D+01**BETP1
!
!           COMPUTE RI, RJ USING A FORWARD RECURRENCE RELATION.
!
  RI(1) = RALF/ALFP1
  RJ(1) = RBET/BETP1
  RI(2) = RI(1)*ALFA/ALFP2
  RJ(2) = RJ(1)*BETA/BETP2
  AN = 0.2D+01
  ANM1 = 0.1D+01
  DO 20 I=3,25
    RI(I) = -(RALF+AN*(AN-ALFP2)*RI(I-1))/(ANM1*(AN+ALFP1))
    RJ(I) = -(RBET+AN*(AN-BETP2)*RJ(I-1))/(ANM1*(AN+BETP1))
    ANM1 = AN
    AN = AN+0.1D+01
   20 CONTINUE
  if ( INTEGR == 1) go to 70
  if ( INTEGR == 3) go to 40
!
!           COMPUTE RG USING A FORWARD RECURRENCE RELATION.
!
  RG(1) = -RI(1)/ALFP1
  RG(2) = -(RALF+RALF)/(ALFP2*ALFP2)-RG(1)
  AN = 0.2D+01
  ANM1 = 0.1D+01
  IM1 = 2
  DO 30 I=3,25
    RG(I) = -(AN*(AN-ALFP2)*RG(IM1)-AN*RI(IM1)+ANM1*RI(I))/ &
    (ANM1*(AN+ALFP1))
    ANM1 = AN
    AN = AN+0.1D+01
    IM1 = I
   30 CONTINUE
  if ( INTEGR == 2) go to 70
!
!           COMPUTE RH USING A FORWARD RECURRENCE RELATION.
!
   40 RH(1) = -RJ(1)/BETP1
  RH(2) = -(RBET+RBET)/(BETP2*BETP2)-RH(1)
  AN = 0.2D+01
  ANM1 = 0.1D+01
  IM1 = 2
  DO 50 I=3,25
    RH(I) = -(AN*(AN-BETP2)*RH(IM1)-AN*RJ(IM1)+ &
    ANM1*RJ(I))/(ANM1*(AN+BETP1))
    ANM1 = AN
    AN = AN+0.1D+01
    IM1 = I
   50 CONTINUE
  DO 60 I=2,25,2
    RH(I) = -RH(I)
   60 CONTINUE
   70 DO 80 I=2,25,2
    RJ(I) = -RJ(I)
   80 CONTINUE
  return
end
