subroutine POLYVL (NDER, XX, YFIT, YP, N, X, C, WORK, IERR)
!
!! POLYVL calculates the value of a polynomial and its first NDER ...
!  derivatives where the polynomial was produced by a previous call to POLINT.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3
!***TYPE      SINGLE PRECISION (POLYVL-S, DPOLVL-D)
!***KEYWORDS  POLYNOMIAL EVALUATION
!***AUTHOR  Huddleston, R. E., (SNLL)
!***DESCRIPTION
!
!     Written by Robert E. Huddleston, Sandia Laboratories, Livermore
!
!     Abstract -
!        Subroutine POLYVL calculates the value of the polynomial and
!     its first NDER derivatives where the polynomial was produced by
!     a previous call to POLINT.
!        The variable N and the arrays X and C must not be altered
!     between the call to POLINT and the call to POLYVL.
!
!     ******  Dimensioning Information *******
!
!     YP   must be dimensioned by at least NDER
!     X    must be dimensioned by at least N (see the abstract )
!     C    must be dimensioned by at least N (see the abstract )
!     WORK must be dimensioned by at least 2*N if NDER is  >  0.
!
!     *** Note ***
!       If NDER=0, neither YP nor WORK need to be dimensioned variables.
!       If NDER=1, YP does not need to be a dimensioned variable.
!
!
!     *****  Input parameters
!
!     NDER - the number of derivatives to be evaluated
!
!     XX   - the argument at which the polynomial and its derivatives
!            are to be evaluated.
!
!     N    - *****
!            *       N, X, and C must not be altered between the call
!     X    - *       to POLINT and the call to POLYVL.
!     C    - *****
!
!
!     *****  Output Parameters
!
!     YFIT - the value of the polynomial at XX
!
!     YP   - the derivatives of the polynomial at XX.  The derivative of
!            order J at XX is stored in  YP(J) , J = 1,...,NDER.
!
!     IERR - Output error flag with the following possible values.
!          = 1  indicates normal execution
!
!     ***** Storage Parameters
!
!     WORK  = this is an array to provide internal working storage for
!             POLYVL.  It must be dimensioned by at least 2*N if NDER is
!              >  0.  If NDER=0, WORK does not need to be a dimensioned
!             variable.
!
!***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!                 Curve fitting by polynomials in one variable, Report
!                 SLA-74-0270, Sandia Laboratories, June 1974.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   740601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  POLYVL
  DIMENSION  YP(*),X(*),C(*),WORK(*)
!***FIRST EXECUTABLE STATEMENT  POLYVL
  IERR=1
     if (NDER > 0) go to 10020
!
!     *****   CODING FOR THE CASE NDER = 0
!
  PIONE=1.0
  PONE=C(1)
  YFIT=PONE
  if (N == 1) RETURN
  DO 10010 K=2,N
  PITWO=(XX-X(K-1))*PIONE
  PIONE=PITWO
  PTWO=PONE+PITWO*C(K)
  PONE=PTWO
10010 CONTINUE
  YFIT=PTWO
  return
!
!     *****   END OF NDER = 0 CASE
!
10020 CONTINUE
     if (N > 1) go to 10040
  YFIT=C(1)
!
!     *****  CODING FOR THE CASE  N=1 AND NDER  >  0
!
  DO 10030 K=1,NDER
  YP(K)=0.0
10030 CONTINUE
  return
!
!     *****  END OF THE CASE  N = 1 AND  NDER  >  0
!
10040 CONTINUE
     if (NDER < N) go to 10050
!
!     *****  SET FLAGS FOR NUMBER OF DERIVATIVES AND FOR DERIVATIVES
!            IN EXCESS OF THE DEGREE (N-1) OF THE POLYNOMIAL.
!
  IZERO=1
  NDR=N-1
     go to 10060
10050 CONTINUE
  IZERO=0
  NDR=NDER
10060 CONTINUE
  M=NDR+1
  MM=M
!
!     *****  START OF THE CASE NDER  >  0  AND N  >  1
!     *****  THE POLYNOMIAL AND ITS DERIVATIVES WILL BE EVALUATED AT XX
!
  DO 10070 K=1,NDR
  YP(K)=C(K+1)
10070 CONTINUE
!
!     *****  THE FOLLOWING SECTION OF CODE IS EASIER TO READ if ONE
!            BREAKS WORK INTO TWO ARRAYS W AND V. THE CODE WOULD THEN
!            READ
!                W(1) = 1.
!                PONE = C(1)
!               *DO   K = 2,N
!               *   V(K-1) =  XX - X(K-1)
!               *   W(K)   =  V(K-1)*W(K-1)
!               *   PTWO   =  PONE + W(K)*C(K)
!               *   PONE   =  PWO
!
!               YFIT = PTWO
!
  WORK(1)=1.0
  PONE=C(1)
  DO 10080 K=2,N
  KM1=K-1
  NPKM1=N+K-1
  WORK(NPKM1)=XX-X(KM1)
  WORK(K)=WORK(NPKM1)*WORK(KM1)
  PTWO=PONE+WORK(K)*C(K)
  PONE=PTWO
10080 CONTINUE
  YFIT=PTWO
!
!     ** AT THIS POINT THE POLYNOMIAL HAS BEEN EVALUATED AND INFORMATION
!        FOR THE DERIVATIVE EVALUATIONS HAVE BEEN STORED IN THE ARRAY
!        WORK
     if (N == 2) go to 10110
  if (M == N) MM=NDR
!
!     ***** EVALUATE THE DERIVATIVES AT XX
!
!                  ******  DO K=2,MM   (FOR MOST CASES, MM = NDER + 1)
!                  *  ******  DO I=2,N-K+1
!                  *  *       W(I) = V(K-2+I)*W(I-1) + W(I)
!                  *  *       YP(K-1) = YP(K-1) + W(I)*C(K-1+I)
!                  ******  CONTINUE
!
  DO 10090 K=2,MM
  NMKP1=N-K+1
  KM1=K-1
  KM2PN=K-2+N
  DO 10090 I=2,NMKP1
  KM2PNI=KM2PN+I
  IM1=I-1
  KM1PI=KM1+I
  WORK(I)=WORK(KM2PNI)*WORK(IM1)+WORK(I)
  YP(KM1)=YP(KM1)+WORK(I)*C(KM1PI)
10090 CONTINUE
     if (NDR == 1) go to 10110
  FAC=1.0
  DO 10100 K=2,NDR
  XK=K
  FAC=XK*FAC
  YP(K)=FAC*YP(K)
10100 CONTINUE
!
!     ***** END OF DERIVATIVE EVALUATIONS
!
10110 CONTINUE
  if (IZERO == 0) RETURN
!
!     *****  SET EXCESS DERIVATIVES TO ZERO.
!
  DO 10120 K=N,NDER
  YP(K)=0.0
10120 CONTINUE
  return
end
