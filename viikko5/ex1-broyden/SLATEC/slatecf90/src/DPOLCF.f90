subroutine DPOLCF (XX, N, X, C, D, WORK)
!
!! DPOLCF computes the coefficients of the polynomial fit (including ...
!  Hermite polynomial fits) produced by a previous call to POLINT.
!
!***LIBRARY   SLATEC
!***CATEGORY  E1B
!***TYPE      DOUBLE PRECISION (POLCOF-S, DPOLCF-D)
!***KEYWORDS  COEFFICIENTS, POLYNOMIAL
!***AUTHOR  Huddleston, R. E., (SNLL)
!***DESCRIPTION
!
!     Abstract
!        Subroutine DPOLCF computes the coefficients of the polynomial
!     fit (including Hermite polynomial fits ) produced by a previous
!     call to DPLINT.  The coefficients of the polynomial, expanded
!     about XX, are stored in the array D. The expansion is of the form
!     P(Z) = D(1) + D(2)*(Z-XX) +D(3)*((Z-XX)**2) + ... +
!                                                  D(N)*((Z-XX)**(N-1)).
!     Between the call to DPLINT and the call to DPOLCF the variable N
!     and the arrays X and C must not be altered.
!
!     *****  INPUT PARAMETERS
!      *** All TYPE REAL variables are DOUBLE PRECISION ***
!
!     XX   - The point about which the Taylor expansion is to be made.
!
!     N    - ****
!            *     N, X, and C must remain unchanged between the
!     X    - *     call to DPLINT and the call to DPOLCF.
!     C    - ****
!
!     *****  OUTPUT PARAMETER
!      *** All TYPE REAL variables are DOUBLE PRECISION ***
!
!     D    - The array of coefficients for the Taylor expansion as
!            explained in the abstract
!
!     *****  STORAGE PARAMETER
!
!     WORK - This is an array to provide internal working storage. It
!            must be dimensioned by at least 2*N in the calling program.
!
!
!     **** Note - There are two methods for evaluating the fit produced
!     by DPLINT. You may call DPOLVL to perform the task, or you may
!     call DPOLCF to obtain the coefficients of the Taylor expansion and
!     then write your own evaluation scheme. Due to the inherent errors
!     in the computations of the Taylor expansion from the Newton
!     coefficients produced by DPLINT, much more accuracy may be
!     expected by calling DPOLVL as opposed to writing your own scheme.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   890213  DATE WRITTEN
!   891006  Cosmetic changes to prologue.  (WRB)
!   891024  Corrected KEYWORD section.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DPOLCF
!
  INTEGER I,IM1,K,KM1,KM1PI,KM2N,KM2NPI,N,NM1,NMKP1,NPKM1
  DOUBLE PRECISION C(*),D(*),PONE,PTWO,X(*),XX,WORK(*)
!***FIRST EXECUTABLE STATEMENT  DPOLCF
  DO 10010 K=1,N
  D(K)=C(K)
10010 CONTINUE
  if (N == 1) RETURN
  WORK(1)=1.0D0
  PONE=C(1)
  NM1=N-1
  DO 10020 K=2,N
  KM1=K-1
  NPKM1=N+K-1
  WORK(NPKM1)=XX-X(KM1)
  WORK(K)=WORK(NPKM1)*WORK(KM1)
  PTWO=PONE+WORK(K)*C(K)
  PONE=PTWO
10020 CONTINUE
  D(1)=PTWO
  if (N == 2) RETURN
  DO 10030 K=2,NM1
  KM1=K-1
  KM2N=K-2+N
  NMKP1=N-K+1
  DO 10030 I=2,NMKP1
  KM2NPI=KM2N+I
  IM1=I-1
  KM1PI=KM1+I
  WORK(I)=WORK(KM2NPI)*WORK(IM1)+WORK(I)
  D(K)=D(K)+WORK(I)*D(KM1PI)
10030 CONTINUE
  return
end
