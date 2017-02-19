subroutine POLCOF (XX, N, X, C, D, WORK)
!
!! POLCOF computes the coefficients of the polynomial fit (including ...
!            Hermite polynomial fits) produced by a previous call to
!            POLINT.
!
!***LIBRARY   SLATEC
!***CATEGORY  E1B
!***TYPE      SINGLE PRECISION (POLCOF-S, DPOLCF-D)
!***KEYWORDS  COEFFICIENTS, POLYNOMIAL
!***AUTHOR  Huddleston, R. E., (SNLL)
!***DESCRIPTION
!
!     Written by Robert E. Huddleston, Sandia Laboratories, Livermore
!
!     Abstract
!        Subroutine POLCOF computes the coefficients of the polynomial
!     fit (including Hermite polynomial fits ) produced by a previous
!     call to POLINT. The coefficients of the polynomial, expanded about
!     XX, are stored in the array D. The expansion is of the form
!     P(Z) = D(1) + D(2)*(Z-XX) +D(3)*((Z-XX)**2) + ... +
!                                                  D(N)*((Z-XX)**(N-1)).
!     Between the call to POLINT and the call to POLCOF the variable N
!     and the arrays X and C must not be altered.
!
!     *****  INPUT PARAMETERS
!
!     XX   - The point about which the Taylor expansion is to be made.
!
!     N    - ****
!            *     N, X, and C must remain unchanged between the
!     X    - *     call to POLINT or the call to POLCOF.
!     C    - ****
!
!     *****  OUTPUT PARAMETER
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
!     by POLINT. You may call POLYVL to perform the task, or you may
!     call POLCOF to obtain the coefficients of the Taylor expansion and
!     then write your own evaluation scheme. Due to the inherent errors
!     in the computations of the Taylor expansion from the Newton
!     coefficients produced by POLINT, much more accuracy may be
!     expected by calling POLYVL as opposed to writing your own scheme.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   890213  DATE WRITTEN
!   891024  Corrected KEYWORD section.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  POLCOF
!
  DIMENSION X(*), C(*), D(*), WORK(*)
!***FIRST EXECUTABLE STATEMENT  POLCOF
  DO 10010 K=1,N
  D(K)=C(K)
10010 CONTINUE
  if (N == 1) RETURN
  WORK(1)=1.0
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
