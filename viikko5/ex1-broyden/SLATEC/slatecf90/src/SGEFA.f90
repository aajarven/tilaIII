subroutine SGEFA (A, LDA, N, IPVT, INFO)
!
!! SGEFA factors a matrix using Gaussian elimination.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2A1
!***TYPE      SINGLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
!***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SGEFA factors a real matrix by Gaussian elimination.
!
!     SGEFA is usually called by SGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for SGECO) = (1 + 9/N)*(Time for SGEFA) .
!
!     On Entry
!
!        A       REAL(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U , where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K)  ==  0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that SGESL or SGEDI will divide by zero
!                     if called.  Use  RCOND  in SGECO for a reliable
!                     indication of singularity.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  ISAMAX, SAXPY, SSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SGEFA
  INTEGER LDA,N,IPVT(*),INFO
  REAL A(LDA,*)
!
  REAL T
  INTEGER ISAMAX,J,K,KP1,L,NM1
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
!***FIRST EXECUTABLE STATEMENT  SGEFA
  INFO = 0
  NM1 = N - 1
  if (NM1  <  1) go to 70
  DO 60 K = 1, NM1
     KP1 = K + 1
!
!        FIND L = PIVOT INDEX
!
     L = ISAMAX(N-K+1,A(K,K),1) + K - 1
     IPVT(K) = L
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
     if (A(L,K)  ==  0.0E0) go to 40
!
!           INTERCHANGE if NECESSARY
!
        if (L  ==  K) go to 10
           T = A(L,K)
           A(L,K) = A(K,K)
           A(K,K) = T
   10       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
        T = -1.0E0/A(K,K)
        call SSCAL(N-K,T,A(K+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
        DO 30 J = KP1, N
           T = A(L,J)
           if (L  ==  K) go to 20
              A(L,J) = A(K,J)
              A(K,J) = T
   20          CONTINUE
           call SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
     go to 50
   40    CONTINUE
        INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
  IPVT(N) = N
  if (A(N,N)  ==  0.0E0) INFO = N
  return
end
