subroutine CCHDD (R, LDR, P, X, Z, LDZ, NZ, Y, RHO, C, S, INFO)
!
!! CCHDD downdates an augmented Cholesky decomposition or the ...
!            triangular factor of an augmented QR decomposition.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D7B
!***TYPE      COMPLEX (SCHDD-S, DCHDD-D, CCHDD-C)
!***KEYWORDS  CHOLESKY DECOMPOSITION, DOWNDATE, LINEAR ALGEBRA, LINPACK,
!             MATRIX
!***AUTHOR  Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     CCHDD downdates an augmented Cholesky decomposition or the
!     triangular factor of an augmented QR decomposition.
!     Specifically, given an upper triangular matrix R of order P,  a
!     row vector X, a column vector Z, and a scalar Y, CCHDD
!     determines a unitary matrix U and a scalar ZETA such that
!
!                        (R   Z )     (RR  ZZ)
!                    U * (      )  =  (      ) ,
!                        (0 ZETA)     ( X   Y)
!
!     where RR is upper triangular.  If R and Z have been obtained
!     from the factorization of a least squares problem, then
!     RR and ZZ are the factors corresponding to the problem
!     with the observation (X,Y) removed.  In this case, if RHO
!     is the norm of the residual vector, then the norm of
!     the residual vector of the downdated problem is
!     SQRT(RHO**2 - ZETA**2).  CCHDD will simultaneously downdate
!     several triplets (Z,Y,RHO) along with R.
!     For a less terse description of what CCHDD does and how
!     it may be applied, see the LINPACK Guide.
!
!     The matrix U is determined as the product U(1)*...*U(P)
!     where U(I) is a rotation in the (P+1,I)-plane of the
!     form
!
!                       ( C(I)  -CONJG(S(I)) )
!                       (                    ) .
!                       ( S(I)       C(I)    )
!
!     the rotations are chosen so that C(I) is real.
!
!     The user is warned that a given downdating problem may
!     be impossible to accomplish or may produce
!     inaccurate results.  For example, this can happen
!     if X is near a vector whose removal will reduce the
!     rank of R.  Beware.
!
!     On Entry
!
!         R      COMPLEX(LDR,P), where LDR  >=  P.
!                R contains the upper triangular matrix
!                that is to be downdated.  The part of R
!                below the diagonal is not referenced.
!
!         LDR    INTEGER.
!                LDR is the leading dimension of the array R.
!
!         p      INTEGER.
!                P is the order of the matrix R.
!
!         X      COMPLEX(P).
!                X contains the row vector that is to
!                be removed from R.  X is not altered by CCHDD.
!
!         Z      COMPLEX(LDZ,NZ), where LDZ  >=  P.
!                Z is an array of NZ P-vectors which
!                are to be downdated along with R.
!
!         LDZ    INTEGER.
!                LDZ is the leading dimension of the array Z.
!
!         NZ     INTEGER.
!                NZ is the number of vectors to be downdated
!                NZ may be zero, in which case Z, Y, and RHO
!                are not referenced.
!
!         Y      COMPLEX(NZ).
!                Y contains the scalars for the downdating
!                of the vectors Z.  Y is not altered by CCHDD.
!
!         RHO    REAL(NZ).
!                RHO contains the norms of the residual
!                vectors that are to be downdated.
!
!     On Return
!
!         R
!         Z      contain the downdated quantities.
!         RHO
!
!         C      REAL(P).
!                C contains the cosines of the transforming
!                rotations.
!
!         S      COMPLEX(P).
!                S contains the sines of the transforming
!                rotations.
!
!         INFO   INTEGER.
!                INFO is set as follows.
!
!                   INFO = 0  if the entire downdating
!                             was successful.
!
!                   INFO =-1  if R could not be downdated.
!                             in this case, all quantities
!                             are left unaltered.
!
!                   INFO = 1  if some RHO could not be
!                             downdated.  The offending RHO's are
!                             set to -1.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CDOTC, SCNRM2
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CCHDD
  INTEGER LDR,P,LDZ,NZ,INFO
  COMPLEX R(LDR,*),X(*),Z(LDZ,*),Y(*),S(*)
  REAL RHO(*),C(*)
!
  INTEGER I,II,J
  REAL A,ALPHA,AZETA,NORM,SCNRM2
  COMPLEX CDOTC,T,ZETA,B,XX
!
!     SOLVE THE SYSTEM CTRANS(R)*A = X, PLACING THE RESULT
!     IN THE ARRAY S.
!
!***FIRST EXECUTABLE STATEMENT  CCHDD
  INFO = 0
  S(1) = CONJG(X(1))/CONJG(R(1,1))
  if (P  <  2) go to 20
  DO 10 J = 2, P
     S(J) = CONJG(X(J)) - CDOTC(J-1,R(1,J),1,S,1)
     S(J) = S(J)/CONJG(R(J,J))
   10 CONTINUE
   20 CONTINUE
  NORM = SCNRM2(P,S,1)
  if (NORM  <  1.0E0) go to 30
     INFO = -1
  go to 120
   30 CONTINUE
     ALPHA = SQRT(1.0E0-NORM**2)
!
!        DETERMINE THE TRANSFORMATIONS.
!
     DO 40 II = 1, P
        I = P - II + 1
        SCALE = ALPHA + ABS(S(I))
        A = ALPHA/SCALE
        B = S(I)/SCALE
        NORM = SQRT(A**2+REAL(B)**2+AIMAG(B)**2)
        C(I) = A/NORM
        S(I) = CONJG(B)/NORM
        ALPHA = SCALE*NORM
   40    CONTINUE
!
!        APPLY THE TRANSFORMATIONS TO R.
!
     DO 60 J = 1, P
        XX = (0.0E0,0.0E0)
        DO 50 II = 1, J
           I = J - II + 1
           T = C(I)*XX + S(I)*R(I,J)
           R(I,J) = C(I)*R(I,J) - CONJG(S(I))*XX
           XX = T
   50       CONTINUE
   60    CONTINUE
!
!        if REQUIRED, DOWNDATE Z AND RHO.
!
     if (NZ  <  1) go to 110
     DO 100 J = 1, NZ
        ZETA = Y(J)
        DO 70 I = 1, P
           Z(I,J) = (Z(I,J) - CONJG(S(I))*ZETA)/C(I)
           ZETA = C(I)*ZETA - S(I)*Z(I,J)
   70       CONTINUE
        AZETA = ABS(ZETA)
        if (AZETA  <=  RHO(J)) go to 80
           INFO = 1
           RHO(J) = -1.0E0
        go to 90
   80       CONTINUE
           RHO(J) = RHO(J)*SQRT(1.0E0-(AZETA/RHO(J))**2)
   90       CONTINUE
  100    CONTINUE
  110    CONTINUE
  120 CONTINUE
  return
end
