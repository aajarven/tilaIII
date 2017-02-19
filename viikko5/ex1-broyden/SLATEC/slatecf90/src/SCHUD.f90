subroutine SCHUD (R, LDR, P, X, Z, LDZ, NZ, Y, RHO, C, S)
!
!! SCHUD updates an augmented Cholesky decomposition of the triangular part ...
!  of an augmented QR decomposition.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D7B
!***TYPE      SINGLE PRECISION (SCHUD-S, DCHUD-D, CCHUD-C)
!***KEYWORDS  CHOLESKY DECOMPOSITION, LINEAR ALGEBRA, LINPACK, MATRIX,
!             UPDATE
!***AUTHOR  Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     SCHUD updates an augmented Cholesky decomposition of the
!     triangular part of an augmented QR decomposition.  Specifically,
!     given an upper triangular matrix R of order P, a row vector
!     X, a column vector Z, and a scalar Y, SCHUD determines a
!     unitary matrix U and a scalar ZETA such that
!
!
!                              (R  Z)     (RR   ZZ )
!                         U  * (    )  =  (        ) ,
!                              (X  Y)     ( 0  ZETA)
!
!     where RR is upper triangular.  If R and Z have been
!     obtained from the factorization of a least squares
!     problem, then RR and ZZ are the factors corresponding to
!     the problem with the observation (X,Y) appended.  In this
!     case, if RHO is the norm of the residual vector, then the
!     norm of the residual vector of the updated problem is
!     SQRT(RHO**2 + ZETA**2).  SCHUD will simultaneously update
!     several triplets (Z,Y,RHO).
!     For a less terse description of what SCHUD does and how
!     it may be applied, see the LINPACK guide.
!
!     The matrix U is determined as the product U(P)*...*U(1),
!     where U(I) is a rotation in the (I,P+1) plane of the
!     form
!
!                       (     C(I)      S(I) )
!                       (                    ) .
!                       (    -S(I)      C(I) )
!
!     The rotations are chosen so that C(I) is real.
!
!     On Entry
!
!         R      REAL(LDR,P), where LDR  >=  P.
!                R contains the upper triangular matrix
!                that is to be updated.  The part of R
!                below the diagonal is not referenced.
!
!         LDR    INTEGER.
!                LDR is the leading dimension of the array R.
!
!         P      INTEGER.
!                P is the order of the matrix R.
!
!         X      REAL(P).
!                X contains the row to be added to R.  X is
!                not altered by SCHUD.
!
!         Z      REAL(LDZ,NZ), where LDZ  >=  P.
!                Z is an array containing NZ P-vectors to
!                be updated with R.
!
!         LDZ    INTEGER.
!                LDZ is the leading dimension of the array Z.
!
!         NZ     INTEGER.
!                NZ is the number of vectors to be updated.
!                NZ may be zero, in which case Z, Y, and RHO
!                are not referenced.
!
!         Y      REAL(NZ).
!                Y contains the scalars for updating the vectors
!                Z.  Y is not altered by SCHUD.
!
!         RHO    REAL(NZ).
!                RHO contains the norms of the residual
!                vectors that are to be updated.  If RHO(J)
!                is negative, it is left unaltered.
!
!     On Return
!
!         RC
!         RHO    contain the updated quantities.
!         Z
!
!         C      REAL(P).
!                C contains the cosines of the transforming
!                rotations.
!
!         S      REAL(P).
!                S contains the sines of the transforming
!                rotations.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SROTG
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SCHUD
  INTEGER LDR,P,LDZ,NZ
  REAL RHO(*),C(*)
  REAL R(LDR,*),X(*),Z(LDZ,*),Y(*),S(*)
!
  INTEGER I,J,JM1
  REAL AZETA,SCALE
  REAL T,XJ,ZETA
!
!     UPDATE R.
!
!***FIRST EXECUTABLE STATEMENT  SCHUD
  DO 30 J = 1, P
     XJ = X(J)
!
!        APPLY THE PREVIOUS ROTATIONS.
!
     JM1 = J - 1
     if (JM1  <  1) go to 20
     DO 10 I = 1, JM1
        T = C(I)*R(I,J) + S(I)*XJ
        XJ = C(I)*XJ - S(I)*R(I,J)
        R(I,J) = T
   10    CONTINUE
   20    CONTINUE
!
!        COMPUTE THE NEXT ROTATION.
!
     call SROTG(R(J,J),XJ,C(J),S(J))
   30 CONTINUE
!
!     if REQUIRED, UPDATE Z AND RHO.
!
  if (NZ  <  1) go to 70
  DO 60 J = 1, NZ
     ZETA = Y(J)
     DO 40 I = 1, P
        T = C(I)*Z(I,J) + S(I)*ZETA
        ZETA = C(I)*ZETA - S(I)*Z(I,J)
        Z(I,J) = T
   40    CONTINUE
     AZETA = ABS(ZETA)
     if (AZETA  ==  0.0E0 .OR. RHO(J)  <  0.0E0) go to 50
        SCALE = AZETA + RHO(J)
        RHO(J) = SCALE*SQRT((AZETA/SCALE)**2+(RHO(J)/SCALE)**2)
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
  return
end
