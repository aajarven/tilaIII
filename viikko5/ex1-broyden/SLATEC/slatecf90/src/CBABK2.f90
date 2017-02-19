subroutine CBABK2 (NM, N, LOW, IGH, SCALE, M, ZR, ZI)
!
!! CBABK2 forms the eigenvectors of a complex general matrix from the ...
!            eigenvectors of matrix output from CBAL.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      COMPLEX (BALBAK-S, CBABK2-C)
!***KEYWORDS  EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure
!     CBABK2, which is a complex version of BALBAK,
!     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
!
!     This subroutine forms the eigenvectors of a COMPLEX GENERAL
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  CBAL.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, ZR and ZI, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix Z=(ZR,ZI).  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        LOW and IGH are INTEGER variables determined by  CBAL.
!
!        SCALE contains information determining the permutations and
!          scaling factors used by  CBAL.  SCALE is a one-dimensional
!          REAL array, dimensioned SCALE(N).
!
!        M is the number of eigenvectors to be back transformed.
!          M is an INTEGER variable.
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the eigenvectors to be back transformed in their first
!          M columns.  ZR and ZI are two-dimensional REAL arrays,
!          dimensioned ZR(NM,M) and ZI(NM,M).
!
!     On OUTPUT
!
!        ZR and ZI contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first M columns.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CBABK2
!
  INTEGER I,J,K,M,N,II,NM,IGH,LOW
  REAL SCALE(*),ZR(NM,*),ZI(NM,*)
  REAL S
!
!***FIRST EXECUTABLE STATEMENT  CBABK2
  if (M  ==  0) go to 200
  if (IGH  ==  LOW) go to 120
!
  DO 110 I = LOW, IGH
     S = SCALE(I)
!     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
!                if THE FOREGOING STATEMENT IS REPLACED BY
!                S=1.0E0/SCALE(I). ..........
     DO 100 J = 1, M
        ZR(I,J) = ZR(I,J) * S
        ZI(I,J) = ZI(I,J) * S
  100    CONTINUE
!
  110 CONTINUE
!     .......... FOR I=LOW-1 STEP -1 UNTIL 1,
!                IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
     I = II
     if (I  >=  LOW .AND. I  <=  IGH) go to 140
     if (I  <  LOW) I = LOW - II
     K = SCALE(I)
     if (K  ==  I) go to 140
!
     DO 130 J = 1, M
        S = ZR(I,J)
        ZR(I,J) = ZR(K,J)
        ZR(K,J) = S
        S = ZI(I,J)
        ZI(I,J) = ZI(K,J)
        ZI(K,J) = S
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
end
