subroutine CORTB (NM, LOW, IGH, AR, AI, ORTR, ORTI, M, ZR, ZI)
!
!! CORTB forms the eigenvectors of a complex general matrix from ...
!  eigenvectors of upper Hessenberg matrix output from CORTH.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      COMPLEX (ORTBAK-S, CORTB-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of a complex analogue of
!     the ALGOL procedure ORTBAK, NUM. MATH. 12, 349-368(1968)
!     by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     This subroutine forms the eigenvectors of a COMPLEX GENERAL
!     matrix by back transforming those of the corresponding
!     upper Hessenberg matrix determined by  CORTH.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, AR, AI, ZR, and ZI, as declared in the
!          calling program dimension statement.  NM is an INTEGER
!          variable.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  CBAL.  If  CBAL  has not been used,
!          set LOW=1 and IGH equal to the order of the matrix.
!
!        AR and AI contain information about the unitary trans-
!          formations used in the reduction by  CORTH  in their
!          strict lower triangles.  AR and AI are two-dimensional
!          REAL arrays, dimensioned AR(NM,IGH) and AI(NM,IGH).
!
!        ORTR and ORTI contain further information about the unitary
!          transformations used in the reduction by  CORTH.  Only
!          elements LOW through IGH are used.  ORTR and ORTI are
!          one-dimensional REAL arrays, dimensioned ORTR(IGH) and
!          ORTI(IGH).
!
!        M is the number of columns of Z=(ZR,ZI) to be back transformed.
!          M is an INTEGER variable.
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the eigenvectors to be back transformed in their first
!          M columns.  ZR and ZI are two-dimensional REAL arrays,
!          dimensioned ZR(NM,M) and ZI(NM,M).
!
!     On OUTPUT
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the transformed eigenvectors in their first M columns.
!
!        ORTR and ORTI have been altered.
!
!     Note that CORTB preserves vector Euclidean norms.
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
!***END PROLOGUE  CORTB
!
  INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
  REAL AR(NM,*),AI(NM,*),ORTR(*),ORTI(*)
  REAL ZR(NM,*),ZI(NM,*)
  REAL H,GI,GR
!
!***FIRST EXECUTABLE STATEMENT  CORTB
  if (M  ==  0) go to 200
  LA = IGH - 1
  KP1 = LOW + 1
  if (LA  <  KP1) go to 200
!     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  DO 140 MM = KP1, LA
     MP = LOW + IGH - MM
     if (AR(MP,MP-1)  ==  0.0E0 .AND. AI(MP,MP-1)  ==  0.0E0) &
        go to 140
!     .......... H BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
     H = AR(MP,MP-1) * ORTR(MP) + AI(MP,MP-1) * ORTI(MP)
     MP1 = MP + 1
!
     DO 100 I = MP1, IGH
        ORTR(I) = AR(I,MP-1)
        ORTI(I) = AI(I,MP-1)
  100    CONTINUE
!
     DO 130 J = 1, M
        GR = 0.0E0
        GI = 0.0E0
!
        DO 110 I = MP, IGH
           GR = GR + ORTR(I) * ZR(I,J) + ORTI(I) * ZI(I,J)
           GI = GI + ORTR(I) * ZI(I,J) - ORTI(I) * ZR(I,J)
  110       CONTINUE
!
        GR = GR / H
        GI = GI / H
!
        DO 120 I = MP, IGH
           ZR(I,J) = ZR(I,J) + GR * ORTR(I) - GI * ORTI(I)
           ZI(I,J) = ZI(I,J) + GR * ORTI(I) + GI * ORTR(I)
  120       CONTINUE
!
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
end
