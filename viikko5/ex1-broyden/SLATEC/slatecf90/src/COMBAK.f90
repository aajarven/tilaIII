subroutine COMBAK (NM, LOW, IGH, AR, AI, INT, M, ZR, ZI)
!
!! COMBAK forms the eigenvectors of a complex general matrix from the ...
!  eigenvectors of a upper Hessenberg matrix output from COMHES.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      COMPLEX (ELMBAK-S, COMBAK-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure COMBAK,
!     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     This subroutine forms the eigenvectors of a COMPLEX GENERAL
!     matrix by back transforming those of the corresponding
!     upper Hessenberg matrix determined by  COMHES.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, AR, AI, ZR and ZI, as declared in the
!          calling program dimension statement.  NM is an INTEGER
!          variable.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  CBAL.  If  CBAL  has not been used,
!          set LOW=1 and IGH equal to the order of the matrix.
!
!        AR and AI contain the multipliers which were used in the
!           reduction by  COMHES  in their lower triangles below
!           the subdiagonal.  AR and AI are two-dimensional REAL
!           arrays, dimensioned AR(NM,IGH) and AI(NM,IGH).
!
!        INT contains information on the rows and columns
!          interchanged in the reduction by  COMHES.  Only
!          elements LOW through IGH are used.  INT is a
!          one-dimensional INTEGER array, dimensioned INT(IGH).
!
!        M is the number of eigenvectors to be back transformed.
!          M is an INTEGER variable.
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the eigenvectors to be back transformed in their first M
!          columns.  ZR and ZI are two-dimensional REAL arrays,
!          dimensioned ZR(NM,M) and ZI(NM,M).
!
!     On OUTPUT
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the transformed eigenvectors in their first M columns.
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
!***END PROLOGUE  COMBAK
!
  INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
  REAL AR(NM,*),AI(NM,*),ZR(NM,*),ZI(NM,*)
  REAL XR,XI
  INTEGER INT(*)
!
!***FIRST EXECUTABLE STATEMENT  COMBAK
  if (M  ==  0) go to 200
  LA = IGH - 1
  KP1 = LOW + 1
  if (LA  <  KP1) go to 200
!     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  DO 140 MM = KP1, LA
     MP = LOW + IGH - MM
     MP1 = MP + 1
!
     DO 110 I = MP1, IGH
        XR = AR(I,MP-1)
        XI = AI(I,MP-1)
        if (XR  ==  0.0E0 .AND. XI  ==  0.0E0) go to 110
!
        DO 100 J = 1, M
           ZR(I,J) = ZR(I,J) + XR * ZR(MP,J) - XI * ZI(MP,J)
           ZI(I,J) = ZI(I,J) + XR * ZI(MP,J) + XI * ZR(MP,J)
  100       CONTINUE
!
  110    CONTINUE
!
     I = INT(MP)
     if (I  ==  MP) go to 140
!
     DO 130 J = 1, M
        XR = ZR(I,J)
        ZR(I,J) = ZR(MP,J)
        ZR(MP,J) = XR
        XI = ZI(I,J)
        ZI(I,J) = ZI(MP,J)
        ZI(MP,J) = XI
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
end
