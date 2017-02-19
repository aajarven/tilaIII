subroutine COMHES (NM, N, LOW, IGH, AR, AI, INT)
!
!! COMHES reduces a complex general matrix to complex upper Hessenberg ...
!            form using stabilized elementary similarity transformations.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B2
!***TYPE      COMPLEX (ELMHES-S, COMHES-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure COMHES,
!     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     Given a COMPLEX GENERAL matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     LOW through IGH to upper Hessenberg form by
!     stabilized elementary similarity transformations.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, AR and AI, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A=(AR,AI).  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  CBAL.  If  CBAL  has not been used,
!          set LOW=1 and IGH equal to the order of the matrix, N.
!
!        AR and AI contain the real and imaginary parts, respectively,
!          of the complex input matrix.  AR and AI are two-dimensional
!          REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
!
!     On OUTPUT
!
!        AR and AI contain the real and imaginary parts, respectively,
!          of the upper Hessenberg matrix.  The multipliers which
!          were used in the reduction are stored in the remaining
!          triangles under the Hessenberg matrix.
!
!        INT contains information on the rows and columns
!          interchanged in the reduction.  Only elements LOW through
!          IGH are used.  INT is a one-dimensional INTEGER array,
!          dimensioned INT(IGH).
!
!     Calls CDIV for complex division.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  CDIV
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COMHES
!
  INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
  REAL AR(NM,*),AI(NM,*)
  REAL XR,XI,YR,YI
  INTEGER INT(*)
!
!***FIRST EXECUTABLE STATEMENT  COMHES
  LA = IGH - 1
  KP1 = LOW + 1
  if (LA  <  KP1) go to 200
!
  DO 180 M = KP1, LA
     MM1 = M - 1
     XR = 0.0E0
     XI = 0.0E0
     I = M
!
     DO 100 J = M, IGH
        if (ABS(AR(J,MM1)) + ABS(AI(J,MM1)) &
            <=  ABS(XR) + ABS(XI)) go to 100
        XR = AR(J,MM1)
        XI = AI(J,MM1)
        I = J
  100    CONTINUE
!
     INT(M) = I
     if (I  ==  M) go to 130
!     .......... INTERCHANGE ROWS AND COLUMNS OF AR AND AI ..........
     DO 110 J = MM1, N
        YR = AR(I,J)
        AR(I,J) = AR(M,J)
        AR(M,J) = YR
        YI = AI(I,J)
        AI(I,J) = AI(M,J)
        AI(M,J) = YI
  110    CONTINUE
!
     DO 120 J = 1, IGH
        YR = AR(J,I)
        AR(J,I) = AR(J,M)
        AR(J,M) = YR
        YI = AI(J,I)
        AI(J,I) = AI(J,M)
        AI(J,M) = YI
  120    CONTINUE
!     .......... END INTERCHANGE ..........
  130    if (XR  ==  0.0E0 .AND. XI  ==  0.0E0) go to 180
     MP1 = M + 1
!
     DO 160 I = MP1, IGH
        YR = AR(I,MM1)
        YI = AI(I,MM1)
        if (YR  ==  0.0E0 .AND. YI  ==  0.0E0) go to 160
        call CDIV(YR,YI,XR,XI,YR,YI)
        AR(I,MM1) = YR
        AI(I,MM1) = YI
!
        DO 140 J = M, N
           AR(I,J) = AR(I,J) - YR * AR(M,J) + YI * AI(M,J)
           AI(I,J) = AI(I,J) - YR * AI(M,J) - YI * AR(M,J)
  140       CONTINUE

        DO J = 1, IGH
           AR(J,M) = AR(J,M) + YR * AR(J,I) - YI * AI(J,I)
           AI(J,M) = AI(J,M) + YR * AI(J,I) + YI * AR(J,I)
        end do
!
  160    CONTINUE
!
  180 CONTINUE
!
  200 RETURN
end
