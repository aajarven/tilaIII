subroutine CSROT (N, CX, INCX, CY, INCY, C, S)
!
!! CSROT applies a plane Givens rotation.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B10
!***TYPE      COMPLEX (SROT-S, DROT-D, CSROT-C)
!***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
!             LINEAR ALGEBRA, PLANE ROTATION, VECTOR
!***AUTHOR  Dongarra, J., (ANL)
!***DESCRIPTION
!
!     CSROT applies the complex Givens rotation
!
!          (X)   ( C S)(X)
!          (Y) = (-S C)(Y)
!
!     N times where for I = 0,...,N-1
!
!          X = CX(LX+I*INCX)
!          Y = CY(LY+I*INCY),
!
!     where LX = 1 if INCX  >=  0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!     Argument Description
!
!        N      (integer)  number of elements in each vector
!
!        CX     (complex array)  beginning of one vector
!
!        INCX   (integer)  memory spacing of successive elements
!               of vector CX
!
!        CY     (complex array)  beginning of the other vector
!
!        INCY   (integer)  memory spacing of successive elements
!               of vector CY
!
!        C      (real)  cosine term of the rotation
!
!        S      (real)  sine term of the rotation.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   810223  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CSROT
  COMPLEX CX(*), CY(*), CTEMP
  REAL C, S
  INTEGER I, INCX, INCY, IX, IY, N
!***FIRST EXECUTABLE STATEMENT  CSROT
  if (N  <=  0) RETURN
  if (INCX == 1 .AND. INCY == 1)go to 20
!
!     Code for unequal increments or equal increments not equal to 1.
!
  IX = 1
  IY = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  if (INCY  <  0) IY = (-N+1)*INCY + 1
  DO 10 I = 1,N
    CTEMP = C*CX(IX) + S*CY(IY)
    CY(IY) = C*CY(IY) - S*CX(IX)
    CX(IX) = CTEMP
    IX = IX + INCX
    IY = IY + INCY
   10 CONTINUE
  return
!
!     Code for both increments equal to 1.
!
   20 DO 30 I = 1,N
    CTEMP = C*CX(I) + S*CY(I)
    CY(I) = C*CY(I) - S*CX(I)
    CX(I) = CTEMP
   30 CONTINUE
  return
end
