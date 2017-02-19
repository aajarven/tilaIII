subroutine DCDOT (N, FM, CX, INCX, CY, INCY, DCR, DCI)
!
!! DCDOT computes the inner product of two vectors with extended ...
!  precision accumulation and result.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      COMPLEX (DSDOT-D, DCDOT-C)
!***KEYWORDS  BLAS, COMPLEX VECTORS, DOT PRODUCT, INNER PRODUCT,
!             LINEAR ALGEBRA, VECTOR
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!    Compute the dot product of 2 complex vectors, CX and CY, e.g.
!    CX DOT CY, or, CXconjugate DOT CY.  The real and imaginary
!    parts of CX and CY are converted to double precision, the dot
!    product accumulation is done in double precision and the output
!    is given as 2 double precision numbers, corresponding to the real
!    and imaginary part of the result.
!     Input
!      N:  Number of complex components of CX and CY.
!      FM: =+1.0   compute CX DOT CY.
!          =-1.0   compute CXconjugate DOT CY.
!      CX(N):
!      CY(N):  Complex arrays of length N.
!      INCX:(Integer)   Spacing of elements of CX to use
!      INCY:(Integer)   Spacing of elements of CY to use.
!     Output
!      DCR:(Double Precision) Real part of dot product.
!      DCI:(Double Precision) Imaginary part of dot product.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DCDOT
  INTEGER I, INCX, INCY, KX, KY, N
  COMPLEX CX(*), CY(*)
  DOUBLE PRECISION DCR, DCI, DT1, DT2, DT3, DT4, FM
!***FIRST EXECUTABLE STATEMENT  DCDOT
  DCR = 0.0D0
  DCI = 0.0D0
  if (N  <=  0) go to 20
!
  KX = 1
  KY = 1
  if (INCX  <  0) KX = 1+(1-N)*INCX
  if (INCY  <  0) KY = 1+(1-N)*INCY
  DO 10 I = 1,N
    DT1 = DBLE(REAL(CX(KX)))
    DT2 = DBLE(REAL(CY(KY)))
    DT3 = DBLE(AIMAG(CX(KX)))
    DT4 = DBLE(AIMAG(CY(KY)))
    DCR = DCR+(DT1*DT2)-FM*(DT3*DT4)
    DCI = DCI+(DT1*DT4)+FM*(DT3*DT2)
    KX = KX+INCX
    KY = KY+INCY
   10 CONTINUE
   20 RETURN
end
