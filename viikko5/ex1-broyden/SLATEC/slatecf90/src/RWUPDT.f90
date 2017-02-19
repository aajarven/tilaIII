subroutine RWUPDT (N, R, LDR, W, B, ALPHA, COS, SIN)
!
!! RWUPDT is subsidiary to SNLS1 and SNLS1E.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (RWUPDT-S, DWUPDT-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an N by N upper triangular matrix R, this subroutine
!     computes the QR decomposition of the matrix formed when a row
!     is added to R. If the row is specified by the vector W, then
!     RWUPDT determines an orthogonal matrix Q such that when the
!     N+1 by N matrix composed of R augmented by W is premultiplied
!     by (Q TRANSPOSE), the resulting matrix is upper trapezoidal.
!     The orthogonal matrix Q is the product of N transformations
!
!           G(1)*G(2)* ... *G(N)
!
!     where G(I) is a Givens rotation in the (I,N+1) plane which
!     eliminates elements in the I-th plane. RWUPDT also
!     computes the product (Q TRANSPOSE)*C where C is the
!     (N+1)-vector (b,alpha). Q itself is not accumulated, rather
!     the information to recover the G rotations is supplied.
!
!     The subroutine statement is
!
!       SUBROUTINE RWUPDT(N,R,LDR,W,B,ALPHA,COS,SIN)
!
!     where
!
!       N is a positive integer input variable set to the order of R.
!
!       R is an N by N array. On input the upper triangular part of
!         R must contain the matrix to be updated. On output R
!         contains the updated triangular matrix.
!
!       LDR is a positive integer input variable not less than N
!         which specifies the leading dimension of the array R.
!
!       W is an input array of length N which must contain the row
!         vector to be added to R.
!
!       B is an array of length N. On input B must contain the
!         first N elements of the vector C. On output B contains
!         the first N elements of the vector (Q TRANSPOSE)*C.
!
!       ALPHA is a variable. On input ALPHA must contain the
!         (N+1)-st element of the vector C. On output ALPHA contains
!         the (N+1)-st element of the vector (Q TRANSPOSE)*C.
!
!       COS is an output array of length N which contains the
!         cosines of the transforming Givens rotations.
!
!       SIN is an output array of length N which contains the
!         sines of the transforming Givens rotations.
!
!***SEE ALSO  SNLS1, SNLS1E
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  RWUPDT
  INTEGER N,LDR
  REAL ALPHA
  REAL R(LDR,*),W(*),B(*),COS(*),SIN(*)
  INTEGER I,J,JM1
  REAL COTAN,ONE,P5,P25,ROWJ,TAN,TEMP,ZERO
  SAVE ONE, P5, P25, ZERO
  DATA ONE,P5,P25,ZERO /1.0E0,5.0E-1,2.5E-1,0.0E0/
!***FIRST EXECUTABLE STATEMENT  RWUPDT
  DO 60 J = 1, N
     ROWJ = W(J)
     JM1 = J - 1
!
!        APPLY THE PREVIOUS TRANSFORMATIONS TO
!        R(I,J), I=1,2,...,J-1, AND TO W(J).
!
     if (JM1  <  1) go to 20
     DO 10 I = 1, JM1
        TEMP = COS(I)*R(I,J) + SIN(I)*ROWJ
        ROWJ = -SIN(I)*R(I,J) + COS(I)*ROWJ
        R(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
!
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES W(J).
!
     COS(J) = ONE
     SIN(J) = ZERO
     if (ROWJ  ==  ZERO) go to 50
     if (ABS(R(J,J))  >=  ABS(ROWJ)) go to 30
        COTAN = R(J,J)/ROWJ
        SIN(J) = P5/SQRT(P25+P25*COTAN**2)
        COS(J) = SIN(J)*COTAN
        go to 40
   30    CONTINUE
        TAN = ROWJ/R(J,J)
        COS(J) = P5/SQRT(P25+P25*TAN**2)
        SIN(J) = COS(J)*TAN
   40    CONTINUE
!
!        APPLY THE CURRENT TRANSFORMATION TO R(J,J), B(J), AND ALPHA.
!
     R(J,J) = COS(J)*R(J,J) + SIN(J)*ROWJ
     TEMP = COS(J)*B(J) + SIN(J)*ALPHA
     ALPHA = -SIN(J)*B(J) + COS(J)*ALPHA
     B(J) = TEMP
   50    CONTINUE
   60    CONTINUE
  return
!
!     LAST CARD OF SUBROUTINE RWUPDT.
!
end
