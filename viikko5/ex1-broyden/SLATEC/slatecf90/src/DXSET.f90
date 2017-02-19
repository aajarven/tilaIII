subroutine DXSET (IRAD, NRADPL, DZERO, NBITS, IERROR)
!
!! DXSET provides double-precision floating-point arithmetic ...
!            with an extended exponent range.
!***LIBRARY   SLATEC
!***CATEGORY  A3D
!***TYPE      DOUBLE PRECISION (XSET-S, DXSET-D)
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
!           Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!
!   SUBROUTINE  DXSET  MUST BE CALLED PRIOR TO CALLING ANY OTHER
! EXTENDED-RANGE SUBROUTINE. IT CALCULATES AND STORES SEVERAL
! MACHINE-DEPENDENT CONSTANTS IN COMMON BLOCKS. THE USER MUST
! SUPPLY FOUR CONSTANTS THAT PERTAIN TO HIS PARTICULAR COMPUTER.
! THE CONSTANTS ARE
!
!          IRAD = THE INTERNAL BASE OF DOUBLE-PRECISION
!                 ARITHMETIC IN THE COMPUTER.
!        NRADPL = THE NUMBER OF RADIX PLACES CARRIED IN
!                 THE DOUBLE-PRECISION REPRESENTATION.
!         DZERO = THE SMALLEST OF 1/DMIN, DMAX, DMAXLN WHERE
!                 DMIN = THE SMALLEST POSITIVE DOUBLE-PRECISION
!                 NUMBER OR AN UPPER BOUND TO THIS NUMBER,
!                 DMAX = THE LARGEST DOUBLE-PRECISION NUMBER
!                 OR A LOWER BOUND TO THIS NUMBER,
!                 DMAXLN = THE LARGEST DOUBLE-PRECISION NUMBER
!                 SUCH THAT LOG10(DMAXLN) CAN BE COMPUTED BY THE
!                 FORTRAN SYSTEM (ON MOST SYSTEMS DMAXLN = DMAX).
!         NBITS = THE NUMBER OF BITS (EXCLUSIVE OF SIGN) IN
!                 AN INTEGER COMPUTER WORD.
!
! ALTERNATIVELY, ANY OR ALL OF THE CONSTANTS CAN BE GIVEN
! THE VALUE 0 (0.0D0 FOR DZERO). if A CONSTANT IS ZERO, DXSET TRIES
! TO ASSIGN AN APPROPRIATE VALUE BY CALLING I1MACH
! (SEE P.A.FOX, A.D.HALL, N.L.SCHRYER, ALGORITHM 528 FRAMEWORK
! FOR A PORTABLE LIBRARY, ACM TRANSACTIONS ON MATH SOFTWARE,
! V.4, NO.2, JUNE 1978, 177-188).
!
!   THIS IS THE SETTING-UP SUBROUTINE FOR A PACKAGE OF SUBROUTINES
! THAT FACILITATE THE USE OF EXTENDED-RANGE ARITHMETIC. EXTENDED-RANGE
! ARITHMETIC ON A PARTICULAR COMPUTER IS DEFINED ON THE SET OF NUMBERS
! OF THE FORM
!
!               (X,IX) = X*RADIX**IX
!
! WHERE X IS A DOUBLE-PRECISION NUMBER CALLED THE PRINCIPAL PART,
! IX IS AN INTEGER CALLED THE AUXILIARY INDEX, AND RADIX IS THE
! INTERNAL BASE OF THE DOUBLE-PRECISION ARITHMETIC.  OBVIOUSLY,
! EACH REAL NUMBER IS REPRESENTABLE WITHOUT ERROR BY MORE THAN ONE
! EXTENDED-RANGE FORM.  CONVERSIONS BETWEEN  DIFFERENT FORMS ARE
! ESSENTIAL IN CARRYING OUT ARITHMETIC OPERATIONS.  WITH THE CHOICE
! OF RADIX WE HAVE MADE, AND THE SUBROUTINES WE HAVE WRITTEN, THESE
! CONVERSIONS ARE PERFORMED WITHOUT ERROR (AT LEAST ON MOST COMPUTERS).
! (SEE SMITH, J.M., OLVER, F.W.J., AND LOZIER, D.W., EXTENDED-RANGE
! ARITHMETIC AND NORMALIZED LEGENDRE POLYNOMIALS, ACM TRANSACTIONS ON
! MATHEMATICAL SOFTWARE, MARCH 1981).
!
!   AN EXTENDED-RANGE NUMBER  (X,IX)  IS SAID TO BE IN ADJUSTED FORM IF
! X AND IX ARE ZERO OR
!
!           RADIX**(-L)  <=  ABS(X)  <  RADIX**L
!
! IS SATISFIED, WHERE L IS A COMPUTER-DEPENDENT INTEGER DEFINED IN THIS
! SUBROUTINE. TWO EXTENDED-RANGE NUMBERS IN ADJUSTED FORM CAN BE ADDED,
! SUBTRACTED, MULTIPLIED OR DIVIDED (IF THE DIVISOR IS NONZERO) WITHOUT
! CAUSING OVERFLOW OR UNDERFLOW IN THE PRINCIPAL PART OF THE RESULT.
! WITH PROPER USE OF THE EXTENDED-RANGE SUBROUTINES, THE ONLY OVERFLOW
! THAT CAN OCCUR IS INTEGER OVERFLOW IN THE AUXILIARY INDEX. if THIS
! IS DETECTED, THE SOFTWARE CALLS XERROR (A GENERAL ERROR-HANDLING
! FORTRAN SUBROUTINE PACKAGE).
!
!   MULTIPLICATION AND DIVISION IS PERFORMED BY SETTING
!
!                 (X,IX)*(Y,IY) = (X*Y,IX+IY)
! OR
!                 (X,IX)/(Y,IY) = (X/Y,IX-IY).
!
! PRE-ADJUSTMENT OF THE OPERANDS IS ESSENTIAL TO AVOID
! OVERFLOW OR  UNDERFLOW OF THE PRINCIPAL PART. SUBROUTINE
! DXADJ (SEE BELOW) MAY BE CALLED TO TRANSFORM ANY EXTENDED-
! RANGE NUMBER INTO ADJUSTED FORM.
!
!   ADDITION AND SUBTRACTION REQUIRE THE USE OF SUBROUTINE DXADD
! (SEE BELOW).  THE INPUT OPERANDS NEED NOT BE IN ADJUSTED FORM.
! HOWEVER, THE RESULT OF ADDITION OR SUBTRACTION IS RETURNED
! IN ADJUSTED FORM.  THUS, FOR EXAMPLE, if (X,IX),(Y,IY),
! (U,IU),  AND (V,IV) ARE IN ADJUSTED FORM, THEN
!
!                 (X,IX)*(Y,IY) + (U,IU)*(V,IV)
!
! CAN BE COMPUTED AND STORED IN ADJUSTED FORM WITH NO EXPLICIT
! CALLS TO DXADJ.
!
!   WHEN AN EXTENDED-RANGE NUMBER IS TO BE PRINTED, IT MUST BE
! CONVERTED TO AN EXTENDED-RANGE FORM WITH DECIMAL RADIX.  SUBROUTINE
! DXCON IS PROVIDED FOR THIS PURPOSE.
!
!   THE SUBROUTINES CONTAINED IN THIS PACKAGE ARE
!
!     SUBROUTINE DXADD
! USAGE
!                  call DXADD(X,IX,Y,IY,Z,IZ,IERROR)
!                  if (IERROR /= 0) RETURN
! DESCRIPTION
!                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) =
!                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED
!                  BEFORE RETURNING. THE INPUT OPERANDS
!                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR
!                  PRINCIPAL PARTS MUST SATISFY
!                  RADIX**(-2L) <= ABS(X) <= RADIX**(2L),
!                  RADIX**(-2L) <= ABS(Y) <= RADIX**(2L).
!
!     SUBROUTINE DXADJ
! USAGE
!                  call DXADJ(X,IX,IERROR)
!                  if (IERROR /= 0) RETURN
! DESCRIPTION
!                  TRANSFORMS (X,IX) SO THAT
!                  RADIX**(-L)  <=  ABS(X)  <  RADIX**L.
!                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
!                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
!                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC.
!
!     SUBROUTINE DXC210
! USAGE
!                  call DXC210(K,Z,J,IERROR)
!                  if (IERROR /= 0) RETURN
! DESCRIPTION
!                  GIVEN K THIS SUBROUTINE COMPUTES J AND Z
!                  SUCH THAT  RADIX**K = Z*10**J, WHERE Z IS IN
!                  THE RANGE 1/10  <=  Z  <  1.
!                  THE VALUE OF Z WILL BE ACCURATE TO FULL
!                  DOUBLE-PRECISION PROVIDED THE NUMBER
!                  OF DECIMAL PLACES IN THE LARGEST
!                  INTEGER PLUS THE NUMBER OF DECIMAL
!                  PLACES CARRIED IN DOUBLE-PRECISION DOES NOT
!                  EXCEED 60. DXC210 IS CALLED BY SUBROUTINE
!                  DXCON WHEN NECESSARY. THE USER SHOULD
!                  NEVER NEED TO call DXC210 DIRECTLY.
!
!     SUBROUTINE DXCON
! USAGE
!                  call DXCON(X,IX,IERROR)
!                  if (IERROR /= 0) RETURN
! DESCRIPTION
!                  CONVERTS (X,IX) = X*RADIX**IX
!                  TO DECIMAL FORM IN PREPARATION FOR
!                  PRINTING, SO THAT (X,IX) = X*10**IX
!                  WHERE 1/10  <=  ABS(X)  <  1
!                  IS RETURNED, EXCEPT THAT IF
!                  (ABS(X),IX) IS BETWEEN RADIX**(-2L)
!                  AND RADIX**(2L) THEN THE REDUCED
!                  FORM WITH IX = 0 IS RETURNED.
!
!     SUBROUTINE DXRED
! USAGE
!                  call DXRED(X,IX,IERROR)
!                  if (IERROR /= 0) RETURN
! DESCRIPTION
!                  IF
!                  RADIX**(-2L)  <=  (ABS(X),IX)  <=  RADIX**(2L)
!                  THEN DXRED TRANSFORMS (X,IX) SO THAT IX=0.
!                  if (X,IX) IS OUTSIDE THE ABOVE RANGE,
!                  THEN DXRED TAKES NO ACTION.
!                  THIS SUBROUTINE IS USEFUL if THE
!                  RESULTS OF EXTENDED-RANGE CALCULATIONS
!                  ARE TO BE USED IN SUBSEQUENT ORDINARY
!                  DOUBLE-PRECISION CALCULATIONS.
!
!***REFERENCES  Smith, Olver and Lozier, Extended-Range Arithmetic and
!                 Normalized Legendre Polynomials, ACM Trans on Math
!                 Softw, v 7, n 1, March 1981, pp 93--105.
!***ROUTINES CALLED  I1MACH, XERMSG
!***COMMON BLOCKS    DXBLK1, DXBLK2, DXBLK3
!***REVISION HISTORY  (YYMMDD)
!   820712  DATE WRITTEN
!   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXSET
  INTEGER IRAD, NRADPL, NBITS
  DOUBLE PRECISION DZERO, DZEROX
  COMMON /DXBLK1/ NBITSF
  SAVE /DXBLK1/
  DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
  INTEGER L, L2, KMAX
  COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
  SAVE /DXBLK2/
  INTEGER NLG102, MLG102, LG102
  COMMON /DXBLK3/ NLG102, MLG102, LG102(21)
  SAVE /DXBLK3/
  INTEGER IFLAG
  SAVE IFLAG
!
  DIMENSION LOG102(20), LGTEMP(20)
  SAVE LOG102
!
!   LOG102 CONTAINS THE FIRST 60 DIGITS OF LOG10(2) FOR USE IN
! CONVERSION OF EXTENDED-RANGE NUMBERS TO BASE 10 .
  DATA LOG102 /301,029,995,663,981,195,213,738,894,724,493,026,768, &
   189,881,462,108,541,310,428/
!
! FOLLOWING CODING PREVENTS DXSET FROM BEING EXECUTED MORE THAN ONCE.
! THIS IS IMPORTANT BECAUSE SOME SUBROUTINES (SUCH AS DXNRMP AND
! DXLEGF) call DXSET TO MAKE SURE EXTENDED-RANGE ARITHMETIC HAS
! BEEN INITIALIZED. THE USER MAY WANT TO PRE-EMPT THIS CALL, FOR
! EXAMPLE WHEN I1MACH IS NOT AVAILABLE. SEE CODING BELOW.
  DATA IFLAG /0/
!***FIRST EXECUTABLE STATEMENT  DXSET
  IERROR=0
  if (IFLAG  /=  0) RETURN
  IRADX = IRAD
  NRDPLC = NRADPL
  DZEROX = DZERO
  IMINEX = 0
  IMAXEX = 0
  NBITSX = NBITS
! FOLLOWING 5 STATEMENTS SHOULD BE DELETED if I1MACH IS
! NOT AVAILABLE OR NOT CONFIGURED TO RETURN THE CORRECT
! MACHINE-DEPENDENT VALUES.
  if (IRADX  ==  0) IRADX = I1MACH (10)
  if (NRDPLC  ==  0) NRDPLC = I1MACH (14)
  if (DZEROX  ==  0.0D0) IMINEX = I1MACH (15)
  if (DZEROX  ==  0.0D0) IMAXEX = I1MACH (16)
  if (NBITSX  ==  0) NBITSX = I1MACH (8)
  if (IRADX == 2) go to 10
  if (IRADX == 4) go to 10
  if (IRADX == 8) go to 10
  if (IRADX == 16) go to 10
  call XERMSG ('SLATEC', 'DXSET', 'IMPROPER VALUE OF IRAD', 201, 1)
  IERROR=201
  return
   10 CONTINUE
  LOG2R=0
  if (IRADX == 2) LOG2R = 1
  if (IRADX == 4) LOG2R = 2
  if (IRADX == 8) LOG2R = 3
  if (IRADX == 16) LOG2R = 4
  NBITSF=LOG2R*NRDPLC
  RADIX = IRADX
  DLG10R = LOG10(RADIX)
  if (DZEROX  /=  0.0D0) go to 14
  LX = MIN ((1-IMINEX)/2, (IMAXEX-1)/2)
  go to 16
   14 LX = 0.5D0*LOG10(DZEROX)/DLG10R
! RADIX**(2*L) SHOULD NOT OVERFLOW, BUT REDUCE L BY 1 FOR FURTHER
! PROTECTION.
  LX=LX-1
   16 L2 = 2*LX
  if (LX >= 4) go to 20
  call XERMSG ('SLATEC', 'DXSET', 'IMPROPER VALUE OF DZERO', 202, 1)
  IERROR=202
  return
   20 L = LX
  RADIXL = RADIX**L
  RAD2L = RADIXL**2
!    IT IS NECESSARY TO RESTRICT NBITS (OR NBITSX) TO BE LESS THAN SOME
! UPPER LIMIT BECAUSE OF BINARY-TO-DECIMAL CONVERSION. SUCH CONVERSION
! IS DONE BY DXC210 AND REQUIRES A CONSTANT THAT IS STORED TO SOME FIXED
! PRECISION. THE STORED CONSTANT (LOG102 IN THIS ROUTINE) PROVIDES
! FOR CONVERSIONS ACCURATE TO THE LAST DECIMAL DIGIT WHEN THE INTEGER
! WORD LENGTH DOES NOT EXCEED 63. A LOWER LIMIT OF 15 BITS IS IMPOSED
! BECAUSE THE SOFTWARE IS DESIGNED TO RUN ON COMPUTERS WITH INTEGER WORD
! LENGTH OF AT LEAST 16 BITS.
  if (15 <= NBITSX .AND. NBITSX <= 63) go to 30
  call XERMSG ('SLATEC', 'DXSET', 'IMPROPER VALUE OF NBITS', 203, 1)
  IERROR=203
  return
   30 CONTINUE
  KMAX = 2**(NBITSX-1) - L2
  NB = (NBITSX-1)/2
  MLG102 = 2**NB
  if (1 <= NRDPLC*LOG2R .AND. NRDPLC*LOG2R <= 120) go to 40
  call XERMSG ('SLATEC', 'DXSET', 'IMPROPER VALUE OF NRADPL', 204, &
               1)
  IERROR=204
  return
   40 CONTINUE
  NLG102 = NRDPLC*LOG2R/NB + 3
  NP1 = NLG102 + 1
!
!   AFTER COMPLETION OF THE FOLLOWING LOOP, IC CONTAINS
! THE INTEGER PART AND LGTEMP CONTAINS THE FRACTIONAL PART
! OF LOG10(IRADX) IN RADIX 1000.
  IC = 0
  DO 50 II=1,20
    I = 21 - II
    IT = LOG2R*LOG102(I) + IC
    IC = IT/1000
    LGTEMP(I) = MOD(IT,1000)
   50 CONTINUE
!
!   AFTER COMPLETION OF THE FOLLOWING LOOP, LG102 CONTAINS
! LOG10(IRADX) IN RADIX MLG102. THE RADIX POINT IS
! BETWEEN LG102(1) AND LG102(2).
  LG102(1) = IC
  DO 80 I=2,NP1
    LG102X = 0
    DO 70 J=1,NB
      IC = 0
      DO 60 KK=1,20
        K = 21 - KK
        IT = 2*LGTEMP(K) + IC
        IC = IT/1000
        LGTEMP(K) = MOD(IT,1000)
   60     CONTINUE
      LG102X = 2*LG102X + IC
   70   CONTINUE
    LG102(I) = LG102X
   80 CONTINUE
!
! CHECK SPECIAL CONDITIONS REQUIRED BY SUBROUTINES...
  if (NRDPLC < L) go to 90
  call XERMSG ('SLATEC', 'DXSET', 'NRADPL  >=  L', 205, 1)
  IERROR=205
  return
   90 if (6*L <= KMAX) go to 100
  call XERMSG ('SLATEC', 'DXSET', '6*L  >  KMAX', 206, 1)
  IERROR=206
  return
  100 CONTINUE
  IFLAG = 1
  return
end
