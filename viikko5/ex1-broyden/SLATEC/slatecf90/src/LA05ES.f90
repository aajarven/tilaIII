subroutine LA05ES (A, IRN, IP, N, IW, IA, REALS)
!
!! LA05ES is subsidiary to SPLP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LA05ES-S, LA05ED-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
!     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
!     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
!     THE FINAL LETTER =S= IN THE NAMES USED HERE.
!     REVISED SEP. 13, 1979.
!
!     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
!     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
!     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
!     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
!     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
!
!***SEE ALSO  SPLP
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    LA05DS
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  LA05ES
  LOGICAL REALS
  REAL A(*)
  INTEGER IRN(*), IW(*)
  INTEGER IP(*)
  COMMON /LA05DS/ SMALL, LP, LENL, LENU, NCP, LROW, LCOL
!***FIRST EXECUTABLE STATEMENT  LA05ES
  NCP = NCP + 1
!     COMPRESS FILE OF POSITIVE INTEGERS. ENTRY J STARTS AT IRN(IP(J))
!  AND CONTAINS IW(J) INTEGERS,J=1,N. OTHER COMPONENTS OF IRN ARE ZERO.
!  LENGTH OF COMPRESSED FILE PLACED IN LROW if REALS IS .TRUE. OR LCOL
!  OTHERWISE.
!  if REALS IS .TRUE. ARRAY A CONTAINS A REAL FILE ASSOCIATED WITH IRN
!  AND THIS IS COMPRESSED TOO.
!  A,IRN,IP,IW,IA ARE INPUT/OUTPUT VARIABLES.
!  N,REALS ARE INPUT/UNCHANGED VARIABLES.
!
  DO 10 J=1,N
! STORE THE LAST ELEMENT OF ENTRY J IN IW(J) THEN OVERWRITE IT BY -J.
     NZ = IW(J)
     if (NZ <= 0) go to 10
     K = IP(J) + NZ - 1
     IW(J) = IRN(K)
     IRN(K) = -J
   10 CONTINUE
! KN IS THE POSITION OF NEXT ENTRY IN COMPRESSED FILE.
  KN = 0
  IPI = 0
  KL = LCOL
  if (REALS) KL = LROW
! LOOP THROUGH THE OLD FILE SKIPPING ZERO (DUMMY) ELEMENTS AND
!     MOVING GENUINE ELEMENTS FORWARD. THE ENTRY NUMBER BECOMES
!     KNOWN ONLY WHEN ITS END IS DETECTED BY THE PRESENCE OF A NEGATIVE
!     INTEGER.
  DO 30 K=1,KL
     if (IRN(K) == 0) go to 30
     KN = KN + 1
     if (REALS) A(KN) = A(K)
     if (IRN(K) >= 0) go to 20
! END OF ENTRY. RESTORE IRN(K), SET POINTER TO START OF ENTRY AND
!     STORE CURRENT KN IN IPI READY FOR USE WHEN NEXT LAST ENTRY
!     IS DETECTED.
     J = -IRN(K)
     IRN(K) = IW(J)
     IP(J) = IPI + 1
     IW(J) = KN - IPI
     IPI = KN
   20    IRN(KN) = IRN(K)
   30 CONTINUE
  if (REALS) LROW = KN
  if (.NOT.REALS) LCOL = KN
  return
end
