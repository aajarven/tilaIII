subroutine XSETUN (IUNIT)
!
!! XSETUN sets output file to which XERROR messages are to be sent.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3B
!***TYPE      ALL (XSETUN-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XSETUN sets the output file to which error messages are to
!        be sent.  Only one file will be used.  See XSETUA for
!        how to declare more than one file.
!
!     Description of Parameter
!      --Input--
!        IUNIT - an input parameter giving the logical unit number
!                to which error messages are to be sent.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XSETUN
!***FIRST EXECUTABLE STATEMENT  XSETUN
  JUNK = J4SAVE(3,IUNIT,.TRUE.)
  JUNK = J4SAVE(5,1,.TRUE.)
  return
end
