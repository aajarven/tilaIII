subroutine XSETUA (IUNITA, N)
!
!! XSETUA sets logical unit numbers to which error messages are to be sent.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3B
!***TYPE      ALL (XSETUA-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XSETUA may be called to declare a list of up to five
!        logical units, each of which is to receive a copy of
!        each error message processed by this package.
!        The purpose of XSETUA is to allow simultaneous printing
!        of each error message on, say, a main output file,
!        an interactive terminal, and other files such as graphics
!        communication files.
!
!     Description of Parameters
!      --Input--
!        IUNIT - an array of up to five unit numbers.
!                Normally these numbers should all be different
!                (but duplicates are not prohibited.)
!        N     - the number of unit numbers provided in IUNIT
!                must have 1  <=  N  <=  5.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900510  Change call to XERRWV to XERMSG.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XSETUA
  DIMENSION IUNITA(5)
  CHARACTER *8 XERN1
!***FIRST EXECUTABLE STATEMENT  XSETUA
!
  if (N < 1 .OR. N > 5) THEN
     WRITE (XERN1, '(I8)') N
     call XERMSG ('SLATEC', 'XSETUA', &
        'INVALID NUMBER OF UNITS, N = ' // XERN1, 1, 2)
     return
  end if
!
  DO 10 I=1,N
     INDEX = I+4
     if (I == 1) INDEX = 3
     JUNK = J4SAVE(INDEX,IUNITA(I),.TRUE.)
   10 CONTINUE
  JUNK = J4SAVE(5,N,.TRUE.)
  return
end
