      SUBROUTINE LMEMCHK
*
* Check memory allocated  with the memory manager
*
      CALL MEMMAN(IDUM,IDUM,'CHECK ',IDUM,'IDUM  ')
      WRITE(6,*) ' Memory check passed '
*
      RETURN
      END
