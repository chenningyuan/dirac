      subroutine abend1(CHARINP)
      implicit none
      character*(*) CHARINP
#ifdef LUCIAREL_STANDALONE
      write(6,*) CHARINP
      stop 'Abend called'
#else
      call quit(charinp)
#endif
      return
      end
