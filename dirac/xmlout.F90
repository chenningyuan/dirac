! Read file from file_path and put it in dirac.xml
! This is not very robust, one should use base64 encoding
subroutine xmlout_embed_formatted(file_path,file_type)
  implicit none
  character file_path *(*), file_type *(*)
  character linebuf*(256)
  integer ios, lu
#include "xmlout.h"
  write(LUXML,*) '    <file name="'//trim(file_path)//'" type="'//trim(file_type)//'" encoding="none">'
  lu = 38
  OPEN(lu,FILE = file_path,IOSTAT=IOS)
  IF (IOS.NE.0) THEN
     CALL QUIT('Error in opening file for xml embed!')
  ENDIF
  do while (.true.)
     read(lu,*,end=666) linebuf
     write(luxml,*),trim(linebuf)
  enddo
666 continue
  write(LUXML,*) '    </file>'
  close(lu)
end subroutine xmlout_embed_formatted

subroutine xmlout_begin
  implicit none
#include "xmlout.h"
!     Start the xml <calculation> and echo the input to the xml output
  print *,'ulf remove: xml output started'
  write(LUXML,*) '<calculation>'
  write(LUXML,*) '  <input>'
  write(LUXML,*) '    <program>dirac</program>'
  call xmlout_embed_formatted('DIRAC.INP','dirac-input')
  call xmlout_embed_formatted('MOLECULE.MOL','dirac-mol')
  write(luxml,*) '   </input>'
end subroutine xmlout_begin

subroutine xmlout_end
  implicit none
#include "xmlout.h"
  write(LUXML,*) '</calculation>'
  close(LUXML,STATUS='KEEP')
end subroutine xmlout_end
