* Copy vector to defined partition of matrix
* (Simply the way the vector is arranged in memory)
*
      subroutine vectomat(VEC,WMAT,ILEN,NROW,NCOL,
     &                    NPARTROW,NPARTCOL,IROWOFF,ICOLOFF)
*
      implicit real*8 (A-H,O-Z)
*
      dimension WMAT(NROW,NCOL)
      dimension VEC(ILEN)
*
*   VEC: Input matrix in vector format
*   WMAT: Output matrix (can be a partition of a bigger matrix!)
*   ILEN: Length of input vector
*   NROW: Rows of WMAT
*   NCOL: Columns of WMAT
*   NPARTROW: Number of rows in output partition
*   NPARTCOL: Number of columns in output partition
*   IROWOFF: Output row offset
*   ICOLOFF: Output column offset
*
      K = 1
      do I=ICOLOFF,(ICOLOFF+NPARTCOL-1),1
         do J=IROWOFF,(IROWOFF+NPARTROW-1),1
            WMAT(J,I) = VEC(K)
            K = K + 1
         end do
      end do
      return
      end
