*
*   Routine for multiplying two rectangular matrices A and B
*
      subroutine matmult(A,NROWSA,B,NCOLSB,C,IDIM)
      implicit real*8 (A-H,O-Z)
      dimension C(NROWSA,NCOLSB),A(NROWSA,IDIM),
     &          B(IDIM,NCOLSB)
*
*   put zeroes in matrix C
*
      do 5 M=1,NCOLSB,1
         do 7 I=1,NROWSA,1
            C(I,M) = 0.
7        continue
5     continue
      do 10 L=1,NCOLSB,1
         do 20 K=1,NROWSA,1
            do 30 J=1,IDIM,1
               C(K,L) = C(K,L) + A(K,J) * B(J,L)
30          continue
20       continue
10    continue
      return
      end
