C Output from Public domain Ratfor, version 1.01
C Modified to remove compiler warnings in pedantic mode. SPG, 11/3/16
      subroutine robcovf(n, p, nc, start, len, u, s, v, w)
      implicit none
C      implicit double precision (a-h,o-z)
      integer n,i,j,k,nc
      integer p, start(nc), len(nc)
      double precision u(n,p), s(p), v(p,p), w(p,p)
      do23000 i=1,p 
      do23002 j=1,p 
      w(i,j)=0d0 
23002 continue
      continue
23000 continue
      continue
      do23004 k=1,nc 
      do23006 i=1,p 
      s(i)=0d0
      do23008 j=1,p 
      v(i,j)=0d0 
23008 continue
      continue
23006 continue
      continue
      do23010 i=start(k),start(k)+len(k)-1 
      do23012 j=1,p 
      s(j)=s(j)+u(i,j) 
23012 continue
      continue
23010 continue
      continue
      do23014 i=1,p 
      do23016 j=1,p 
      v(i,j)=v(i,j)+s(i)*s(j)
23016 continue
      continue
23014 continue
      continue
      do23018 i=1,p 
      do23020 j=1,p 
      w(i,j)=w(i,j)+v(i,j)
23020 continue
      continue
23018 continue
      continue
23004 continue
      continue
      return
      end
