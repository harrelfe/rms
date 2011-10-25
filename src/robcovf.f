      subroutine robcovf(n, p, nc, start, len, u, s, v, w)
      implicit DOUBLE PRECISION (a-h,o-z)
      integer p, start(p), len(p)
      DOUBLE PRECISION u(n,p), s(p), v(p,p), w(p,p)
      do 23000 i=1,p
      do 23002 j=1,p
      w(i,j)=0d0
23002 continue
23000 continue
      do 23004 k=1,nc
      do 23006 i=1,p
      s(i)=0d0
      do 23008 j=1,p
      v(i,j)=0d0
23008 continue
23006 continue
      do 23010 i=start(k),start(k)+len(k)-1
      do 23012 j=1,p
      s(j)=s(j)+u(i,j)
23012 continue
23010 continue
      do 23014 i=1,p
      do 23016 j=1,p
      v(i,j)=v(i,j)+s(i)*s(j)
23016 continue
23014 continue
      do 23018 i=1,p
      do 23020 j=1,p
      w(i,j)=w(i,j)+v(i,j)
23020 continue
23018 continue
23004 continue
      return
      end

