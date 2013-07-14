C Output from Public domain Ratfor, version 1.01
      subroutine robcovf(n, p, nc, start, len, u, s, v, w)
      implicit real*8 (a-h,o-z)
      integer p, start(nc), len(nc)
      real*8 u(n,p), s(p), v(p,p), w(p,p)
      do23000 i=1,p 
      do23002 j=1,p 
      w(i,j)=0d0 
23002 continue
23003 continue
23000 continue
23001 continue
      do23004 k=1,nc 
      do23006 i=1,p 
      s(i)=0d0
      do23008 j=1,p 
      v(i,j)=0d0 
23008 continue
23009 continue
23006 continue
23007 continue
      do23010 i=start(k),start(k)+len(k)-1 
      do23012 j=1,p 
      s(j)=s(j)+u(i,j) 
23012 continue
23013 continue
23010 continue
23011 continue
      do23014 i=1,p 
      do23016 j=1,p 
      v(i,j)=v(i,j)+s(i)*s(j)
23016 continue
23017 continue
23014 continue
23015 continue
      do23018 i=1,p 
      do23020 j=1,p 
      w(i,j)=w(i,j)+v(i,j)
23020 continue
23021 continue
23018 continue
23019 continue
23004 continue
23005 continue
      return
      end
