C Output from Public domain Ratfor, version 1.01
      subroutine ormuv(n, p, kint, nx, x, y, pr, fpa, fpb, fppa, fppb, u
     *, v, ja, ia, l, lia, kk)
      implicit real*8 (a-h,o-z)
      integer p, y(n), ja(l), ia(lia), z, kk(p)
      real*8 x(n,nx), pr(n), fpa(n), fpb(n), fppa(n), fppb(n), u(p), v(l
     *), ld
      do23000 k=1,kint 
      uk = 0d0
      do23002 j=1, n 
      z = y(j)
      a = 0d0
      if(z - 1 .eq. k)then
      a = a + fpa(j)
      else
      if(z .eq. k)then
      a = a - fpb(j)
      endif
      endif
      uk = uk + a / pr(j)
23002 continue
23003 continue
      u(k) = uk
23000 continue
23001 continue
      if(nx .gt. 0)then
      do23010 k = (kint + 1), p 
      uk = 0d0
      do23012 j=1, n 
      uk = uk + (fpa(j) - fpb(j)) * x(j, k-kint) / pr(j)
23012 continue
23013 continue
      u(k) = uk
23010 continue
23011 continue
      endif
      iv = 0
      do23014 m = 1,p 
      if(kint .gt. 1)then
      if(m .eq. 1)then
      nkk = 2
      kk(1) = 1
      kk(2) = 2
      else
      if(m .gt. 1 .and. m .lt. kint)then
      nkk = 3
      kk(1) = m - 1
      kk(2) = m
      kk(3) = m + 1
      else
      if(m .eq. kint)then
      nkk = 2
      kk(1) = m-1
      kk(2) = m
      else
      nkk = kint
      do23024 mm=1, kint 
      kk(mm) = mm
23024 continue
23025 continue
      endif
      endif
      endif
      do23026 mm=(kint+1),p 
      nkk = nkk + 1
      kk(nkk) = mm
23026 continue
23027 continue
      else
      nkk = p
      endif
      do23028 ik = 1,nkk 
      if(kint .eq. 1)then
      k = ik
      else
      k = kk(ik)
      endif
      vmk = 0d0
      do23032 j = 1,n 
      z = y(j)
      pa = fpa(j)
      pb = fpb(j)
      ppa = fppa(j)
      ppb = fppb(j)
      w = 1/(pr(j)*pr(j))
      if(m .le. kint .and. k .le. kint)then
      a = - w * (pa*ld(z - 1 .eq. m) - pb*ld(z .eq. m)) * (pa*ld(z - 1 .
     *eq. k) - pb*ld(z .eq. k)) + (ppa*ld(z - 1 .eq. m)*ld(m .eq. k) - p
     *pb*ld(z .eq. m)*ld(m .eq. k))/pr(j)
      else
      if(m .gt. kint .and. k .gt. kint)then
      a = x(j,m-kint) * x(j,k-kint) / pr(j) * (-1/pr(j) * (pa - pb) * (p
     *a - pb) + ppa - ppb)
      else
      mi = max(m, k)
      ki = min(m, k)
      a = x(j, mi - kint) / pr(j) * (-1/pr(j) * (pa - pb) * (pa*ld(z - 1
     * .eq. ki) - pb*ld(z .eq. ki)) + ppa*ld(z - 1 .eq. ki) - ppb*ld(z .
     *eq. ki))
      endif
      endif
      vmk = vmk + a
23032 continue
23033 continue
      iv = iv + 1
      v(iv) = - vmk
      if(kint .gt. 1)then
      ja(iv) = k
      if(ik .eq. 1)then
      ia(m) = iv
      endif
      endif
23028 continue
23029 continue
23014 continue
23015 continue
      if(kint .gt. 1)then
      ia(p+1) = iv + 1
      endif
      return
      end
      real*8 function ld(a)
      logical a
      if(a)then
      ld = 1d0
      else
      ld = 0d0
      endif
      return
      end
