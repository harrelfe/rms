## Usage: ratfor -o ../ormuv.f ormuv.r
## Computes the score vector and compressed information matrix for
## an ordinal regression model with possibly very many intercepts
## l= nx^2 + 2nx*kint + 3kint - 2
## If kint=1 just use a regular square matrix, l=p^2
SUBROUTINE ormuv(n, p, kint, nx, x, y, pr, fpa, fpb, fppa, fppb, u, v,
                 ja, ia, l, lia, kk)
IMPLICIT REAL*8 (a-h,o-z)
INTEGER p, y(n), ja(l), ia(lia), z, kk(p)
REAL*8 X(n,nx), pr(n), fpa(n), fpb(n), fppa(n), fppb(n), u(p), v(l), ld
do k=1,kint {
  uk = 0d0
  do j=1, n {
    z = y(j)
    a = 0d0
    if(z - 1 == k)  a = a + fpa(j)
    else if(z == k) a = a - fpb(j)
    uk = uk + a / pr(j)
  }
  u(k) = uk
}
if(nx > 0) do k = (kint + 1), p {
  uk = 0d0
  do j=1, n {
    uk = uk + (fpa(j) - fpb(j)) * x(j, k-kint) / pr(j)
  }
  u(k) = uk
}

iv = 0
do m = 1,p {
  if(kint > 1) {
    ## Compute column numbers for nonzero elements: kk
    if(m == 1) {
      nkk = 2
      kk(1) = 1
      kk(2) = 2
    } else if(m > 1 & m < kint) {
      nkk = 3
      kk(1) = m - 1
      kk(2) = m
      kk(3) = m + 1
    } else if(m == kint) {
      nkk = 2
      kk(1) = m-1
      kk(2) = m
    } else {
      nkk = kint
      do mm=1, kint {
        kk(mm) = mm
      }
    }
    do mm=(kint+1),p {
      nkk = nkk + 1
      kk(nkk) = mm
    }
  } else nkk = p
  # call intpr('nkk', 3, nkk, 1)
  # call intpr('kk',  2, kk, nkk)
  do ik = 1,nkk {
    if(kint == 1) k = ik
    else k = kk(ik)
    vmk = 0d0
    do j = 1,n {
      z = y(j)
      pa  = fpa(j)
      pb  = fpb(j)
      ppa = fppa(j)
      ppb = fppb(j)
      w = 1/(pr(j)*pr(j))
      if(m <= kint & k <= kint) a =
        - w * (pa*ld(z - 1 == m) - pb*ld(z == m)) *
          (pa*ld(z - 1 == k) - pb*ld(z == k)) +
            (ppa*ld(z - 1 == m)*ld(m == k) -
             ppb*ld(z == m)*ld(m == k))/pr(j)
      else if(m > kint & k > kint) a =
        x(j,m-kint) * x(j,k-kint) / pr(j) *
          (-1/pr(j) * (pa - pb) * (pa - pb) + ppa - ppb)
      else {
        mi = max(m, k)
        ki = min(m, k)
        a = x(j, mi - kint) / pr(j) *
          (-1/pr(j) * (pa - pb) * (pa*ld(z - 1 == ki) - pb*ld(z == ki)) +
           ppa*ld(z - 1 == ki) - ppb*ld(z == ki))
      }
      vmk = vmk + a
    }
    iv = iv + 1
    v(iv) = - vmk
    if(kint > 1) {
      ja(iv) = k
      if(ik == 1) ia(m) = iv
    }
  }
}
if(kint > 1) ia(p+1) = iv + 1
return
end

real*8 function ld(a)
logical a
if(a) ld = 1d0
else ld = 0d0
return
end
