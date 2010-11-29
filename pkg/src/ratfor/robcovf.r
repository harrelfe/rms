# Computes sum of (within cluster sum of U)(within cluster sum of U)'
#
SUBROUTINE robcovf(n, p, nc, start, len, u, s, v, w)
IMPLICIT REAL*8 (a-h,o-z)
INTEGER p, start(p), len(p)
REAL*8 u(n,p), s(p), v(p,p), w(p,p)

do i=1,p		{
   do j=1,p		{
   w(i,j)=0d0		}}

do k=1,nc					{
   do i=1,p		{
     s(i)=0d0
     do j=1,p		{
     v(i,j)=0d0		}}
   do i=start(k),start(k)+len(k)-1	{
     do j=1,p		{
       s(j)=s(j)+u(i,j)	}
					}
   do i=1,p		{
     do j=1,p		{
     v(i,j)=v(i,j)+s(i)*s(j)
			}}

   do i=1,p		{
     do j=1,p		{
       w(i,j)=w(i,j)+v(i,j)
			}}
						}

return
end

