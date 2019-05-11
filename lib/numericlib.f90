module numericlib
  use kind_parameters,ONLY:&
       sp
  implicit none
contains
  
  
REAL(KIND=SP) FUNCTION LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
  use kind_parameters,ONLY:&
       sp
  REAL(kind=sp)        LAMS,PHIS,POLPHI,POLLAM
  REAL(kind=sp) ZRPI18,ZPIR18,ZSINPOL,ZCOSPOL,ZLAMPOL,ZPHIS,ZLAMS
  REAL(kind=sp) ZARG1,ZARG2
  DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
  
  ZSINPOL = SIN(ZPIR18*POLPHI)
  ZCOSPOL = COS(ZPIR18*POLPHI)
  ZLAMPOL = ZPIR18*POLLAM
  ZPHIS   = ZPIR18*PHIS
  ZLAMS   = LAMS
  IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
  ZLAMS   = ZPIR18*ZLAMS
  
  ZARG1   = SIN(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +&
       ZCOSPOL*           SIN(ZPHIS)) -&
       COS(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
  ZARG2   = COS(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +&
       ZCOSPOL*           SIN(ZPHIS)) +&
       SIN(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
  IF (ABS(ZARG2).LT.1.E-30) THEN
     IF (ABS(ZARG1).LT.1.E-30) THEN
        LMSTOLM =   0.0
     ELSEIF (ZARG1.GT.0.) THEN
        LMSTOLM =  90.0
     ELSE
        LMSTOLM = -90.0
     ENDIF
  ELSE
     LMSTOLM = ZRPI18*ATAN2(ZARG1,ZARG2)
  ENDIF
  
  RETURN
END FUNCTION LMSTOLM


REAL(kind=sp) FUNCTION PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
  use kind_parameters,ONLY:&
       sp

  REAL(KIND=SP)        LAMS,PHIS,POLPHI,POLLAM
  REAL(KIND=SP) ZRPI18,ZPIR18,SINPOL,COSPOL,ZPHIS,ZLAMS,ARG
  DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
  
  SINPOL = SIN(ZPIR18*POLPHI)
  COSPOL = COS(ZPIR18*POLPHI)
  ZPHIS  = ZPIR18*PHIS
  ZLAMS  = LAMS
  IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
  ZLAMS  = ZPIR18*ZLAMS
  ARG     = COSPOL*COS(ZPHIS)*COS(ZLAMS) + SINPOL*SIN(ZPHIS)
  
  PHSTOPH = ZRPI18*ASIN(ARG)
  
  RETURN
END FUNCTION PHSTOPH


 SUBROUTINE spline(x,y,n,yp1,ypn,y2)
!  use,intrinsic::ieee_arithmetic,only:IEEE_SELECTED_REAL_KIND
!  use,intrinsic::ieee_exceptions
   use ieee_arithmetic
   use ieee_features
   use kind_parameters,ONLY:&
        sp
   implicit none
!   integer,parameter:: sp=IEEE_SELECTED_REAL_KIND(6,37)
   logical :: underflow_support,gradual,underflow
   INTEGER n,NMAX
   logical :: should_halt, was_flagged
   
   REAL(kind=sp) yp1,ypn,x(n),y(n),y2(n)
   PARAMETER (NMAX=500)
   INTEGER i,k
   REAL(kind=sp) p,qn,sig,un,u(NMAX)


   call ieee_set_underflow_mode(.false.)
   call ieee_get_underflow_mode(gradual)
   
   if (.not.gradual) then 
      !     print *,'Able to set abrupt underflow mode'
   else
      !     stop 'error setting underflow mode'
   end if
   

   
   if (yp1.gt..99e30) then
      y2(1)=0.
      u(1)=0.
   else
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
   endif
   do  i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+&
           1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*&
           u(i-1))/p
      call ieee_get_flag(ieee_underflow,underflow) 
  end do
  if (ypn.gt..99e30) then
     qn=0.
     un=0.
  else
     qn=0.5
     un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do  k=n-1,1,-1
     y2(k)=y2(k)*y2(k+1)+u(k)
  end do

END SUBROUTINE spline
 
SUBROUTINE splint(xa,ya,y2a,n,x,y)
  
  INTEGER n
  REAL(kind=sp) x,y,xa(n),y2a(n),ya(n)
  INTEGER k,khi,klo
  REAL(kind=sp) a,b,h
  klo=1
  khi=n
1 if (khi-klo.gt.1) then
     k=(khi+klo)/2
     if(xa(k).gt.x)then
        khi=k
     else
        klo=k
     endif
     goto 1
  endif
  h=xa(khi)-xa(klo)
  if (h.eq.0.) pause 'bad xa input in splint'
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  return
END SUBROUTINE splint


!   ---------------------------------------------------------

  
!   -----------------------------------------------------------------
!   Rene's KINK function (for smoothing at bounadries)
!   -----------------------------------------------------------------

real(kind=sp) function kink (x,a)
  
  implicit none

  !   declaration of parameters
  real(kind=sp)   x,a
  
      !   parameters
  real(kind=sp)   pi
  parameter  (pi=3.1415926535)
  
  if (x.lt.0.) then
     kink=0.
  elseif (x.gt.a) then
     kink=1.
  else
     kink=(1.+cos(pi*(x-a)/a))/2.
  endif
  
  return
end function kink

!   -----------------------------------------------------------------
!   Rene's KINK function (for smoothing at bounadries)
!   -----------------------------------------------------------------

subroutine mdvfill (out,inp,flag,nx,ny,nz,maxiter)
  
  implicit none
  
  !   Declaration of subroutine parameters
  integer nx,ny,nz
  real(kind=sp)    inp (nx,ny,nz)
  real(kind=sp)    out (nx,ny,nz)
  integer flag(nx,ny,nz)
  integer maxiter
  
  !   Parameters
  real(kind=sp)             omega        ! Omega fopr SOR
  parameter        (omega=1.5)
  
  !   Auxiliary variables
  integer i,j,k
  integer iter
  real(kind=sp)    tmp0(nx,ny,nz)
  real(kind=sp)    tmp1(nx,ny,nz)
  integer il,ir,ju,jd
  integer count
  real(kind=sp)    mean(nz)
  
  !   Calculate mean of variable for all levels
  do k=1,nz
     mean(k) = 0.
     count   = 0
     do i=1,nx
        do j=1,ny
           if ( flag(i,j,k).eq.0 ) then
              count   = count + 1
              mean(k) = mean(k) + inp(i,j,k)
           endif
        enddo
     enddo
     if ( count.ne.0 ) then
        mean(k) = mean(k)/real(count)
     else
        mean(k) = 0.
     endif
  enddo

  !   Create first guess
  do k=1,nz
     do i=1,nx
        do j=1,ny
           if ( flag(i,j,k).eq.0 ) then
              tmp0(i,j,k) = inp(i,j,k)
           else
              tmp0(i,j,k) = mean(k)
        endif
     enddo
  enddo
enddo
   
!   SOR iterations
iter = 0
122 continue

!   Loop over all points
do k=1,nz
   do i=1,nx
      do j=1,ny
         
         !       Apply the updating only for specified points
          if ( flag(i,j,k).ne.1 ) goto 121
          
          !       Get neighbouring points - no handling of periodic domains!
          il = i-1
          if (il.lt.1) il=1
          ir = i+1
          if ( ir.gt.nx ) ir=nx
          jd = j-1
          if (jd.lt.1) jd=1
          ju = j+1
          if ( ju.gt.ny ) ju=ny

!       Update field
          tmp1(i,j,k) = 0.25 * ( tmp0(il,j,k) + tmp0(ir,j,k) +&
               tmp0(i,ju,k) + tmp0(i,jd,k) )
          
          tmp0(i,j,k) = omega * tmp1(i,j,k) +&
               (1. - omega) * tmp1(i,j,k)
          
          !       Exit point for loop
 121      continue
          
       enddo
    enddo
 enddo
 
!   Decide whether further iterations are needed
 iter = iter + 1
 if ( iter.lt.maxiter) goto 122
 
 !   Save output
 do i=1,nx
    do j=1,ny
       do k=1,nz
          if ( flag(i,j,k).eq.1 ) then
             out(i,j,k) = tmp0(i,j,k)
          endif
       enddo
    enddo
 enddo
 
end subroutine mdvfill



end module numericlib


