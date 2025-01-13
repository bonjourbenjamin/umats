! common utilities for umats

module utils
   implicit none
   private
   public gradsh_2d, gradtemp_2d, gradsh_3d, gradtemp_3d

   contains

   pure function matinv2(a) result(b)
      ! performs a direct calculation of the inverse of a 2×2 matrix.
      ! reference https://fortranwiki.org/fortran/show/matrix+inversion

      real(kind=8), intent(in) :: a(2,2)   !! matrix
      real(kind=8)             :: b(2,2)   !! inverse matrix
      real(kind=8)             :: detinv

      ! calculate the inverse determinant of the matrix
      detinv = 1/(a(1,1)*a(2,2) - a(1,2)*a(2,1))

      ! calculate the inverse of the matrix
      b(1,1) = +detinv * a(2,2)
      b(2,1) = -detinv * a(2,1)
      b(1,2) = -detinv * a(1,2)
      b(2,2) = +detinv * a(1,1)

   end function matinv2

   pure function matinv3(a) result(b)
      ! performs a direct calculation of the inverse of a 3×3 matrix.
      ! reference https://fortranwiki.org/fortran/show/matrix+inversion

      real(kind=8), intent(in) :: a(3,3)   !! matrix
      real(kind=8)             :: b(3,3)   !! inverse matrix
      real(kind=8)             :: detinv

      ! calculate the inverse determinant of the matrix
      detinv = 1/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2)&
         - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
         + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

      ! calculate the inverse of the matrix
      b(1,1) = +detinv * (a(2,2)*a(3,3) - a(2,3)*a(3,2))
      b(2,1) = -detinv * (a(2,1)*a(3,3) - a(2,3)*a(3,1))
      b(3,1) = +detinv * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
      b(1,2) = -detinv * (a(1,2)*a(3,3) - a(1,3)*a(3,2))
      b(2,2) = +detinv * (a(1,1)*a(3,3) - a(1,3)*a(3,1))
      b(3,2) = -detinv * (a(1,1)*a(3,2) - a(1,2)*a(3,1))
      b(1,3) = +detinv * (a(1,2)*a(2,3) - a(1,3)*a(2,2))
      b(2,3) = -detinv * (a(1,1)*a(2,3) - a(1,3)*a(2,1))
      b(3,3) = +detinv * (a(1,1)*a(2,2) - a(1,2)*a(2,1))

   end function matinv3

   subroutine gradsh_2d(noel, coort, sht, gradsh)
      ! compute the gradient of the hydrostatic stress in 2d element (CPE8RT)

      implicit none

      integer, intent(in) :: noel
      real(kind=8), intent(in) :: coort(:,:,:), sht(:,:)
      real(kind=8), intent(inout) :: gradsh(:,:,:)

      integer :: k1
      real(kind=8) :: h, g, a1, a2, a3, a4, b1, b2, b3, b4
      real(kind=8), dimension(2,2) :: xjacm, xjaci
      real(kind=8), dimension(2,4) :: deriv

      do k1=1,4 ! hard coded for 4 integration points
         if (k1==1) then
            g=-1.d0
            h=-1.d0
         elseif (k1==2) then
            g= 1.d0
            h=-1.d0
         elseif (k1==3) then
            g=-1.d0
            h= 1.d0
         elseif (k1==4) then
            g=1.d0
            h=1.d0
         end if

         ! d(ni)/d(g), shape function derivatives w.r.t natural coordinates g-direction
         deriv(1,1)=-(1.d0/4.0)*(1-h)
         deriv(1,2)= (1.d0/4.0)*(1-h)
         deriv(1,3)=-(1.d0/4.0)*(1+h)
         deriv(1,4)= (1.d0/4.0)*(1+h)

         ! d(ni)/g(h), shape function derivatives w.r.t natural coordinates h-direction
         deriv(2,1)=-(1.d0/4.0)*(1-g)
         deriv(2,2)=-(1.d0/4.0)*(1+g)
         deriv(2,3)= (1.d0/4.0)*(1-g)
         deriv(2,4)= (1.d0/4.0)*(1+g)

         ! j, jacobian matrix
         xjacm(1,1)=deriv(1,1)*coort(noel,1,1)+deriv(1,2)*coort(noel,2,1)+deriv(1,3)*coort(noel,3,1)+deriv(1,4)*coort(noel,4,1)
         xjacm(1,2)=deriv(1,1)*coort(noel,1,2)+deriv(1,2)*coort(noel,2,2)+deriv(1,3)*coort(noel,3,2)+deriv(1,4)*coort(noel,4,2)
         xjacm(2,1)=deriv(2,1)*coort(noel,1,1)+deriv(2,2)*coort(noel,2,1)+deriv(2,3)*coort(noel,3,1)+deriv(2,4)*coort(noel,4,1)
         xjacm(2,2)=deriv(2,1)*coort(noel,1,2)+deriv(2,2)*coort(noel,2,2)+deriv(2,3)*coort(noel,3,2)+deriv(2,4)*coort(noel,4,2)

         ! j^-1, inverse of the jacobian matrix
         xjaci = matinv2(xjacm)

         !djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
         !xjaci(1,1)=xjacm(2,2)/djacb
         !xjaci(1,2)=-xjacm(1,2)/djacb
         !xjaci(2,1)=-xjacm(2,1)/djacb
         !xjaci(2,2)=xjacm(1,1)/djacb

         ! d(ni)/d(x), shape function derivatives w.r.t global coordinates x-direction
         a1=xjaci(1,1)*deriv(1,1)+xjaci(1,2)*deriv(2,1)
         a2=xjaci(1,1)*deriv(1,2)+xjaci(1,2)*deriv(2,2)
         a3=xjaci(1,1)*deriv(1,3)+xjaci(1,2)*deriv(2,3)
         a4=xjaci(1,1)*deriv(1,4)+xjaci(1,2)*deriv(2,4)

         ! d(ni)/d(y), shape function derivatives w.r.t global coordinates y-direction
         b1=xjaci(2,1)*deriv(1,1)+xjaci(2,2)*deriv(2,1)
         b2=xjaci(2,1)*deriv(1,2)+xjaci(2,2)*deriv(2,2)
         b3=xjaci(2,1)*deriv(1,3)+xjaci(2,2)*deriv(2,3)
         b4=xjaci(2,1)*deriv(1,4)+xjaci(2,2)*deriv(2,4)

         ! calculate gradient of hydrostatic stress
         ! gradsh x-direction
         gradsh(noel,k1,1)=a1*sht(noel,1)+a2*sht(noel,2)+a3*sht(noel,3)+a4*sht(noel,4)

         ! gradsh y-direction
         gradsh(noel,k1,2)=b1*sht(noel,1)+b2*sht(noel,2)+b3*sht(noel,3)+b4*sht(noel,4)

      end do

   end subroutine gradsh_2d


   subroutine gradtemp_2d(noel, coort, temp, gradtemp)
      ! compute the gradient of temperature in 2d element (CPE8RT)

      implicit none

      integer, intent(in) :: noel
      real(kind=8), intent(in) :: coort(:,:,:), temp(:,:)
      real(kind=8), intent(inout) :: gradtemp(:,:,:)

      ! use same calculation as gradiant of hydrostatic stress
      call gradsh_2d(noel, coort, temp, gradtemp)

   end subroutine gradtemp_2d


   subroutine gradsh_3d(noel, coort, sht, gradsh)
      ! compute the gradient of the hydrostatic stress in 3d element (C3D20RT)

      implicit none

      integer, intent(in) :: noel
      real(kind=8), intent(in) :: coort(:,:,:), sht(:,:)
      real(kind=8), intent(inout) :: gradsh(:,:,:)

      integer :: k1
      real(kind=8) :: h, g, r
      real(kind=8) :: a1, a2, a3, a4, a5, a6, a7, a8
      real(kind=8) :: b1, b2, b3, b4, b5, b6, b7, b8
      real(kind=8) :: c1, c2, c3, c4, c5, c6, c7, c8
      real(kind=8), dimension(3,3) :: xjacm, xjaci
      real(kind=8), dimension(3,8) :: deriv

      do k1=1,8 ! hard coded for 8 integration points
         if (k1==1) then
            g=-1.d0
            h=-1.d0
            r=-1.d0
         elseif (k1==2) then
            g= 1.d0
            h=-1.d0
            r=-1.d0
         elseif (k1==3) then
            g= 1.d0
            h= 1.d0
            r=-1.d0
         elseif (k1==4) then
            g=-1.d0
            h= 1.d0
            r=-1.d0
         elseif (k1==5) then
            g=-1.d0
            h=-1.d0
            r= 1.d0
         elseif (k1==6) then
            g= 1.d0
            h=-1.d0
            r= 1.d0
         elseif (k1==7) then
            g=1.d0
            h=1.d0
            r=1.d0
         elseif (k1==8) then
            g=-1.d0
            h= 1.d0
            r= 1.d0
         end if

         ! d(ni)/d(g), shape function derivatives w.r.t natural coordinates g-direction
         deriv(1,1)=-(1.d0/8.0)*(1-h)*(1-r)
         deriv(1,2)= (1.d0/8.0)*(1-h)*(1-r)
         deriv(1,3)= (1.d0/8.0)*(1+h)*(1-r)
         deriv(1,4)=-(1.d0/8.0)*(1+h)*(1-r)
         deriv(1,5)=-(1.d0/8.0)*(1-h)*(1+r)
         deriv(1,6)= (1.d0/8.0)*(1-h)*(1+r)
         deriv(1,7)= (1.d0/8.0)*(1+h)*(1+r)
         deriv(1,8)=-(1.d0/8.0)*(1+h)*(1+r)

         ! d(ni)/g(h), shape function derivatives w.r.t natural coordinates h-direction
         deriv(2,2)=-(1.d0/8.0)*(1-g)*(1-r)
         deriv(2,3)=-(1.d0/8.0)*(1+g)*(1-r)
         deriv(2,1)= (1.d0/8.0)*(1+g)*(1-r)
         deriv(2,4)= (1.d0/8.0)*(1-g)*(1-r)
         deriv(2,5)=-(1.d0/8.0)*(1-g)*(1+r)
         deriv(2,6)=-(1.d0/8.0)*(1+g)*(1+r)
         deriv(2,7)= (1.d0/8.0)*(1+g)*(1+r)
         deriv(2,8)= (1.d0/8.0)*(1-g)*(1+r)

         ! d(ni)/g(h), shape function derivatives w.r.t natural coordinates r-direction
         deriv(3,2)=-(1.d0/8.0)*(1-g)*(1-h)
         deriv(3,3)=-(1.d0/8.0)*(1+g)*(1-h)
         deriv(3,1)=-(1.d0/8.0)*(1+g)*(1+h)
         deriv(3,4)=-(1.d0/8.0)*(1-g)*(1+h)
         deriv(3,5)= (1.d0/8.0)*(1-g)*(1-h)
         deriv(3,6)= (1.d0/8.0)*(1+g)*(1-h)
         deriv(3,7)= (1.d0/8.0)*(1+g)*(1+h)
         deriv(3,8)= (1.d0/8.0)*(1-g)*(1+h)

         ! j, jacobian matrix
         xjacm(1,1)=deriv(1,1)*coort(noel,1,1)+deriv(1,2)*coort(noel,2,1)+deriv(1,3)*coort(noel,3,1)+deriv(1,4)*coort(noel,4,1) &
            +deriv(1,5)*coort(noel,5,1)+deriv(1,6)*coort(noel,6,1)+deriv(1,7)*coort(noel,7,1)+deriv(1,8)*coort(noel,8,1)

         xjacm(1,2)=deriv(1,1)*coort(noel,1,2)+deriv(1,2)*coort(noel,2,2)+deriv(1,3)*coort(noel,3,2)+deriv(1,4)*coort(noel,4,2) &
            +deriv(1,5)*coort(noel,5,2)+deriv(1,6)*coort(noel,6,2)+deriv(1,7)*coort(noel,7,2)+deriv(1,8)*coort(noel,8,2)

         xjacm(1,3)=deriv(1,1)*coort(noel,1,3)+deriv(1,2)*coort(noel,2,3)+deriv(1,3)*coort(noel,3,3)+deriv(1,4)*coort(noel,4,3) &
            +deriv(1,5)*coort(noel,5,3)+deriv(1,6)*coort(noel,6,3)+deriv(1,7)*coort(noel,7,3)+deriv(1,8)*coort(noel,8,3)

         xjacm(2,1)=deriv(2,1)*coort(noel,1,1)+deriv(2,2)*coort(noel,2,1)+deriv(2,3)*coort(noel,3,1)+deriv(2,4)*coort(noel,4,1) &
            +deriv(2,5)*coort(noel,5,1)+deriv(2,6)*coort(noel,6,1)+deriv(2,7)*coort(noel,7,1)+deriv(2,8)*coort(noel,8,1)

         xjacm(2,2)=deriv(2,1)*coort(noel,1,2)+deriv(2,2)*coort(noel,2,2)+deriv(2,3)*coort(noel,3,2)+deriv(2,4)*coort(noel,4,2) &
            +deriv(2,5)*coort(noel,5,2)+deriv(2,6)*coort(noel,6,2)+deriv(2,7)*coort(noel,7,2)+deriv(2,8)*coort(noel,8,2)

         xjacm(2,3)=deriv(2,1)*coort(noel,1,3)+deriv(2,2)*coort(noel,2,3)+deriv(2,3)*coort(noel,3,3)+deriv(2,4)*coort(noel,4,3) &
            +deriv(2,5)*coort(noel,5,3)+deriv(2,6)*coort(noel,6,3)+deriv(2,7)*coort(noel,7,3)+deriv(2,8)*coort(noel,8,3)

         xjacm(3,1)=deriv(3,1)*coort(noel,1,1)+deriv(3,2)*coort(noel,2,1)+deriv(3,3)*coort(noel,3,1)+deriv(3,4)*coort(noel,4,1) &
            +deriv(3,5)*coort(noel,5,1)+deriv(3,6)*coort(noel,6,1)+deriv(3,7)*coort(noel,7,1)+deriv(3,8)*coort(noel,8,1)

         xjacm(3,2)=deriv(3,1)*coort(noel,1,2)+deriv(3,2)*coort(noel,2,2)+deriv(3,3)*coort(noel,3,2)+deriv(3,4)*coort(noel,4,2) &
            +deriv(3,5)*coort(noel,5,2)+deriv(3,6)*coort(noel,6,2)+deriv(3,7)*coort(noel,7,2)+deriv(3,8)*coort(noel,8,2)

         xjacm(3,3)=deriv(3,1)*coort(noel,1,3)+deriv(3,2)*coort(noel,2,3)+deriv(3,3)*coort(noel,3,3)+deriv(3,4)*coort(noel,4,3) &
            +deriv(3,5)*coort(noel,5,3)+deriv(3,6)*coort(noel,6,3)+deriv(3,7)*coort(noel,7,3)+deriv(3,8)*coort(noel,8,3)

         ! j^-1, inverse of the jacobian matrix
         xjaci = matinv3(xjacm)

         ! d(n)/d(x), shape function derivatives w.r.t global coordinates x-direction
         a1=xjaci(1,1)*deriv(1,1)+xjaci(1,2)*deriv(2,1)+xjaci(1,3)*deriv(3,1)
         a2=xjaci(1,1)*deriv(1,2)+xjaci(1,2)*deriv(2,2)+xjaci(1,3)*deriv(3,2)
         a3=xjaci(1,1)*deriv(1,3)+xjaci(1,2)*deriv(2,3)+xjaci(1,3)*deriv(3,3)
         a4=xjaci(1,1)*deriv(1,4)+xjaci(1,2)*deriv(2,4)+xjaci(1,3)*deriv(3,4)
         a5=xjaci(1,1)*deriv(1,5)+xjaci(1,2)*deriv(2,5)+xjaci(1,3)*deriv(3,5)
         a6=xjaci(1,1)*deriv(1,6)+xjaci(1,2)*deriv(2,6)+xjaci(1,3)*deriv(3,6)
         a7=xjaci(1,1)*deriv(1,7)+xjaci(1,2)*deriv(2,7)+xjaci(1,3)*deriv(3,7)
         a8=xjaci(1,1)*deriv(1,8)+xjaci(1,2)*deriv(2,8)+xjaci(1,3)*deriv(3,8)

         ! d(n)/d(y), shape function derivatives w.r.t global coordinates y-direction
         b1=xjaci(2,1)*deriv(1,1)+xjaci(2,2)*deriv(2,1)+xjaci(2,3)*deriv(3,1)
         b2=xjaci(2,1)*deriv(1,2)+xjaci(2,2)*deriv(2,2)+xjaci(2,3)*deriv(3,2)
         b3=xjaci(2,1)*deriv(1,3)+xjaci(2,2)*deriv(2,3)+xjaci(2,3)*deriv(3,3)
         b4=xjaci(2,1)*deriv(1,4)+xjaci(2,2)*deriv(2,4)+xjaci(2,3)*deriv(3,4)
         b5=xjaci(2,1)*deriv(1,5)+xjaci(2,2)*deriv(2,5)+xjaci(2,3)*deriv(3,5)
         b6=xjaci(2,1)*deriv(1,6)+xjaci(2,2)*deriv(2,6)+xjaci(2,3)*deriv(3,6)
         b7=xjaci(2,1)*deriv(1,7)+xjaci(2,2)*deriv(2,7)+xjaci(2,3)*deriv(3,7)
         b8=xjaci(2,1)*deriv(1,8)+xjaci(2,2)*deriv(2,8)+xjaci(2,3)*deriv(3,8)

         ! d(n)/d(z), shape function derivatives w.r.t global coordinates z-direction
         c1=xjaci(3,1)*deriv(1,1)+xjaci(3,2)*deriv(2,1)+xjaci(3,3)*deriv(3,1)
         c2=xjaci(3,1)*deriv(1,2)+xjaci(3,2)*deriv(2,2)+xjaci(3,3)*deriv(3,2)
         c3=xjaci(3,1)*deriv(1,3)+xjaci(3,2)*deriv(2,3)+xjaci(3,3)*deriv(3,3)
         c4=xjaci(3,1)*deriv(1,4)+xjaci(3,2)*deriv(2,4)+xjaci(3,3)*deriv(3,4)
         c5=xjaci(3,1)*deriv(1,5)+xjaci(3,2)*deriv(2,5)+xjaci(3,3)*deriv(3,5)
         c6=xjaci(3,1)*deriv(1,6)+xjaci(3,2)*deriv(2,6)+xjaci(3,3)*deriv(3,6)
         c7=xjaci(3,1)*deriv(1,7)+xjaci(3,2)*deriv(2,7)+xjaci(3,3)*deriv(3,7)
         c8=xjaci(3,1)*deriv(1,8)+xjaci(3,2)*deriv(2,8)+xjaci(3,3)*deriv(3,8)

         ! calculate gradient of hydrostatic stress
         ! gradsh x-direction
         gradsh(noel,k1,1)=a1*sht(noel,1)+a2*sht(noel,2)+a3*sht(noel,3)+a4*sht(noel,4)+a5*sht(noel,5)+a6*sht(noel,6)+a7*sht(noel,7) &
            +a8*sht(noel,8)

         ! gradsh y-direction
         gradsh(noel,k1,2)=b1*sht(noel,1)+b2*sht(noel,2)+b3*sht(noel,3)+b4*sht(noel,4)+b5*sht(noel,5)+b6*sht(noel,6)+b7*sht(noel,7) &
            +b8*sht(noel,8)

         ! gradsh z-direction
         gradsh(noel,k1,3)=c1*sht(noel,1)+c2*sht(noel,2)+c3*sht(noel,3)+c4*sht(noel,4)+c5*sht(noel,5)+c6*sht(noel,6)+c7*sht(noel,7) &
            +c8*sht(noel,8)

      end do

   end subroutine gradsh_3d


   subroutine gradtemp_3d(noel, coort, temp, gradtemp)
      ! compute the gradient of temperature in 3d element (C3D20RT)

      implicit none

      integer, intent(in) :: noel
      real(kind=8), intent(in) :: coort(:,:,:), temp(:,:)
      real(kind=8), intent(inout) :: gradtemp(:,:,:)

      ! use same calculation as gradiant of hydrostatic stress
      call gradsh_3d(noel, coort, temp, gradtemp)

   end subroutine gradtemp_3d


end module utils
