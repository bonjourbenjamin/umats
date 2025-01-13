program name
   !implicit none
   !integer :: a = 10
   !real :: b = 1.2344

   !print *, "hello world", a, b
   use utils

   ! 2D
   ! integer, parameter :: noels=1
   ! integer, parameter :: nopts=4
   ! integer, parameter :: ndims=2

   ! integer :: noel

   ! real(kind=8), dimension(noels,nopts) :: ShT           ! hydrostatic stress at each integration point
   ! real(kind=8), dimension(noels,nopts,ndims) :: coorT  ! integration point corrds
   ! real(kind=8), dimension(noels,nopts,ndims) :: gradSh ! graident of hydrostatic stress

   ! ! define these to use subroutine
   ! noel=1

   ! ShT(noel,:)=(/2.d0,  1.d0, 1.d0, 2.d0/)

   ! coorT(noel,:,1)=(/ 0.d0,  1.d0,  1.d0, 0.d0/)
   ! coorT(noel,:,2)=(/ 0.d0,  0.d0,  1.d0, 1.d0/)

   ! gradSh=0.d0

   ! call gradSh_2D(noel,coorT,ShT,gradSh)

   ! call gradTemp_2D(noel,coorT,ShT,gradSh)

   ! print *, gradSh(:,:,1)
   ! print *, gradSh(:,:,2)

   ! 3D
   integer, parameter :: noels=1
   integer, parameter :: nopts=8
   integer, parameter :: ndims=3

   integer :: noel

   real(kind=8), dimension(noels,nopts) :: ShT           ! hydrostatic stress at each integration point
   real(kind=8), dimension(noels,nopts,ndims) :: coorT  ! integration point corrds
   real(kind=8), dimension(noels,nopts,ndims) :: gradSh ! graident of hydrostatic stress

   ! define these to use subroutine
   noel=1

   ShT(noel,:)=(/2.d0, 1.d0, 1.d0, 1.d0, 1.d0, 2.d0, 1.d0, 2.d0/)

   coorT(noel,:,1)=(/-1.d0,  1.d0,  1.d0, -1.d0, -1.d0,  1.d0,  1.d0, -1.d0/)
   coorT(noel,:,2)=(/-1.d0, -1.d0,  1.d0,  1.d0, -1.d0, -1.d0,  1.d0,  1.d0/)
   coorT(noel,:,3)=(/-1.d0, -1.d0, -1.d0, -1.d0,  1.d0,  1.d0,  1.d0,  1.d0/)

   gradSh=0.d0

   call gradSh_3D(noel,coorT,ShT,gradSh)

   call gradTemp_3D(noel,coorT,ShT,gradSh)

   print *, gradSh(:,:,1)
   print *, gradSh(:,:,2)
   print *, gradSh(:,:,3)

end program name
