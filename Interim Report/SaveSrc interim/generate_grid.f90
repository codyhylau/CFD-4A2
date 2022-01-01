      subroutine generate_grid

      use common_block


      implicit none

! Local variables
      integer :: i, j
      real :: a1,a2,b1,b2, L1, L2  

! Calculate x and y values  x(i,j),y(i,j) of the grid nodes.

! For each value of "i" the i-grid line joins (xlow(i),ylow(i)) to
! (xhigh(i),yhigh(i)). for each value of "i" grid points (nodes) should be
! linearly interpolated between these values for j between 1 and nj.
! i.e.  x(i,1) should be xlow(i), x(i,nj) should be xhigh(i), etc.

! INSERT your code here

  do i=1,ni
    do j=1,nj
    
      ! Linear interpolate x coordinates using number of j divisions on each i-line, starting from xlow
      x(i,j) = xlow(i) + (j-1) * (xhigh(i)-xlow(i)) / (nj-1) 
      ! Linear interpolate y coordinates using number of j divisions on each i-line, starting from ylow
      y(i,j) = ylow(i) + (j-1) * (yhigh(i)-ylow(i)) / (nj-1)
      
    end do
  end do  

! Calculate the areas of the cells area(i,j)
! (N.B. there are (ni-1) x (nj-1) cells.

! The area of a quadrilateral (regular or irregular) can be shown to be
! half of the cross product of the vectors forming the diagonals.
! see Hirsch volume 1, section 6.2.1. (or lecture).
! Make sure that the area comes out positive!

! INSERT your code here

  do i=1,ni-1
    do j=1,nj-1
      
      a1 = x(i,j+1) - x(i+1,j)
      a2 = x(i+1,j+1)-x(i,j)
      b1 = y(i,j+1) - y(i+1,j)
      b2 = y(i+1,j+1) - y(i,j)
      area(i,j) = abs(a2*b1-a1*b2)*0.5  ! Area calculation using 0.5 VD1 cross VD2
      
    end do
  end do

! Calculate the x and y components of the length vector of the i-faces
! (i.e. those corresponding to i = constant).
! The length vector of a face is a vector normal to the face with
! magnitude equal to the length of the face.
! It is positive in the direction of an inward normal to the cell i,j .
! Call these lengths dlix(i,j) and dliy(i,j)

! INSERT your code here

  dmin = 100.0
  do i=1,ni
    do j=1,nj-1
    
      dlix(i,j) = y(i,j+1)-y(i,j)
      dliy(i,j) = x(i,j)-x(i,j+1)
      
      ! Calculate magnitude of cell edge lengths and set dmin if smaller
      L1 = sqrt(dlix(i,j)**2.0+dliy(i,j)**2.0)
    !  if(i==1.AND.j==1) dmin = L1
      if(L1<dmin) dmin = L1

    end do
  end do


! Now calculate the x and y components of the length vector of the j-faces. (i.e. those corresponding to j = constant)
! Call these lengths dljx(i,j) and dljy(i,j)

! INSERT your code here

  do i = 1,ni-1
    do j = 1,nj
      dljx(i,j) = y(i,j)-y(i+1,j)
      dljy(i,j) = x(i+1,j)-x(i,j)
      
      ! Calculate magnitude of cell edge lengths and set dmin if smaller
      L2 = sqrt(dljx(i,j)**2.0+dljy(i,j)**2.0)
      if(L2<dmin) dmin = L2
      
    end do
  end do
  
! Now find "dmin" the minimum length scale of any element. This is
! defined as the length of the shortest side of the element.
! Call this minimum "dmin". it is used to set the time step from the cfl no.

! Insert your code here (or in the do loops above).
    ! done in above DO loops

      write(6,*)
      write(6,*)  'Overall minimum element size = ', dmin
      write(6,*)
      
      end
