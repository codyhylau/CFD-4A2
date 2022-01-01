      subroutine check_grid


      use common_block


      implicit none

! Local stuff
      integer :: i, j
      real :: xSum, ySum

! Check your grid and areas for both the "bump" and the "bend"
! test data.

! First check that all areas are positive (by program or writing out)

! INSERT your code here

  do i=1,ni-1
    do j=1,nj-1
      if(area(i,j)<0) then
        write(6,*) 'Area not positive!'
        stop
      end if
    end do
  end do  

! Next check that the sum of the length vectors of the 4 faces
! of every element is very nearly zero in each coordinate direction.
! It is absolutely essential that this is correct !
! If not go back and check your subroutine "generate_grid".

! Careful with a test of the form
!          if( a == 0.0 ) then .....
! This will probably never be true.  Computers work to a finite number of
! Significant figures and "a" will probably be +0.0000001 or -0.0000001.
! Test for something like
!          if( abs(a) <= small_number ) then ...

! Insert your code here

  ! Compare area for case 0 with hand calculation:
  ! write(6,*)
  ! write(6,*) 'Top left area =', area(1,nj-1), ' Hand Calculated =', 3.0/19.0*10.0/59.0
  ! write(6,*) 'Top right area =', area(nj-1,nj-1), ' Hand Calculated =', 3.0/19.0*10.0/59.0
  ! write(6,*)
  
  do i=1,ni-1
    do j=1,nj-1
    
        xSum = dlix(i,j)+dljx(i,j)-dlix(i+1,j)-dljx(i,j+1) ! Sum of vectors in x
        ySum = dliy(i,j)+dljy(i,j)-dliy(i+1,j)-dljy(i,j+1) ! Sum of vectors in y
      
      if(abs(xSum)>0.0000001) then
        write(6,*) 'Grid vector sum in x not zero! Sum is: ', xSum, 'Coordinates = ', i,', ', j
        stop
      end if
      if(abs(ySum)>0.0000001) then
        write(6,*) 'Grid vector sum in y not zero! Sum is: ', ySum, 'Coordinates = ', i,', ', j
        stop
      end if
      
    end do
  end do

! Any other tests that you can think of. For example you could plot
! contours of  "area(i,j)" by using  --  call output_hg(area) .
!  call output_hg(area,1)

      end
