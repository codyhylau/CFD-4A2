      subroutine set_timestep

      use common_block


      implicit none

! Local stuff
      integer :: i,j
      real  ::  umax, U, cloc
      real, dimension(ni,nj) :: vel, tstat, c

! Tthis subroutine sets the length of the time step based on the
! stagnation speed of sound "astag" and the minimum length scale
! of any element, "dmin". The timestep must be called "deltat"

! An assumption that the maximum flow speed will be equal to "astag"
! is also made. This will be pessimistic for subsonic flows
! but may be optimistic for supersonic flows. In the latter case the
! length of the time step as determined by "cfl" may need to be reduced.

! The cfl number was input as data in data set "flow"
  
  if(var_time /= 1) then
      astag  = sqrt(gamma*rgas*tstagin)
      umax   = astag
      deltat = cfl * dmin/(umax+astag)
  else  
    do i = 1,ni
      do j = 1,nj
        vel(i,j) = sqrt(vx(i,j)**2.0 + vy(i,j)**2.0)
        tstat(i,j) = tstagin - vel(i,j)**2.0/(2.0*cp)
        c(i,j) = sqrt(gamma*rgas*tstat(i,j))
      end do
    end do
    
    do i = 1,ni-1
      do j = 1,nj-1
        U = 0.25 * (vel(i,j)+vel(i,j+1)+vel(i+1,j)+vel(i+1,j+1))
        cloc = 0.25* (c(i,j)+c(i,j+1)+c(i+1,j)+c(i+1,j+1))
              
        step(i,j) = cfl*dmin1(i,j) / (cloc + U)
    
      end do
    end do
  end if
end