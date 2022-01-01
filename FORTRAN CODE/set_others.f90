      subroutine set_others

      use common_block


      implicit none

! Local stuff
      integer  ::  i, j
! This routine calculates secondary flow variables from the primary ones
! at every grid point.

! The primary variables are ro, rovx, rovy and roe

! The secondary variables are the velocity components vx(i,j) and vy(i,j),
! the static pressure p(i,j) and the stagnation enthalpy hstag(i,j).
! Note:  "hstag"  not  "ho".

! INSERT your code here
    
    do i=1,ni
      do j=1,nj
        vx(i,j) = rovx(i,j)/ro(i,j)
        vy(i,j) = rovy(i,j)/ro(i,j)
        vxy(i,j) = sqrt(vx(i,j)**2.0+vy(i,j)**2.0)  
        t(i,j) = tstagin - vxy(i,j)**2.0/(2.0*cp)  
        
        if(const_hstag /= 1) then  
          ! p(i,j) = ro(i,j) * rgas * t(i,j) 
          p(i,j) = (gamma-1.0)*(roe(i,j)-0.5*ro(i,j)*vxy(i,j)**2.0)
          hstag(i,j) = (roe(i,j)+p(i,j))/ro(i,j)
        else                          ! constant stagnation enthalpy extension
          hstag(i,j) = cp * tstagin               
          p(i,j) = ro(i,j) * rgas * t(i,j)
        end if
        
      end do
    end do
    
      end
