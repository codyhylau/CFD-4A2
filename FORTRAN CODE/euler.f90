      program euler

      use common_block

      implicit none



! Local stuff
      integer i,j, sup
      real :: U, cloc, M, Mmax, Mmin, Mloc

!***********************************************************************    
      ! set which extensions to use, 1 = true, 0 = false
!***********************************************************************  
      rkut        = 1            ! runge kutta (class A)
      def_corr    = 1        ! deferred correction (class B)
      const_hstag = 1     ! constant stagnation enthalpy (class C)
      var_time    = 1        ! spatially variable time step (class C)
!***********************************************************************       

! Open files to store the convergence history. Plotting is done via a separate
! program. "euler.log" is for use by pltconv. "pltcv.csv" is for use by paraview.

      open(unit=3,file='euler.log')
      open(unit=31,file='pltcv.csv')

! "read_data": to read in the data on the duct and geometry and flow conditions.

      call read_data

! "generate_grid": to set up the grid coordinates, element areas and projected
! lengths of the sides of the elements.

      call generate_grid

! You can call subroutines "output*" here to plot out the grid you have generated.
! "output" writes "euler.csv" for paraview, "ouput_hg" write "euler.plt" for eulplt,
! "output_mat" writes "euler.mat" for matlab.

      call output(0)
      call output_hg(p,0)
      call output_mat(0)

! "check_grid": to check that the areas and projected lengths are correct.

      call check_grid

! "crude_guess" is what its name says. it enables you to
! start a calculation and obtain a solution but it will take longer than
! necessary. when your program is working you should replace it
! with "flow_guess" to obtain a better guess and a faster solution.

    ! call crude_guess
!*********************************************************************************
      ! SET FLOW GUESS / NEW GUESS / SUPERSONIC GUESS
    sup = 0
    if(title==' Mach 3 Corner ') then    ! TEST CASE 5
      print *, 'Calling new guess'
      call new_guess
      
    else if(title == 'supersonic wedge with M1 = 1.6 ') then    ! TEST CASE 4
      print *, 'Calling supersonic flow guess'
      call flow_guess_sup
      
      sup = 1  ! for using supersonic apply_bconds later
    else
      print *, 'Calling flow guess'    ! TEST CASES 0,1,2,3
      call flow_guess
    end if
!*********************************************************************************

! You can call "output" here to plot out your initial guess of
! the flow field.

      call output(1)
      call output_mat(1)
      if(title==' 180 degree bend  with bump and extension') var_time=0
      
! "set_timestep": to set the length of the timestep.
! initially this is a constant time step based on a conservative guess
! of the mach number.

      call set_timestep


!************************************************************************
!     start the time stepping do loop for "nsteps" loops.
!************************************************************************

      do nstep = 1, nsteps

        do i=1,ni
          do j=1,nj
            ro_start(i,j) = ro(i,j)
            if(const_hstag /= 1) then
              roe_start(i,j) = roe(i,j)
            end if
            rovx_start(i,j) = rovx(i,j)
            rovy_start(i,j) = rovy(i,j)
          end do
        end do



      if(rkut == 1) then
        do nrkut = 1,4
          frkut = 1.0/(5.0-nrkut)
          call set_others
          if(sup == 1) then
            call apply_bconds_sup
          else
            call apply_bconds
          end if
          call set_fluxes
          
          call sum_fluxes_rkut(fluxi_mass,fluxj_mass,delro  , ro_inc, frkut)
          if(const_hstag /= 1) then
            call sum_fluxes_rkut(fluxi_enth,fluxj_enth,delroe ,roe_inc, frkut)
          end if
          call sum_fluxes_rkut(fluxi_xmom,fluxj_xmom,delrovx,rovx_inc, frkut)
          call sum_fluxes_rkut(fluxi_ymom,fluxj_ymom,delrovy,rovy_inc, frkut)
      
! Update solution

          do i=1,ni
            do j=1,nj
              ro  (i,j) = ro_start  (i,j) + ro_inc  (i,j)
              if(const_hstag /= 1) then
                roe (i,j) = roe_start (i,j) + roe_inc (i,j)
              end if  
              rovx(i,j) = rovx_start(i,j) + rovx_inc(i,j)
              rovy(i,j) = rovy_start(i,j) + rovy_inc(i,j)
            end do
          end do


          if(def_corr == 1) then
            call smooth_corr(ro, corr_ro)
            call smooth_corr(rovx, corr_rovx)
            call smooth_corr(rovy, corr_rovy)
            if(const_hstag /= 1) then
              call smooth_corr(roe, corr_roe)
            end if
          else
            call smooth(ro)
            call smooth(rovx)
            call smooth(rovy)
            if(const_hstag /= 1) then
              call smooth(roe)
            end if
          end if
        end do    ! RUNGE KUTTA END

      else    ! if not using RUNGE KUTTA
        call set_others
        if(sup == 1) then
            call apply_bconds_sup
        else
            call apply_bconds
        end if
        call set_fluxes
        call sum_fluxes(fluxi_mass,fluxj_mass,delro  , ro_inc)
        if(const_hstag /= 1) then
          call sum_fluxes(fluxi_enth,fluxj_enth,delroe ,roe_inc)
        end if
        call sum_fluxes(fluxi_xmom,fluxj_xmom,delrovx,rovx_inc)
        call sum_fluxes(fluxi_ymom,fluxj_ymom,delrovy,rovy_inc)
        
        do i=1,ni
          do j=1,nj
            ro  (i,j) = ro_start  (i,j) + ro_inc  (i,j)
            if(const_hstag /= 1) then
              roe (i,j) = roe_start (i,j) + roe_inc (i,j)
            end if  
            rovx(i,j) = rovx_start(i,j) + rovx_inc(i,j)
            rovy(i,j) = rovy_start(i,j) + rovy_inc(i,j)
          end do
        end do
        
        if(def_corr == 1) then
          call smooth_corr(ro, corr_ro)
          call smooth_corr(rovx, corr_rovx)
          call smooth_corr(rovy, corr_rovy)
          if(const_hstag /= 1) then
            call smooth_corr(roe, corr_roe)
          end if
        else
          call smooth(ro)
          call smooth(rovx)
          call smooth(rovy)
          if(const_hstag /= 1) then
            call smooth(roe)
          end if
        end if
          
      end if
        


! Check convergence and write out summary every 5 steps


        if(mod(nstep,5)==0) then
          call check_conv          
        end if

        if(var_time == 1 .AND. mod(nstep,5)==0) then    ! update spatially variablt time step every 5 steps
        
          do i = 1,ni
            do j = 1,nj
              sound_speed(i,j) = sqrt(gamma*rgas*t(i,j))
            end do
          end do
        
          do i = 1,ni-1
            do j = 1,nj-1
              U = 0.25 * (vxy(i,j)+vxy(i,j+1)+vxy(i+1,j)+vxy(i+1,j+1))
              cloc = 0.25* (sound_speed(i,j)+sound_speed(i,j+1)+sound_speed(i+1,j)+sound_speed(i+1,j+1))
                    
              step(i,j) = cfl*dmin1(i,j) / (cloc + U)
         
            end do
          end do
         end if   
       
! Stop looping if converged to the input tolerance "conlim"

        if( emax < conlim .and.  eavg < (0.5*conlim) ) then
          write(6,*) ' Calculation converged in ',nstep,' iterations'
          write(6,*) ' To a convergence limit of ', conlim
          exit
        endif

      end do
      
      !call output(1)
      !call output_mat(1)

!************************************************************************
!  end of time stepping do loop for "nsteps" loops.
!************************************************************************

! Calculation finished. call "output" to write the plotting file.
! N.B. Solution hasn't necessarily converged.

      call output(1)
      call output_hg(p,1)
      call output_mat(1)
      
!*************************************************************************
!   PRINT OUT AVERAGE OUTLET MACH NUMBER FOR VERIFICATION
      M=0
      Mmax=0
      Mmin=10.
      do j=1,nj
        Mloc = vxy(ni,j)/sqrt(gamma*rgas*t(ni,j))
        if(Mloc>Mmax) Mmax= Mloc
        if(Mloc<Mmin) Mmin= Mloc
        M = M + Mloc
      end do
      M = M/nj
      print *, '**********************************************'
      print *, 'Average Outlet M = ', M
      print *, 'Mmax = ', Mmax
      print *, 'Mmin = ', Mmin
!*************************************************************************

      close(3)
      close(31)


      end
