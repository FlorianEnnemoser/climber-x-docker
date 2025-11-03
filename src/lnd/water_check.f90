!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : w a t e r _ c h e c k _ m o d
!
!  Purpose : water conservation checks
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
!
! This file is part of CLIMBER-X.
!
! CLIMBER-X is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! CLIMBER-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with CLIMBER-X.  If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module water_check_mod

   use precision, only : wp
   use lnd_grid, only : nsurf, nveg, nsoil, nl, is_veg, is_ice, is_lake, i_ice, i_lake
   use lnd_params, only : dt

   implicit none

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  w a t e r _ b a l a n c e
  !   Purpose    :  compute soil temperature 
  !              :  by solving the tridiagonal system
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine water_check(i,j,frac_surf,f_veg,rain,snow,et,runoff_sur,calving,drainage,icemelt, &
                          w_w,w_i,w_snow,w_can,s_can,w_w_old,w_i_old,w_snow_old,w_can_old,s_can_old,lake_water_tendency, &
                          water_cons)

    implicit none

    real(wp), intent(in) :: f_veg
    real(wp), dimension(:), intent(in) :: rain, snow
    real(wp), dimension(:), intent(in) :: runoff_sur, calving, drainage, icemelt
    real(wp), dimension(:), intent(in) :: frac_surf, et, w_can, w_can_old, s_can, s_can_old
    real(wp), dimension(:), intent(in) :: w_w, w_i, w_w_old, w_i_old
    real(wp), dimension(:), intent(in) :: w_snow, w_snow_old
    real(wp), intent(in) :: lake_water_tendency
    real(wp), dimension(:), intent(out) :: water_cons

    integer :: n,i,j
    real(wp) :: et_sum, dw_can, rain_sum, snow_sum


    water_cons = 0._wp

    if( f_veg .gt. 0._wp ) then

     rain_sum = 0._wp
     snow_sum = 0._wp
     et_sum = 0._wp
     dw_can = 0._wp
     do n=1,nveg
      rain_sum = rain_sum + rain(n) * frac_surf(n)/f_veg
      snow_sum = snow_sum + snow(n) * frac_surf(n)/f_veg
      et_sum = et_sum + et(n) * frac_surf(n)/f_veg
      dw_can = dw_can + (w_can(n) - w_can_old(n) + s_can(n) - s_can_old(n)) * frac_surf(n)/f_veg
     enddo

     water_cons(is_veg) = rain_sum*dt + snow_sum*dt &
                       - runoff_sur(is_veg)*dt - calving(is_veg)*dt - drainage(is_veg)*dt &
                       - et_sum*dt - dw_can &
                       - sum(w_w-w_w_old + w_i-w_i_old) &
                       - (w_snow(is_veg)-w_snow_old(is_veg))

    !  if( drainage(is_veg)*dt.lt.-0.1_wp ) then
    !    print *,''
    !    print *,'negative drainage',drainage(is_veg)*dt
    !    print *,'i,j',i,j
    !    print *,'water balance over SOIL',water_cons(is_veg) !,t_skin_veg
    !    print *,'fveg',f_veg
    !    print *,'sum(frac_surf)',sum(frac_surf)
    !    print *,'frac_surf',frac_surf
    !    print *,'dw_soil',sum(w_w-w_w_old + w_i-w_i_old)
    !    print *,'dw_snow_if',(w_snow(is_veg)-w_snow_old(is_veg))
    !    print *,'dw_can',dw_can
    !    print *,'rain+snow,rain,snow',rain_sum*dt + snow_sum*dt,rain_sum*dt,snow_sum*dt
    !    print *,'run_w,run_i',runoff_sur(is_veg)*dt,calving(is_veg)*dt
    !    print *,'drain',drainage(is_veg)*dt
    !    print *,'et',et_sum*dt
    !    print *,'et(n)',et*dt
    !    print *,'w_w old',w_w_old
    !    print *,'w_w final',w_w
    !    print *,'w_i old',w_i_old
    !    print *,'w_i final',w_i
    !    print *,'w_snow, w_snow_old if',w_snow(is_veg),w_snow_old(is_veg)
    !    print *,'w_snow, w_snow_old i',w_snow(is_ice),w_snow_old(is_ice)
    !    print *,'w_can_old',w_can_old
    !    print *,'w_can',w_can
    !    print *,'s_can_old',s_can_old
    !    print *,'s_can',s_can
    !  endif

      !if(abs(water_cons(is_veg)).gt.1.d-10 ) then
      ! if(abs(water_cons(is_veg)).gt.1.d-3 ) then
      !  !print *,'mask_snow,i,j',mask_snow,i,j
      !  print *,''
      !  print *,'i,j',i,j
      !  print *,'water balance over SOIL',water_cons(is_veg) !,t_skin_veg
      !  print *,'fveg',f_veg
      !  print *,'sum(frac_surf)',sum(frac_surf)
      !  print *,'frac_sur',frac_surf
      !  print *,'dw_soil',sum(w_w-w_w_old + w_i-w_i_old)
      !  print *,'dw_snow_if',(w_snow(is_veg)-w_snow_old(is_veg))
      !  print *,'dw_can',dw_can
      !  print *,'rain+snow,rain,snow',rain_sum*dt + snow_sum*dt,rain_sum*dt,snow_sum*dt
      !  print *,'run_w,run_i',runoff_sur(is_veg)*dt,calving(is_veg)*dt
      !  print *,'drain',drainage(is_veg)*dt
      !  print *,'et',et_sum*dt
      !  print *,'w_w old',w_w_old
      !  print *,'w_w final',w_w
      !  print *,'w_i old',w_i_old
      !  print *,'w_i final',w_i
      !  print *,'w_snow, w_snow_old if',w_snow(is_veg),w_snow_old(is_veg)
      !  print *,'w_snow, w_snow_old i',w_snow(is_ice),w_snow_old(is_ice)
      !  print *,'w_can_old',w_can_old
      !  print *,'w_can',w_can
      !  print *,'s_can_old',s_can_old
      !  print *,'s_can',s_can
      !  if(abs(water_cons(is_veg)).gt.0.1 ) stop
      ! endif

    endif 

    if( frac_surf(i_ice) .gt. 0._wp ) then

     water_cons(is_ice) = rain(i_ice)*dt + snow(i_ice)*dt + icemelt(is_ice)*dt &
                       - runoff_sur(is_ice)*dt - calving(is_ice)*dt - drainage(is_ice)*dt &
                       - et(i_ice)*dt &
                       - (w_snow(is_ice)-w_snow_old(is_ice))

      ! if(abs(water_cons(is_ice)).gt.1.d-10) then
       !print *,'mask_snow,i,j',mask_snow,i,j
      !  print *,''
      !  print *,'water balance over ICE',water_cons(is_ice) !,t_skin_veg
       !print *,'sum(frac_surf)',sum(frac_surf)
       !print *,'frac_sur',frac_surf
       !print *,'dw_snow_i',(w_snow(is_ice)-w_snow_old(is_ice))
       !print *,'rain+snow,rain,snow',rain*dt + snow*dt,rain*dt,snow*dt
       !print *,'run_w,run_i',runoff_sur(is_ice)*dt,calving(is_ice)*dt
       !print *,'drain',drainage(is_ice)*dt
       !print *,'et',et(i_ice)*dt
       !print *,'w_snow, w_snow_old i',w_snow(is_ice),w_snow_old(is_ice)
       !stop
      ! endif

    endif

! lake water is not conserved    
!    if( frac_surf(i_lake) .gt. 0._wp ) then
!
!     water_cons(is_lake) = rain*dt + snow*dt &
!                         - runoff_sur(is_lake)*dt - calving(is_lake)*dt - drainage(is_lake)*dt &
!                         - et(i_lake)*dt &
!                         - (w_snow(is_lake)-w_snow_old(is_lake)) &
!                         - lake_water_tendency*dt
!
!      if (water_cons(is_lake).gt.1.e-10_wp .or. water_cons(is_lake).lt.-1.e-10_wp) then ! avoid warning message because of ice melt!
!       !print *,'mask_snow,i,j',mask_snow,i,j
!       print *,''
!       print *,'water balance over LAKE',water_cons(is_lake) !,t_skin_veg
!       print *,'lake_water_tendency',lake_water_tendency
!       print *,'sum(frac_surf)',sum(frac_surf)
!       print *,'frac_surf',frac_surf
!       print *,'rain+snow,rain,snow',rain*dt + snow*dt,rain*dt,snow*dt
!       print *,'run_w,run_i',runoff_sur(is_lake)*dt,calving(is_lake)*dt
!       print *,'drain',drainage(is_lake)*dt
!       print *,'et',et(i_lake)*dt
!       print *,'dw_snow',(w_snow(is_lake)-w_snow_old(is_lake))
!       print *,'w_snow, w_snow_old i',w_snow(is_lake),w_snow_old(is_lake)
!       stop
!      endif
!
!    endif


    return

  end subroutine water_check

end module water_check_mod

