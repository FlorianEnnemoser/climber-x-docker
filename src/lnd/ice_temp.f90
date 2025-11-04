!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i c e _ t e m p _ m o d
!
!  Purpose : ice sheet temperature
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
module ice_temp_mod

  use precision, only : wp
  use constants, only : Lf, T0, cap_i
  use control, only : check_energy, check_water
  use lnd_grid, only : z, dz, nl, rdz_pos, rdz_neg
  use lnd_params, only : dt, rdt, hydro_par 
  use tridiag, only : tridiag_solve

  implicit none

  private
  public :: ice_temp

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i c e _ t e m p
  !   Purpose    :  compute ice temperature
  !              :  by solving the tridiagonal system
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_temp(mask_snow,h_snow,cap_ice,lambda_int,flx_g,dflxg_dT,flx_melt, &
                     t_ice,w_snow,snowmelt,icemelt,t_ice_old,w_snow_old, &
                     energy_cons_ice,i,j)

    implicit none

    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: h_snow
    real(wp), intent(in) :: flx_g, dflxg_dT, flx_melt
    real(wp), dimension(0:), intent(in) :: cap_ice, lambda_int
    real(wp), intent(inout) :: w_snow
    real(wp), dimension(0:), intent(inout) :: t_ice
    real(wp), intent(out) :: snowmelt, icemelt
    real(wp), intent(inout) :: w_snow_old
    real(wp), dimension(0:), intent(out) :: t_ice_old
    real(wp), intent(out) :: energy_cons_ice

    integer :: i,j
    integer :: k
    real(wp), dimension(0:nl) :: a, b, c, r, x
    real(wp), dimension(0:nl) :: cap, rcap
    real(wp), dimension(0:nl) :: dz_loc, z_loc, rdz_loc_pos, rdz_loc_neg
    real(wp) :: w_ice, H, w_melt, H_star, dw_snow, dw_ice, w_snow_tmp, energy_warm_snow, energy_warm_ice, flx_excess
    real(wp) :: i_p, j_p, n_p
    logical :: i_print

    i_print = .false.
    n_p = 1
    i_p = 55
    j_p = 31

    ! set snowmelt to zero
    snowmelt = 0._wp
    icemelt = 0._wp
    energy_warm_snow = 0._wp
    energy_warm_ice = 0._wp

    dz_loc(0) = h_snow
    z_loc(0) = -0.5_wp * dz_loc(0)
    dz_loc(1:nl) = dz(1:nl)
    z_loc(1:nl) = z(1:nl)

    rdz_loc_pos(1:nl-1) = rdz_pos(1:nl-1)
    rdz_loc_neg(2:nl) = rdz_neg(2:nl)
    rdz_loc_pos(0) = 1._wp / (z_loc(1) - z_loc(0))
    rdz_loc_neg(1) = 1._wp / (z_loc(1) - z_loc(0))

    ! assign some arrays
    cap = cap_ice*dz_loc
    rcap(1:nl) = 1._wp/cap(1:nl)
    if( cap(0) .gt. 0._wp ) rcap(0) = 1._wp/cap(0)

    ! save old prognostic variables
    t_ice_old = t_ice
    w_snow_old = w_snow


    if( mask_snow .eq. 1 ) then

      ! top layer, snow, k=0
      k = 0
      a(k) = 0._wp
      b(k) = 1._wp + dt * rcap(k) * ( lambda_int(k) * rdz_loc_pos(k) - dflxg_dT )
      c(k) = - dt * rcap(k) * lambda_int(k) * rdz_loc_pos(k)
      r(k) = t_ice(k) + dt * rcap(k) * ( flx_g - dflxg_dT*t_ice(k) )

      ! intermediate layers
      do k=1,nl-1
        a(k) = - dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
        b(k) = 1._wp + dt * rcap(k) * ( lambda_int(k-1) * rdz_loc_neg(k) + lambda_int(k) * rdz_loc_pos(k) )
        c(k) = - dt * rcap(k) * lambda_int(k) * rdz_loc_pos(k)
        r(k) = t_ice(k)
      enddo

      ! bottom layer, k=nl
      k = nl
      a(k) = - dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
      b(k) = 1._wp + dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
      c(k) = 0._wp
      r(k) = t_ice(k)

      ! solve tridiagonal system
      call tridiag_solve(a,b,c,r,x,nl+1)
      ! assign new soil temperature
      t_ice = x

      ! PHASE CHANGE

      k = 0

      w_ice = w_snow

      ! freezing or melting occuring
      if( (t_ice(k) .gt. T0 .and. w_ice .gt. 0._wp) ) then
        H = dflxg_dT * (T0 - t_ice(k)) - cap(k)*rdt*(T0 -t_ice(k))
        w_melt = H*dt/Lf
        ! melting
        if( w_melt .gt. 0._wp ) then
          snowmelt = min(w_ice, w_melt) * rdt  ! kg/m2/s
          w_snow = max(0._wp, w_ice-w_melt)
          ! determine energy not used during melt and adjust temperature if necessary
          ! add excess snow energy to first ice layer
          t_ice(k) = T0
          H_star = H - Lf * (w_ice - w_snow) * rdt
          t_ice(1) = t_ice(1) + dt*rcap(1) * H_star
          if(i_print .and. k.eq.0) print *,'H,w_melt,H_star snow layer'
          if(i_print .and. k.eq.0) print *,H,w_melt,H_star
        endif
      endif

      ! use flx_melt to directly melt part of the snow layer 
      ! without the need for the whole snow layer to reach melting point
      if (flx_melt.gt.0._wp) then
        dw_snow = flx_melt / (rdt*((T0-t_ice(0))*cap_i+Lf))
        if (dw_snow.le.w_snow) then
          ! enough snow to melt
          snowmelt = snowmelt + dw_snow*rdt  ! kg/m2/s
          w_snow = w_snow-dw_snow
          energy_warm_snow = rdt*(T0-t_ice(0))*cap_i*dw_snow ! W/m2
        else
          ! not enough snow to melt
          ! melt all
          snowmelt = snowmelt + w_snow*rdt  ! kg/m2/s
          ! determine energy not used during melt and adjust temperature
          flx_excess = flx_melt + dflxg_dT*(T0-t_ice(0)) - w_snow*rdt*Lf - rdt*(T0-t_ice(0))*cap_i*w_snow  ! W/m2
          w_snow = 0._wp
          t_ice(0) = T0
          t_ice(1) = t_ice(1) + dt*rcap(1) * flx_excess
          ! distribute the excess energy to all soil layers
          !t_ice(1:nl) = t_ice(1:nl) + dt*rcap(1:nl) * flx_excess * dz_loc(1:nl) / sum(dz_loc(1:nl))
        endif
      endif

      ! ice melt
      H = 0._wp
      do k=1,nl
        if (t_ice(k).gt.T0) then
          H = H - cap(k)*rdt*(T0-t_ice(k))
          t_ice(k) = T0
        endif
      enddo
      icemelt = H/Lf ! kg/m2/s
      
    else ! w_snow lower than w_snow_crit, no explicit snow layer

      ! top layer, k=1
      k = 1
      a(k) = 0._wp
      b(k) = 1._wp + dt * rcap(k) * ( lambda_int(k) * rdz_loc_pos(k) - dflxg_dT )
      c(k) = - dt * rcap(k) * lambda_int(k) * rdz_loc_pos(k)
      r(k) = t_ice(k) + dt * rcap(k) * ( flx_g - dflxg_dT*t_ice(k) )

      ! intermediate layers
      do k=2,nl-1
        a(k) = - dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
        b(k) = 1._wp + dt * rcap(k) * ( lambda_int(k-1) * rdz_loc_neg(k) + lambda_int(k) * rdz_loc_pos(k) )
        c(k) = - dt * rcap(k) * lambda_int(k) * rdz_loc_pos(k)
        r(k) = t_ice(k) 
      enddo

      ! bottom layer, k=nl
      k = nl
      a(k) = - dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
      b(k) = 1._wp + dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
      c(k) = 0._wp
      r(k) = t_ice(k) 

      ! solve tridiagonal system
      call tridiag_solve(a(1:nl),b(1:nl),c(1:nl),r(1:nl),x(1:nl),nl)
      ! assign new soil temperature
      t_ice(1:nl) = x(1:nl)

      ! top layer ice melt
      k = 1
      if (t_ice(k).gt.T0) then
        H = dflxg_dT * (T0-t_ice(k)) - cap(k)*rdt*(T0-t_ice(k))
        icemelt = H/Lf
        t_ice(k) = T0
      endif

      ! use flx_melt to directly melt (not-explicit) snow (if present)
      flx_excess = 0._wp
      if (flx_melt.gt.0._wp) then
        dw_snow = flx_melt / (rdt*Lf) ! snow that can be melted
        if (dw_snow.le.w_snow) then
          ! enough snow to melt
          snowmelt = snowmelt + dw_snow*rdt  ! kg/m2/s
          w_snow = w_snow-dw_snow
          flx_excess = 0._wp
        else
          ! not enough snow to melt
          ! melt all
          snowmelt = snowmelt + w_snow*rdt  ! kg/m2/s
          ! determine energy not used during melt and adjust temperature
          flx_excess = flx_melt - w_snow*rdt*Lf ! W/m2
          w_snow = 0._wp
        endif
      endif
      ! use flx_excess to directly melt ice 
      ! without the need for the whole ice layer to reach melting point
      if (flx_excess.gt.0._wp) then
        dw_ice = flx_excess / (rdt*(max(0._wp,T0-t_ice(1))*cap_i+Lf)) ! kg/m2, ice that can be melted
        icemelt = icemelt + dw_ice*rdt  ! kg/m2/s
        energy_warm_ice = rdt*max(0._wp,T0-t_ice(1))*cap_i*dw_ice ! W/m2
      endif

      ! additionally melt ice if temperature of any layer above freezing point
      H = 0._wp
      do k=1,nl
        if (t_ice(k).gt.T0) then
          H = H - cap(k)*rdt*(T0-t_ice(k))
          t_ice(k) = T0
        endif
      enddo
      icemelt = icemelt + H/Lf ! kg/m2/s

      ! if snow present but no explicit snow layer, use possible first ice layer excess energy to melt snow
      if( w_snow.gt.0._wp ) then
        H = dflxg_dT * (T0 - t_ice(1)) - cap(1)*rdt*(T0 - t_ice(1))

        if( H .gt. 0._wp ) then

          w_melt = H*dt/Lf ! kg/m2
          w_snow_tmp = w_snow
          w_snow = max( 0._wp, w_snow-w_melt )  ! kg/m2
          H_star = H - Lf*rdt * (w_snow_tmp - w_snow)

          t_ice(1) = T0 + dt*rcap(1) * H_star / (1._wp - dt*rcap(1) * dflxg_dT)

          snowmelt = snowmelt + (w_snow_tmp - w_snow) * rdt  ! kg/m2/s

        endif
      endif

      ! set temperature of non existing snow layer equal to first ice layer
      t_ice(0) = min(T0, t_ice(1))

    endif

    if( check_energy ) then
      ! ice energy balance
      if( mask_snow .eq. 1 ) then
        energy_cons_ice =  flx_g + dflxg_dT*(t_ice(0)-t_ice_old(0)) + flx_melt &
        - Lf*(snowmelt+icemelt) - energy_warm_snow - energy_warm_ice &
        - sum(cap*rdt * (t_ice-t_ice_old))
      else
        energy_cons_ice = flx_g + dflxg_dT*(t_ice(1)-t_ice_old(1)) + flx_melt &
        - Lf*(snowmelt+icemelt) - energy_warm_snow - energy_warm_ice &
        - sum(cap(1:nl)*rdt * (t_ice(1:nl)-t_ice_old(1:nl)))
      endif


      ! if(abs(energy_cons_ice).gt.1.d-10 ) then
        ! print *,' '
        ! print *,'energy conservation ICE',energy_cons_ice
        ! print *,'i,j,mask_snow',i,j,mask_snow
        ! print *,'t_soil_old',t_ice_old
        ! print *,'t_soil_after_cond',x
        ! print *,'t_soil_new',t_ice
        ! if(mask_snow.eq.1) then
        !   print *,'g,dg_dT',flx_g,dflxg_dT*(t_ice(0)-t_ice_old(0)),flx_g+dflxg_dT*(t_ice(0)-t_ice_old(0))
        ! else
        !   print *,'g,dg_dT',flx_g,dflxg_dT*(t_ice(1)-t_ice_old(1)),flx_g+dflxg_dT*(t_ice(1)-t_ice_old(1))
        ! endif
        ! print *,'flx_melt',flx_melt
        ! print *,'energy_warm_snow',energy_warm_snow
        ! print *,'energy_warm_ice',energy_warm_ice
        ! print *,'energy_internal',cap*rdt * (t_ice-t_ice_old)
        ! print *,'snowmelt energy',Lf*snowmelt
        ! print *,'icemelt energy',Lf*icemelt
        ! print *,'snow_new,snow_old,snowmelt',w_snow,w_snow_old,snowmelt*dt
        ! print *,'cap',cap
        ! print *,'lambda_int',lambda_int
        ! print *,' '
        ! stop
      ! endif
    endif

    ! if( check_water ) then
    !   ! ice water balance
    !   if( abs(w_snow - w_snow_old+snowmelt*dt) .gt. 1.d-10 ) then
    !     print *,'water conservation ',snowmelt*dt+w_snow - w_snow_old 
    !     print *,w_snow,w_snow_old,snowmelt
    !     stop
    !   endif
    ! endif


    ! check temperature range
    !     if(maxval(t_ice) .gt. 320._wp ) print *,'WARNING t_ice > 20°C ',t_ice

    if(maxval(t_ice) .gt. 350._wp .or. minval(t_ice) .lt. 100._wp) then
      print *,''
      print *,'WARNING t_ice out of range',t_ice
      print *,'i,j',i,j
      print *,'masksnow',mask_snow
      k = 0
      print *,'SNOW layer'
      print *,'mask_snow,snowmelt',mask_snow, snowmelt*dt
      print *,'t_snow, before, after',t_ice_old(k),t_ice(k)
      print *,'w_snow, before, after',w_snow_old,w_snow

      k=1
      print *,'top ice layer'
      print *,'t_ice, before, after',t_ice_old(k),t_ice(k)
      print *,'g,dg_dT,dg_dT*dT',flx_g,dflxg_dT,dflxg_dT*(t_ice(1)-t_ice_old(1))
      print *,'cap,lambda',cap, lambda_int
      print *,''
      !if(i.eq.54.and.j.eq.25) stop
    endif

    if(i_print .and. i.eq.54.and.j.eq.25) then
      print *,'WARNING t_soil out of range',t_ice
      print *,'i,j',i,j
      print *,'masksnow',mask_snow
      k = 0
      print *,'SNOW layer'
      print *,'mask_snow,snowmelt',mask_snow, snowmelt*dt
      print *,'t_snow, before, after',t_ice_old(k),t_ice(k)
      print *,'w_snow, before, after',w_snow_old,w_snow

      k=1
      print *,'top ice layer'
      print *,'t_ice, before, after',t_ice_old(k),t_ice(k)
      print *,'g,dg_dT',flx_g,dflxg_dT*(t_ice(1)-t_ice_old(1))
      print *,''
      print *,'cap,lambda',cap, lambda_int
      print *,''
    endif


     return

  end subroutine ice_temp


end module ice_temp_mod

