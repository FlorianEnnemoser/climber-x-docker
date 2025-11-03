!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : e b a l _ o c n _ m o d
!
!  Purpose : energy balance over open ocean
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
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
module ebal_ocn_mod

   use precision, only : wp
   use control, only: check_energy
   use constants, only : Le, Lf, q_sat_w
   use constants, only : cap_a, cap_w, rho_a, rho_w, emis_w, sigma, T0
   use sic_grid, only : h1
   use sic_params, only : rho_sic, dt

   implicit none

   private
   public :: ebal_ocn

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e b a l _ o c n
  !   Purpose    :  compute surface energy fluxes
  !              :  over ocean water by solving the surface energy balance equation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ebal_ocn(t_ocn, t_freeze, t_air, qair, pressure, swnet, lwdown, snow, rain, wind, Cde, Cdh, &
                      t_skin, &
                      flx_ocn, flx_sh, flx_lwu, flx_lh, fw_ocn, evp, dh_sic, i, j)

    implicit none

    real(wp), intent(in   ) :: t_ocn, t_freeze, t_air, qair, swnet
    real(wp), intent(in   ) :: pressure, lwdown, snow, rain, wind, Cde, Cdh
    real(wp), intent(inout) :: t_skin
    real(wp), intent(out  ) :: flx_ocn, flx_sh, flx_lwu, flx_lh, fw_ocn, evp
    real(wp), intent(out  ) :: dh_sic 

    integer :: i, j
    real(wp) :: sh, lh, g, lw, flux_freeze
    real(wp) :: emiss, rhoa
    real(wp) :: f_sh, f_e, f_lh, qsat
    real(wp) :: energy_cons


    ! top ocean layer temperature is used to compute the ocean-atmosphere fluxes
    t_skin = t_ocn

    emiss = emis_w
    qsat = q_sat_w(t_skin,pressure)

    rhoa = rho_a(t_air,pressure)
    f_sh = Cdh*wind*rhoa*cap_a
    f_lh = Le*Cde*wind*rhoa
    f_e  = Cde*wind*rhoa

    sh = f_sh * (t_skin - t_air)
    lh = f_lh * (qsat  - qair)
    lw = (1._wp-emiss)*lwdown + emiss*sigma * t_skin**4 
    g  = snow * Lf

    ! new virtual top ocean layer temperature
    t_skin = t_skin + (swnet + lwdown - lw - sh - lh - snow*Lf)*dt/(rho_w*cap_w*h1)

    dh_sic = 0._wp
    flux_freeze = 0._wp

    ! limit t_skin to be above freezing point.
    ! form sea ice if t_skin < t_freeze
    if (t_skin .lt. t_freeze) then
      ! reset t_skin to t_freeze

      ! energy released from sea ice formation, positive
      flux_freeze = rho_w*cap_w*h1 * (t_freeze-t_skin)/dt

      ! reset temperature
      t_skin = t_freeze

      ! thickness of sea ice which is forming
      dh_sic = flux_freeze / (rho_sic*Lf) * dt

    endif

    flx_ocn = swnet + lwdown - lw - sh - lh - snow*Lf + flux_freeze
    flx_sh  = sh
    flx_lwu = lw
    flx_lh  = lh
    evp     = lh/Le

    ! frewshwater flux to the ocean, except runoff which is added in the coupler
    fw_ocn = rain + snow - evp - dh_sic*rho_sic/dt

    if (check_energy) then
      ! energy conservation check
      energy_cons = swnet &
      + lwdown        &
      - lw &
      - sh &
      - lh &
      - snow*Lf &
      + flux_freeze &
      - flx_ocn
      ! if( abs(energy_cons) .gt. 1.d-9 ) then
      !   print *,''
      !   print *,'energy balance ocn',energy_cons
      !   print *,'tskin',t_skin
      !   print *,'t_ocn,t_air',t_ocn,t_air
      !   print *,'sw,lw_d,lw_u,sh,lh,snowmelt',swnet,emiss*lwdown,lw,sh,lh,snow*Lf
      !   print *,'flux_freeze',flux_freeze
      !   print *,'dh_sic',dh_sic
      !   print *,''
      ! endif
    endif

    ! if( t_skin.gt.T0+50. .or. t_skin.lt.200._wp ) then
    !   print *,''
    !   print *,'i,j',i,j
    !   print *,'t_skin over ocn out of range!!!',t_skin
    !   print *,'sw,lw_d,lw_u,sh,lh,g',swnet,emiss*lwdown,lw,sh,lh,flx_ocn
    !   print *,'t_ocn',t_ocn
    !   print *,'t_air',t_air
    !   print *,'qair',qair
    !   print *,'Cde,Cdh',Cde,Cdh
    !   print *,''
    ! endif

    return

  end subroutine ebal_ocn

end module ebal_ocn_mod

