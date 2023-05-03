      function scour_uset(den,dia,visc,rhof,gtot,totvfsp)
!
!    ***************************************************************
!    **                          notice                           **
!    **  this subprogram contains flow science, inc. proprietary  **
!    **       trade secret and confidential information.          **
!    **                                                           **
!    **                unauthorized use prohibited                **
!    **           copyright 1985-2016 flow science, inc.          **
!    ***************************************************************
!
!     Routine for modeling scour settling velocity in shallow water
!     This routine can be user customized.
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/const.f90'
      include '../comdeck/scala.f90'
!
! ------ The following list describes parameters needed for customization
!
!      variable          description
!      --------        --------------------------------
!       scour_uset       settling velocity
!       gtot             Magnitude of gravity vector
!       den              Density of sediment species,
!                        is input parameter rhfss(k), k is species index
!       dia              Grain diameter
!                        is input parameter ssdia(k), k is species index
!       visc             Dynamic viscosity of fluid
!       rhof             FLuid density
!       totvfsp          Total volume fraction of suspended sediments
!       ssrzmlt          User-defined multiplier for Richarson-Zaki coefficient, default is 1.0
!
! --------------------------------------------------------------------------
!
!   Calculate dimensionless mean particle diameter
!
      dstar=dia*(rhof*max(zero,den-rhof)*gtot/visc**2)**ob3
!
!   Calculate settling velocity (in length/time),
!   From Soulsby (1997), Dynamics of Marine Sands, p134
!
      scour_uset=(sqrt(10.36d0**2+1.049d0*dstar**3)-10.36d0)*visc/ & 
     &                (rhof*dia)
!
!   Add Richardson-Zaki correlation
!   (Ref: Trans.Instn.Chem.Engrs.,Vol.32, p.35, 1954)
!   (See also Felice & Kehlenbeck, Chem.Eng.Technol.23,(2000)12)
!
      repart=rhof*dia*scour_uset/visc
      if(repart.le.0.2d0) then
          rzexp=4.65d0
      elseif(repart.le.one) then
          rzexp=4.35d0/repart**0.03d0
      elseif(repart.le.500.d0) then
          rzexp=4.45d0/repart**0.1d0
      else
          rzexp=2.39d0
      endif
      rzexp=ssrzmlt*rzexp
      if(totvfsp.lt.scrfcr) then
        scour_uset=scour_uset*(one-min(half,totvfsp))**rzexp
      else
        scour_uset=zero   ! at the critical fraction (scrfcr),
                          ! the suspended sediments become fully packed,
                          ! thus no settling can happen
      end if
!
      return
      end
!
