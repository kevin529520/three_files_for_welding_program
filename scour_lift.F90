      function scour_lift(shields,crtshields,alfa,den,dia,visc, & 
     &         rhof,gtot)
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
!     Routine for modeling scour -- This routine can be user customized.
!     This routine computes the entrainment rate in terms of a volumetric
!     flux (i.e. velocity) of the packed sediment away from the packed
!     bed interface.  The model currently implemented here is that by
!     Mastbergen, D.R. and Von den Berg, J.H. "Breaching in fine sands
!     and the generation of sustained turbidity currents in submarine
!     canyons", Sedimentology (50) 625-637, 2003.
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/const.f90'
!
! ------ The following list describes parameters needed for customization
!
!      variable        description
!      --------        --------------------------------
!       scour_lift     Lift velocity
!       gtot           Magnitude of gravitational acceleration
!       den            Density of sediment species,
!                      is input parameter rhfss(k), k is species index
!       shields        Shields number
!       crtshields     Critical Shields number
!       dia            Grain diameter
!                      is input parameter ssdia(k), k is species index
!       visc           Dynamic viscosity of fluid
!       rhof           FLuid density
!       alfa           Sediment entrainment parameter for sediment species,
!                      is input parameter ssentrf(k),k is species index
! --------------------------------------------------------------------------
!     Initialize lift velocity scour_lift
!
      scour_lift=zero
!
      if(shields.gt.crtshields) then
!
!       Calculate dimensionless mean particle diameter
        dstar=dia*(rhof*max(zero,den-rhof)*gtot/visc**2)**ob3
!
        xmparam=1.5d0
        xnparam=0.3d0
        reldens=max(zero,den-rhof)/rhof
!
!       calculate lift velocity (in length/time)
        scour_lift=alfa*(shields-crtshields)**xmparam & 
     &      *dstar**xnparam*sqrt(reldens*gtot*dia)
!
      endif
!
      return
      end function scour_lift
!
