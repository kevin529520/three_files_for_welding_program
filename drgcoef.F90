!> \brief 
!
!> \param[in]     df  depth of fluid
!> \param[in]     rgh roughness
!
      function drgcoef(df,rgh)
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
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/const.f90'
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     This subroutine calculates the drag coefficient for quadratic law of bottom shear stress
!     in shallow water model.
!
!     The subroutine is open for customization by users
!
!     The equation is from "Dynamics of Marine Sands" by Soulsby, 1997, page 55, equation(36).
!     The upper limit value of cd is cd_max.
!
!     bdrg: =0.71, defined in const.F. It is B in equation (36)
!     ckar: =0.4, defined as const.F. It is k in equation (36)
!     cd_max: =0.01, upper limit of drag coefficient defined in gsetup.F.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      z0=rgh/30.d0
      z0df=z0/df
      z0df=min(z0df,1.d-2)
      z0df=max(z0df,1.d-7)
      drgcoef=(ckar/(bdrg+log(z0df)))**2
!
      drgcoef=min(drgcoef,cd_max)
!
      return
      end function
!
