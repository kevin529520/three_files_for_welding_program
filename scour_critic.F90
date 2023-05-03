!> \brief 
!>    Calculate the critical Shields paramter
!
!> \param[in]      crtshields     Critical Shields number
!> \param[in]      srfang         Locally computed slope angle
!> \param[in]      phiang         Locally computed angle between upslope direction and tangential flow direction
!> \param[in]      repang         Angle of repose
!> \param[in]      dia            Diameter of sediment species
!> \param[in]      den            Density of sediment species
!> \param[in]      d50            Median grain diameter over all sediment species
!> \param[in]      visc           Dynamic viscosity of fluid
!> \param[in]      rhof           Fluid density
!> \param[in]      gtot           Magnitude of gravity vector
!
      subroutine scour_critic(crtshields,srfang,phiang,repang,dia,den, & 
     &                 d50,visc,rhof,gtot)
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
      include '../comdeck/cntrl.f90'
      include '../comdeck/scala.f90'
!
! ------ The following list describes parameters needed for customization
!
!      variable        description
!      --------        --------------------------------
!       nsspc          Total number of sediment species
!       isscrt         =1, calculate critical shields number using Shields-Rouse equation
!                      =0 (default), do not calculate shields number. use prescribed or default value
!       issslp         =1, adjust critical Shields number for slope
!                      =0, do not adjust shields number for slope
!
! --------------------------------------------------------------------------
!
      dstar=dia*(rhof*max(zero,den-rhof)*gtot/visc**2)**ob3
!
!     Use Soulsby-Whitehouse equation (1997)
!     Presented in "Dynamics of Marine Sands" by Soulsby, 1997, page 106, equation(77).
      if(isscrt.eq.1) crtshields=0.3d0/(1.d0+1.2d0*dstar)+ & 
     &             0.055d0*(1.d0-exp(-0.02d0*dstar))
!
!     Adjust for slope. From Soulsby, R., Dynamics of Marine Sands, Thomas Telford Publications,
!     London, 1997, p107.
!
      if(issslp.eq.1) then
        if(repang.gt.ztest .and. srfang.gt.ztest) then
          critadj=(cos(phiang)*sin(srfang) & 
     &         + sqrt(max(zero,cos(srfang)**2*tan(repang)**2 & 
     &                        -sin(phiang)**2*sin(srfang)**2)) & 
     &          )/tan(repang)
          crtshields=crtshields*critadj
        end if
      endif
!
      return
      end
!
