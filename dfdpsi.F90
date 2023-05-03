!> \brief 
!>    This function calculates the derivative of charge density with respect to zeta-potential. 
!
!>    The Boltzman distribution is used for charge density and is calculated in function fpsi.F with single valence assumed.
!
      function dfdpsi(sclrcell)
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
      include '../comdeck/scala.f90'
      include '../comdeck/state.f90'
!
! Users are allowed to provide other different distributions
! in function fpsi.F. If a different distribution is provided
! in function fpsi.F, the corresponding derivative of charge
! density with respect to zeta-potential must be
! provided here. The descriptions of related variables in
! this function are:
!
! dfdpsi---the derivative of charge density with respect
!          to zeta potential
!
! sclrcell---zeta potential
!
! other variables appearing here are input parameters and
! can be found in user's manual
!
      dfdpsi=-concen*tempsol*(exp(concen*sclrcell) & 
     &       +exp(-concen*sclrcell))
!
! start user customization
!
!     dfspsi=.........
!
! end of user customization
!
      return
      end function dfdpsi
