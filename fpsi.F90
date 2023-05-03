      function fpsi(sclrcell)
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
! This function calculates the charge density distribution
! used in electro-osmosis model. Here Boltzmann distribution
! with single valence is used to solve for zeta-potential.
! Users are allowed to provide other different distribution
! here.
!
! If a different distribution is provided here, please
! modify function dfdpsi.F accordingly, where the
! derivative of charge density with respect to
! zeta potential is calculated.
!
! The descriptions on related variables in this function
! are:
!
! fpsi---charge density
!
! sclrcell---zeta-potential
!
! other variabless appearing here are input parameters
! and their definitions can be found in user's manual
!
      fpsi=-tempsol*(exp(concen*sclrcell)-exp(-concen*sclrcell))
!
! start user customization
!
!     fpsi=..................
!
! end of user customization
!
      return
      end function fpsi
