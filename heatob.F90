!> \brief 
!>    Add heat sources to obstacles
!> \return Thermal energy per unit volume of obstacle material per unit time (or power per unit volume)
!> \param[in]   i         Current cell index in x-direction
!> \param[in]   j         Current cell index in y-direction
!> \param[in]   k         Current cell index in z-direction
!> \param[in]  nob        Obstacle number in cell (i,j,k)
!
      function heatob(i,j,k,nob)
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
      use arrays_module
!
      use meshcb_module
!
      use obsijk_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/diag.f90'
      include '../comdeck/logvar.f90'
      include '../comdeck/obsd.f90'
      include '../comdeck/bcinfo.f90'
      include '../comdeck/acctrm.f90'
      include '../comdeck/scala.f90'
      include '../comdeck/state.f90'
      include '../comdeck/dumn.f90'
!
      heatob=zero
!
!cc     include '../comdeck/ijk.f90'
!cc     include '../comdeck/ijksd.f90'
!
      return
      end
