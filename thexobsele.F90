!> \brief 
!>    Evaluate element thermal expansion coeff for FSI component nn
!
!> \note
!>    Used ONLY for Fluid Structure interaction (FSI) model
!
!> \return Thermal expansion value needed in model (density/Temperature)
!> \param[in]    nn     Current component number (1<=nn<=nobs)
!> \param[in]    temp   Local temperature
!
      function thexobsele(nn,temp)
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
! thexobs(nn) - Thermal expansion (1/Temperature)
! irhobs(nn)  - Flag: =0, use constant; =2, use user-defined temperature-dependent table
!
      use arrays_module
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
      include '../comdeck/state.f90'
!
!            evaluate cell fluid thermal conductivity
!
      if(irhobs(nn).eq.1) then
          thexobsele=dlookup(temp,rhobs_var1(nn),rhobs_var2(nn), & 
     &                rhobs_deltar(nn),rhobs_t(1,nn),ntabmxo)
      else
          thexobsele=-thexobs(nn)*rhobs(nn)
      endif
!
      return
      end
