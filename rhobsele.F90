!> \brief 
!>    Evaluate element density for FSI component nn
!
!> \note
!>    Used ONLY for Fluid Structure interaction (FSI) model. 
!
!>    This routine is called both from the preprocessor to set up temp-dependent tables,
!>    and from the solverto compute local value of the density.
!
!> \param[in] nn    Current component number (1<=nn<=nobs)
!> \param[in] temp  Local temperature
!
      function rhobsele(nn,temp)
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
! rhobs(nn)  - Constant user-specified bulk modulus of component nn (mass/volume)
! rhobsele   - Output local bulk modulus (mass/volume)
! irhobs(nn) - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate FSI component density
!
      if(irhobs(nn).eq.1) then
          rhobsele=alookup(temp,rhobs_var1(nn),rhobs_var2(nn), & 
     &rhobs_delta(nn),rhobs_deltar(nn),rhobs_t(1,nn),ntabmxo)
      else
          rhobsele=rhobs(nn)
      endif
!
      return
      end
