!> \brief 
!>    Evaluate element poisson ratio for FSI component nn
!
!>    Used ONLY for Fluid Structure interaction (FSI) model. This routine is called
!>    only from the preprocessor to set up temp-dependent tables of shear and bulk moduli.
!
!> \return Local Poisson ratio (D'less)
!> \param[in]   nn    Current component number (1<=nn<=nobs)
!> \param[in]   temp  Local temperature
!
      function poissonrele(nn,temp)
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
! poissonr(nn)  - Constant user-specified Poisson ratio of component nn (D'less)
! ipoissonr(nn) - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate FSI component poisson ratio
!
      if(ipoissonr(nn).eq.1) then
          poissonrele=alookup(temp,poissonr_var1(nn),poissonr_var2(nn), & 
     &poissonr_delta(nn),poissonr_deltar(nn),poissonr_t(1,nn),ntabmxo)
      else
          poissonrele=poissonr(nn)
      endif
!
      return
      end
