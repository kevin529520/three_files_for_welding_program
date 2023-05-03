!> \brief 
!>    Evaluate element shear modulus for FSI component nn
!> \note
!>    Used ONLY for Fluid Structure interaction (FSI) model.
!
!>    Called both from the preprocessor to set up temp-dependent tables, and from the solver to compute local value of the shear modulus.
!
!> \return  Local shear modulus (Force/area)
!> \param[in]   nn        Current component number (1<=nn<=nobs)
!> \param[in] temp        Local temperature
!
      function shrmodele(nn,temp)
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
! shrmod(nn)  - Constant user-specified shear modulus of component nn (Force/area)
! ishrmod(nn) - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate FSI component shear modulus
!
      if(ishrmod(nn).eq.1) then
          shrmodele=alookup(temp,shrmod_var1(nn),shrmod_var2(nn), & 
     &shrmod_delta(nn),shrmod_deltar(nn),shrmod_t(1,nn),ntabmxo)
      else
          shrmodele=shrmod(nn)
      endif
!
      return
      end
