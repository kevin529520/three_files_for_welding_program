!> \brief 
!>    Evaluate element elastic modulus for FSI component nn
!
!> \note
!>    Used ONLY for Fluid Structure interaction (FSI) model. 
!
!>    This routine is called from the preprocessor to set up temp-dependent tables of shear and bulk moduli.
!
!> \return   Local Young's elastic modulus (Force/area)
!
!> \param[in] nn    Current component number (1<=nn<=nobs)
!> \param[in] temp  Local temperature
!
      function elastmodele(nn,temp)
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
! elastmod(nn)  - Constant user-specified Young's elastic modulus of component nn (Force/area)
! ielastmod(nn) - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate FSI component elastic modulus
!
      if(ielastmod(nn).eq.1) then
          elastmodele=alookup(temp,elastmod_var1(nn),elastmod_var2(nn), & 
     &elastmod_delta(nn),elastmod_deltar(nn),elastmod_t(1,nn),ntabmxo)
      else
          elastmodele=elastmod(nn)
      endif
!
      return
      end
