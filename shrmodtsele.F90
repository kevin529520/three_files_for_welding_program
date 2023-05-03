!> \brief 
!>    Evaluate element shear modulus for TSE solidified fluid region
!> \note
!>    Used ONLY for Thermal stress evoluation (TSE) model.
!>    Called both from the preprocessor to set up temp-dependent tables, and from the solver to compute local value of the shear modulus.
!
!> \return Local shear modulus (Force/area)
!> \param[in] temp  Local temperature
!
      function shrmodtsele(temp)
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
! shrmodts    - Constant user-specified shear modulus of solidified fluid (Force/area)
! ishrmodts   - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate solidified fluid shear modulus
!
      if(ishrmodts.eq.1) then
          shrmodtsele=alookup(temp,shrmodts_var1,shrmodts_var2, & 
     &shrmodts_delta,shrmodts_deltar,shrmodts_t(1),ntabmx)
      else
          shrmodtsele=shrmodts
      endif
!
      return
      end
