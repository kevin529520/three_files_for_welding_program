!> \brief 
!>    Evaluate element elastic modulus for TSE solidified fluid region
!
!> \note
!>    Used ONLY for Thermal stress evoluation (TSE) model. 
!
!>    This routine is called only from the preprocessor to set up temp-dependent tables of shear and bulk moduli.
!
!> \return Local Young's elastic modulus (Force/area)
!
!> \param[in] temp  - Local temperature
!
      function elastmodtsele(temp)
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
! elastmodts    - Constant user-specified Young's elastic modulus of solidified fluid (Force/area)
! ielastmodts   - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate solidified fluid elastic modulus
!
      if(ielastmodts.eq.1) then
          elastmodtsele=alookup(temp,elastmodts_var1,elastmodts_var2, & 
     &elastmodts_delta,elastmodts_deltar,elastmodts_t(1),ntabmx)
      else
          elastmodtsele=elastmodts
      endif
!
      return
      end
