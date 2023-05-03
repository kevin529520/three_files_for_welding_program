!> \brief 
!>    Evaluate element poisson ratio for TSE solidified fluid region
!
!>    Used ONLY for Thermal stress evoluation (TSE) model. This routine is called
!>    only from the preprocessor to set up temp-dependent tables of shear and bulk moduli.
!
!> \return  Local Poisson ratio (D'less)
!> \param[in]   temp    Local temperature
!
      function poissonrtsele(temp)
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
! poissonrts    - Constant user-specified Poisson ratio of solidified fluid (D'less)
! ipoissonrts   - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate solidified fluid poisson ratio
!
      if(ipoissonrts.eq.1) then
          poissonrtsele=alookup(temp,poissonrts_var1,poissonrts_var2, & 
     &poissonrts_delta,poissonrts_deltar,poissonrts_t(1),ntabmx)
      else
          poissonrtsele=poissonrts
      endif
!
      return
      end
