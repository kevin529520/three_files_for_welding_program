!> \brief 
!>    Evaluate element yield stress limit for TSE solidified fluid region
!
!>    Used ONLY for Fluid Structure interaction (FSI) model. This routine is 
!>    Called both from the preprocessor to set up temp-dependent tables, and from the solver to compute local value of the yield stress limit.
!
!> \return Local yield stress limit (force/area)
!> \param[in] temp    Local temperature
!
      function yldstrtsele(temp)
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
! yldstrts    - Constant user-specified yield stress limit of solidified fluid (Force/area)
! iyldstrts   - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate solidified fluid yield stress limit
!
      if(iyldstrts.eq.1) then
          yldstrtsele=alookup(temp,yldstrts_var1,yldstrts_var2, & 
     &yldstrts_delta,yldstrts_deltar,yldstrts_t(1),ntabmx)
      else
          yldstrtsele=yldstrts
      endif
!
      return
      end
