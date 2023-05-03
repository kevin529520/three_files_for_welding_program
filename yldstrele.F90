!> \brief 
!>    Evaluate element yield stress limit for FSI component nn
!> \note
!>    Used ONLY for Fluid Structure interaction (FSI) model. 
!>    Calledboth from the preprocessor to set up temp-dependent tables, and from the solver to compute local value of the yield stress limit.
!
!> \return Local yield stress limit (force/area)
!> \param[in] nn      Current component number (1<=nn<=nobs)
!> \param[in] temp    Local temperature
!
      function yldstrele(nn,temp)
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
! yldstr(nn)  - Constant user-specified yield stress limit of component nn (force/area)
! iyldstr(nn) - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate yield stress limit for FSI component
!
      if(iyldstr(nn).eq.1) then
          yldstrele=alookup(temp,yldstr_var1(nn),yldstr_var2(nn), & 
     &     yldstr_delta(nn),yldstr_deltar(nn),yldstr_t(1,nn),ntabmxo)
      else
          yldstrele=yldstr(nn)
      endif
!
      return
      end
