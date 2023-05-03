!> \brief 
!>    Evaluate heat transfer coefficient between component nn and fluid #2
!
!> \param[in]    nn     obstacle number
!> \param[in]    ijk    cell ijk index
!> \param[in]    twall    
!> \param[in]    tfluid     
!
      function hobs2cl(nn,ijk,twall,tfluid)
!
!  evaluate heat transfer coefficient between component #nn and fluid #2
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
      use arrays_module
!
      use meshcb_module
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
      include '../comdeck/diag.f90'
      include '../comdeck/logvar.f90'
      include '../comdeck/obsd.f90'
      include '../comdeck/bcinfo.f90'
      include '../comdeck/acctrm.f90'
      include '../comdeck/scala.f90'
      include '../comdeck/state.f90'
      include '../comdeck/dumn.f90'
!
!    hobs2(nn)   = input htc for obstacle nn
!
      if(ihobs2(nn).eq.0) then
!
          hobs2cl=hobs2(nn)
!
      else
! tabular definition as a function of temperature tfluid
          hobs2cl=alookup(tfluid,hobs2_var1(nn), & 
     &            hobs2_var2(nn),hobs2_delta(nn), & 
     &            hobs2_deltar(nn),hobs2_t(1,nn),ntabmxo)
      endif
!
      return
      end
