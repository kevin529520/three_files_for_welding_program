!> \brief 
!>    Evaluate heat transfer coefficient between component nn and fluid #1
!
!> \param[in]    nn   obstacle number
!> \param[in]    m    
!> \param[in]    ijk  cell ijk index
!> \param[in]    twall  
!> \param[in]    tfluid 
!
      function hobs1cl(nn,m,ijk,twall,tfluid)
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
      use fsidata_module
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
!
!    hobs1(nn)   = input htc for obstacle nn for liquid metal (fluid)
!    hobs1s(nn)  = input htc for obstacle nn for solid metal
!
      if(ihobs1(nn).eq.0) then
!
          hobs1cl=hobs1(nn)
!
! include dependence on solid fraction
!
          if(ifsld.gt.0) & 
     &      hobs1cl=fliq(ijk)*hobs1(nn)+(one-fliq(ijk))*hobs1s(nn)
!
      else
! tabular definition as a function of temperature tfluid
            hobs1cl=alookup(tfluid,hobs1_var1(nn), & 
     &              hobs1_var2(nn),hobs1_delta(nn), & 
     &              hobs1_deltar(nn),hobs1_t(1,nn),ntabmxo)
      endif
!
      return
      end
