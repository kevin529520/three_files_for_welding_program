!> \brief 
!>   Evaluate element bulk modulus for TSE solidified fluid region
!
!>   Used ONLY for Thermal stress evoluation (TSE) model. This routine is called
!>   both from the preprocessor to set up temp-dependent tables, and from the solver
!>   to compute local value of the bulk modulus.
!
!> \param[in] temp Local temperature
!
      function blkmodtsele(temp)
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
! blkmodts    - Constant user-specified bulk modulus of solidified fluid (Force/area)
! blkmodtsele - Output local bulk modulus (Force/area)
! iblkmodts   - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate solidified fluid bulk modulus
!
      if(iblkmodts.eq.1) then
          blkmodtsele=alookup(temp,blkmodts_var1,blkmodts_var2, & 
     &blkmodts_delta,blkmodts_deltar,blkmodts_t(1),ntabmx)
      else
          blkmodtsele=blkmodts
      endif
!
      return
      end
