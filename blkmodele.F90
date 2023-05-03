!> \brief 
!>    Evaluate element bulk modulus for FSI component nn
!
!>   Used ONLY for Fluid Structure interaction (FSI) model. This routine is called
!>   both from the preprocessor to set up temp-dependent tables, and from the solver
!>   to compute local value of the bulk modulus.
!
!> \param[in] nn    Current component number (1<=nn<=nobs)
!> \param[in] temp  Local temperature
!
      function blkmodele(nn,temp)
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
! blkmod(nn)  - Constant user-specified bulk modulus of component nn (Force/area)
! blkmodele   - Output local bulk modulus (Force/area)
! iblkmod(nn) - Flag: =0, use constant; =1, use user-defined temperature-dependent table
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
!     evaluate FSI component bulk modulus
!
      if(iblkmod(nn).eq.1) then
          blkmodele=alookup(temp,blkmod_var1(nn),blkmod_var2(nn), & 
     &blkmod_delta(nn),blkmod_deltar(nn),blkmod_t(1,nn),ntabmxo)
      else
          blkmodele=blkmod(nn)
      endif
!
      return
      end
