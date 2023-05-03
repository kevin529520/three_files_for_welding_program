!> \brief 
!>    Subroutine for mass transfer coefficient between solid solute and fluid
!
!> \param[in]  nn     Obstacle number
!> \param[in]  ijk    Cell ijk index
!> \param[out] htran  Computed mass transfer coefficient
!
      subroutine smtrcl(nn,ijk,htran)
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
      use voids_module
!
      use obsijk_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
!
      include '../comdeck/params.f90'
      include '../comdeck/dparam.f90'
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
!     dissmtc=input mass transfer coffficient for obstacle nn
!     gxg,gyg,gzg = components of gravitational acceleration
!     rough(nn)=characteristic length for heat transfer "roughness"
!
! limit mass transfer by a diffusion length sqrt(D*dt)
!
      diff_coeff=cmsc(idissolve)
      if(ifvis.ne.0) diff_coeff=diff_coeff+rmscs(idissolve)*mu(ijk)
      if(diff_coeff.gt.ztest) then
        htran=min(dissmtc,sqrt(diff_coeff/rho(ijk)*rdelt))
      else
        htran=dissmtc
      endif
!
!
! ------- to include pressure and temperature dependencies:
!
!     p(ijk) = pressure of ijk cell
!     tn(ijk) = fluid temperature of ijk cell
!     xi(i),yj(j),zk(k) = center of ijk cell
!     i,j,k recovered from ijk_str:
!        i=i_str(ijk_str)
!        j=j_str(ijk_str)
!        k=k_str(ijk_str)
!     idum1,...,idum25 = dummy integer parameters for model building
!     dum1,...,dum25 = dummy floating point parameters for model building
!
! -------
!
      return
      end
