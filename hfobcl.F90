!> \brief 
!>    Compute heat transfer coefficient between metal and foam
!
!> \param[in]      nn         Component number 
!> \param[in]       m         Index of the component surface element 
!> \param[in]     ijk         Cell ijk index
!> \param[out]  htran         Computed heat transfer coefficient
!
      subroutine hfobcl(nn,m,ijk,htran)
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
      use topo_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/dparam.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/diag.f90'
      include '../comdeck/logvar.f90'
      include '../comdeck/obsd.f90'
      include '../comdeck/obsdef.f90'
      include '../comdeck/bcinfo.f90'
      include '../comdeck/acctrm.f90'
      include '../comdeck/scala.f90'
      include '../comdeck/state.f90'
!
!     gxg,gyg,gzg = components of gravitational acceleration
!     cgfob=coefficient of proportionality in "roughness" model
!     rough(nn)=characteristic length for heat transfer "roughness"
!
      htran=hobs1cl(nn,m,ijk,zero,tn(ijk))
!
! ------- include pressure and temperature dependencies
!     p(ijk) = pressure of ijk cell
!     tn(ijk) = fluid temperature of ijk cell
!     xi(i),yj(j),zk(k) = center of ijk cell
!     i,j,k recovered from ijk_str by:
!           i=i_str(ijk_str), j=j_str(ijk_str), k=k_str(ijk_str)
!     idum1,...,idum25 = dummy integer parameters for model building
!     dum1,...,dum25 = dummy floating point parameters for model building
! -------
!
!  add gravitational "roughness" effect (foaml is characteristic length)
!
      if(cgfob.gt.ztest .and. rough(nn).gt.ztest) then
!
!  compute outward unit normal to metal/foam interface
!
        call srfnorm(f,ijk,snormx,snormy,snormz,sar,leneax,1)
        if(isurf.eq.1) then
          call srfnew(f,ijk,snormx,snormy,snormz,sar,leneax,1,.true.)
        endif
        ght=gxg*snormx+gyg*snormy+gzg*snormz
!
        if(abs(ght).gt.ztest) then
          gvel=sqrt(abs(ght)*rough(nn))
          htcor=sign(cgfob,ght)*gvel/(gvel+htran/rcobs(nn))
          htran=htran*(one+htcor)
        endif
!
      endif
!
      return
      end
