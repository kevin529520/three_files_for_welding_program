!> \brief 
!>    Evaluate forces on fluid due to electric field effects.
!
!>    Includes forces on charged fluid and dielectrophoresis forces when idep = 1
!
!>    Includes electrosmotic forces when ieof2 > 0
!
!> \param[in]   idir   
!> \verbatim
!>    1: for x-direction
!>    2: for y-direction
!>    3: for z-direction
!> \endverbatim
!
!> \param[in]   i     Cell index in x-direction
!> \param[in]   j     Cell index in y-direction
!> \param[in]   k     Cell index in z-direction
!> \param[out]  efor  Force value
!> \param[in]   term1     
!> \param[in]   term2    
!
      subroutine eforcal(idir,i,j,k,efor,term1,term2)
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
      use avgco_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/const.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/edit.f90'
      include '../comdeck/scala.f90'
      include '../comdeck/state.f90'
!
      efor=zero
      include '../comdeck/ijk.f90'
!
      if(idir.eq.1) then
        ipjk_str=ijk_str+1
        ipjk=ijk_str2unstr(ipjk_str)
!
! include dielectrophoretic and free charge forces
!
        if(iechrg.gt.0) then
          fvchg=zero
          fvchgp=zero
!
          if(nmat.eq.1) then
            if(isclr(iechrg).eq.3 .and. econd1.gt.ztest) then
              fvchg=f(ijk)
              fvchgp=f(ipjk)
            endif
          else
            if((isclr(iechrg).eq.1 .or. isclr(iechrg).eq.3) .and. & 
     &           econd1.gt.ztest) then
                fvchg=f(ijk)
                fvchgp=f(ipjk)
            endif
            if((isclr(iechrg).eq.1 .or. isclr(iechrg).eq.5) .and. & 
     &           econd2.gt.ztest) then
                fvchg=fvchg+one-f(ijk)
                fvchgp=fvchgp+one-f(ipjk)
            endif
          endif
!
          efor=(rdxp(i)*sclr(ijk,iechrg)*fvchg*vf(ijk) & 
     &                   +rdxm(i)*sclr(ipjk,iechrg)*fvchgp*vf(ipjk)) & 
     &                   *(epot(ijk_str)-epot(ipjk_str))
        endif
!
        if(idep.gt.0) then
          dielr=dielc(ijk_str)+dielc(ipjk_str)
          if(dielr.gt.ztest) then
            dielr=two*dielc(ijk_str)*dielc(ipjk_str)/dielr
          else
            dielr=one
          endif
          efor=efor+half*elperm*(dielr-one)* & 
     &         (efield(ipjk_str)-efield(ijk_str))
        endif
!
! include electrosmotic forces
!
        if(ieof2.gt.0) then
          efort=rdxp(i)*sclr(ijk,ieof2)+ & 
     &         rdxm(i)*sclr(ipjk,ieof2)
          efort=-fpsi(efort)
          efor=efor+efort*(epot(ipjk_str)-epot(ijk_str))
        endif
!
! include electrothermal forces
!
        if(iethermo.eq.2) then
          efor=efor+term2*(rdxp(i)*efield(ijk_str)+ & 
     &                     rdxm(i)*efield(ipjk_str))*(tn(ijk)-tn(ipjk)) & 
     &             +term1*(epot(ijk_str)-epot(ipjk_str))
        endif
      elseif(idir.eq.2) then
        ijpk_str=ijk_str+ii1
        ijpk=ijk_str2unstr(ijpk_str)
!
! include dielectrophoretic and free charge forces
!
        if(iechrg.gt.0) then
          fvchg=zero
          fvchgp=zero
!
          if(nmat.eq.1) then
            if(isclr(iechrg).eq.3 .and. econd1.gt.ztest) then
              fvchg=f(ijk)
              fvchgp=f(ijpk)
            endif
          else
            if((isclr(iechrg).eq.1 .or. isclr(iechrg).eq.3) .and. & 
     &           econd1.gt.ztest) then
                fvchg=f(ijk)
                fvchgp=f(ijpk)
            endif
            if((isclr(iechrg).eq.1 .or. isclr(iechrg).eq.5) .and. & 
     &           econd2.gt.ztest) then
                fvchg=fvchg+one-f(ijk)
                fvchgp=fvchgp+one-f(ijpk)
            endif
          endif
!
          efor=(rdyp(j)*sclr(ijk,iechrg)*fvchg*vf(ijk) & 
     &                   +rdym(j)*sclr(ijpk,iechrg)*fvchgp*vf(ijpk)) & 
     &                   *(epot(ijk_str)-epot(ijpk_str))
        endif
!
        if(idep.gt.0) then
          dielbk=dielc(ijk_str)+dielc(ijpk_str)
          if(dielbk.gt.ztest) then
            dielbk=two*dielc(ijk_str)*dielc(ijpk_str)/dielbk
          else
            dielbk=one
          endif
          efor=efor+half*elperm*(dielbk-one)* & 
     &         (efield(ijpk_str)-efield(ijk_str))
        endif
!
! include electrosmotic forces
!
        if(ieof2.gt.0) then
          efort=rdyp(j)*sclr(ijk,ieof2)+ & 
     &          rdym(j)*sclr(ijpk,ieof2)
          efort=-fpsi(efort)
          efor=efor+efort*(epot(ijpk_str)-epot(ijk_str))
        endif
!
! include electrothermal forces
!
        if(iethermo.eq.2) then
          efor=efor+term2*(rdyp(j)*efield(ijk_str)+ & 
     &                     rdym(j)*efield(ijpk_str))*(tn(ijk)-tn(ijpk)) & 
     &             +term1*(epot(ijk_str)-epot(ijpk_str))
        endif
      elseif(idir.eq.3) then
        ijkp_str=ijk_str+ii2
        ijkp=ijk_str2unstr(ijkp_str)
!
! include dielectrophoretic and free charge forces
!
        if(iechrg.gt.0) then
          fvchg=zero
          fvchgp=zero
!
          if(nmat.eq.1) then
            if(isclr(iechrg).eq.3 .and. econd1.gt.ztest) then
              fvchg=f(ijk)
              fvchgp=f(ijkp)
            endif
          else
            if((isclr(iechrg).eq.1 .or. isclr(iechrg).eq.3) .and. & 
     &           econd1.gt.ztest) then
                fvchg=f(ijk)
                fvchgp=f(ijkp)
            endif
            if((isclr(iechrg).eq.1 .or. isclr(iechrg).eq.5) .and. & 
     &           econd2.gt.ztest) then
                fvchg=fvchg+one-f(ijk)
                fvchgp=fvchgp+one-f(ijkp)
            endif
          endif
!
          efor=(rdzp(k)*sclr(ijk,iechrg)*fvchg*vf(ijk) & 
     &                   +rdzm(k)*sclr(ijkp,iechrg)*fvchgp*vf(ijkp)) & 
     &                   *(epot(ijk_str)-epot(ijkp_str))
        endif
!
        if(idep.gt.0) then
          dielt=dielc(ijk_str)+dielc(ijkp_str)
          if(dielt.gt.ztest) then
            dielt=two*dielc(ijk_str)*dielc(ijkp_str)/dielt
          else
            dielt=one
          endif
          efor=efor+half*elperm*(dielt-one)* & 
     &         (efield(ijkp_str)-efield(ijk_str))
        endif
!
! include electrosmotic forces
!
        if(ieof2.gt.0)then
          efort=rdzp(k)*sclr(ijk,ieof2)+ & 
     &          rdzm(k)*sclr(ijkp,ieof2)
          efort=-fpsi(efort)
          efor=efor+efort*(epot(ijkp_str)-epot(ijk_str))
        endif
!
! include electrothermal forces
!
        if(iethermo.eq.2) then
          efor=efor+term2*(rdzp(k)*efield(ijk_str)+ & 
     &                     rdzm(k)*efield(ijkp_str))*(tn(ijk)-tn(ijkp)) & 
     &             +term1*(epot(ijk_str)-epot(ijkp_str))
        endif
      endif
!
      return
!
      end
