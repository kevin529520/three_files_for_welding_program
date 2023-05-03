!> \brief 
!>    Evaluate cell derivative d(rhoe)/d(temp) using 'n' time level variables, assuming density is constant
!
      function rhocvn(ijk)
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
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/state.f90'
      include '../comdeck/phiou.f90'
!
      ff=fn(ijk)
      cv1_loc=cvcal(tn(ijk),1)
!
      if(nmat.eq.1) then
        rhocvn=ff*rhof*cv1_loc
      else
        cv2_loc=cvcal(tn(ijk),2)
        if(icmprs.eq.1) then
          if(incg.gt.0) then
            cvnc_loc=cvcal(tn(ijk),3)
            rho2t=rhon(ijk)-ff*rhof
            if(rho2t.lt.ztest) rho2t=rho2cl(p(ijk),tn(ijk), & 
     &             rf2t(ijk))*(one-ff)
            rhonct=sclrn(ijk,incg)*(one-ff)
            rhocvn=rho2t*cv2_loc+max(zero,min(rho2t,rhonct)) & 
     &             *(cvnc_loc-cv2_loc)+ff*rhof*cv1_loc
          else
            rhocvn=rhon(ijk)*cv2_loc+ff*rhof*(cv1_loc-cv2_loc)
          endif
        else
          rhocvn=ff*rhof*cv1_loc+(one-ff)*rhof2*cv2_loc
        endif
!        
        if(rhocvn.lt.ztest) then
          include '../comdeck/ijk2str.f90'
          ii=i_str(ijk_str)
          jj=j_str(ijk_str)
          kk=k_str(ijk_str)
          write(iout,1000) cycle,ii,jj,kk,rho(ijk),ff, & 
     &       tn(ijk),p(ijk),rhoe(ijk)
          write(ierrfl,1000) cycle,ii,jj,kk,rho(ijk),ff, & 
     &       tn(ijk),p(ijk),rhoe(ijk)
          rhocvn=rhof*cv1*ff+em6
        endif
      endif
!
! add latent heat release
!
      if(ifsld.gt.0) then
        dfsolid_dt=dlookup(tn(ijk),fs1_var1,fs1_var2, & 
     &                      fs1_deltar,fs1_t,ntabmx)
        rhocvn=rhocvn-ff*rhof*dfsolid_dt*clht1
      endif
!
 1000 format(//,5x,'**** negative rhocvn ****',/, & 
     &  5x,'cycle=',i6,' i=',i4,' j=',i4,' k=',i4,' rho=',1p,e12.5, & 
     &     ' fn=',1p,e12.5,' tn=',1p,e12.5,' p=',1p,e12.5, & 
     &     ' rhoe=',1p,e12.5,//)
!
      return
      end
!
! #############################################################################
!
      function rhocv(ijk)
!
!    evaluate cell derivative d(rhoe)/d(temp) using 'n+1'
!    time level variables, assuming density is constant
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
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/state.f90'
      include '../comdeck/phiou.f90'
!
      ff=f(ijk)
      cv1_loc=cvcal(tn(ijk),1)
!
      if(nmat.eq.1) then
        rhocv=ff*rhof*cv1_loc
      else
        cv2_loc=cvcal(tn(ijk),2)
        if(icmprs.eq.1) then
          if(incg.gt.0) then
            cvnc_loc=cvcal(tn(ijk),3)
            rho2t=rho(ijk)-ff*rhof
            if(rho2t.lt.ztest) rho2t=rho2cl(p(ijk),tn(ijk), & 
     &             rf2t(ijk))*(one-ff)
            rhonct=sclr(ijk,incg)*(one-ff)
            rhocv=rho2t*cv2_loc+max(zero,min(rho2t,rhonct)) & 
     &             *(cvnc_loc-cv2_loc)+ff*rhof*cv1_loc
          else
            rhocv=rho(ijk)*cv2_loc+ff*rhof*(cv1_loc-cv2_loc)
          endif
        else
          rhocv=ff*rhof*cv1_loc+(one-ff)*rhof2*cv2_loc
        endif
!        
        if(rhocv.lt.ztest) then
          include '../comdeck/ijk2str.f90'
          ii=i_str(ijk_str)
          jj=j_str(ijk_str)
          kk=k_str(ijk_str)
          write(iout,1000) cycle,ii,jj,kk,rho(ijk),ff, & 
     &       tn(ijk),p(ijk),rhoe(ijk)
          write(ierrfl,1000) cycle,ii,jj,kk,rho(ijk),ff, & 
     &       tn(ijk),p(ijk),rhoe(ijk)
          rhocv=rhof*cv1*ff+em6
        endif
      endif
!
! add latent heat release
!
      if(ifsld.gt.0) then
        dfsolid_dt=dlookup(tn(ijk),fs1_var1,fs1_var2, & 
     &                      fs1_deltar,fs1_t,ntabmx)
        rhocv=rhocv-ff*rhof*dfsolid_dt*clht1
      endif
!
 1000 format(//,5x,'**** negative rhocv ****',/, & 
     &  5x,'cycle=',i6,' i=',i4,' j=',i4,' k=',i4,' rho=',1p,e12.5, & 
     &     ' fn=',1p,e12.5,' tn=',1p,e12.5,' p=',1p,e12.5, & 
     &     ' rhoe=',1p,e12.5,//)
!
      return
      end
