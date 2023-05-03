!> \brief 
!>    Evaluate turbulent energy and dissipation from the logarithmic form of the law of the wall
!
      subroutine tkewbc(i,j,k)
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
      use mblock_module
!
      use topo_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/bcinfo.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/obsd.f90'
      include '../comdeck/obsdef.f90'
      include '../comdeck/phiou.f90'
      include '../comdeck/state.f90'
      include '../comdeck/trbcom.f90'
      include '../comdeck/scala.f90'
!
      real(kind(zzzz)) :: mucal
!
!               evaluate local normal to the wall
!
      m=0
      ucs=zero
      vcs=zero
      wcs=zero
      ucsi=zero
      vcsi=zero
      wcsi=zero
!
      include '../comdeck/ijk.f90'
      include '../comdeck/mijk.f90'
!
      bal=one-afr(imjk)
      bar=one-afr(ijk)
      if(i.eq.2 .and. wl.eq.1) bal=bar
      if(i.eq.im1 .and. wr.eq.1) bar=bal
      dbadx=rdx(i)*(bar-bal)
!
      baf=one-afb(ijmk)
      babk=one-afb(ijk)
      if(j.eq.2 .and. wf.eq.1) baf=babk
      if(j.eq.jm1 .and. wbk.eq.1) babk=baf
      dbady=rdy(j)*rri(i)*(babk-baf)
!
      bab=one-aft(ijkm)
      bat=one-aft(ijk)
      if(k.eq.2 .and. wb.eq.1) bab=bat
      if(k.eq.km1 .and. wt.eq.1) bat=bab
      dbadz=rdz(k)*(bat-bab)
!
!              reciprocal normalization number
!
      banrms=dbadx*dbadx+dbady*dbady+dbadz*dbadz
      if(banrms.lt.ztest) go to 1000
      banrms=one/banrms
!
!              calculate tangential velocity components
!              from   utan=u-(n dot u)n
!
      uc=half*(u(imjk)+u(ijk))
      vc=half*(v(ijmk)+v(ijk))
      wc=half*(w(ijkm)+w(ijk))
! ------ evaluate fluid velocities at sliding walls
      if(i.eq.2 .and. wl.eq.2) then
         vcsi=vcsi+one
         wcsi=wcsi+one
         vcs=vcs+vbc(1,nbl)/rr(1)
         wcs=wcs+wbc(1,nbl)
      endif
      if(i.eq.im1 .and. wr.eq.2) then
         vcsi=vcsi+one
         wcsi=wcsi+one
         vcs=vcs+vbc(2,nbl)
         wcs=wcs+wbc(2,nbl)
      endif
      if(j.eq.2 .and. wf.eq.2) then
         ucsi=ucsi+one
         wcsi=wcsi+one
         ucs=ucs+ubc(3,nbl)
         wcs=wcs+wbc(3,nbl)
      endif
      if(j.eq.jm1 .and. wbk.eq.2) then
         ucsi=ucsi+one
         wcsi=wcsi+one
         ucs=ucs+ubc(4,nbl)
         wcs=wcs+wbc(4,nbl)
      endif
      if(k.eq.2 .and. wb.eq.2) then
         ucsi=ucsi+one
         vcsi=vcsi+one
         ucs=ucs+ubc(5,nbl)
         vcs=vcs+vbc(5,nbl)/rri(i)
      endif
      if(k.eq.km1 .and. wt.eq.2) then
         ucsi=ucsi+one
         vcsi=vcsi+one
         ucs=ucs+ubc(6,nbl)
         vcs=vcs+vbc(6,nbl)/rri(i)
      endif
!
      ishf=1
      if(abs(dbady).gt.abs(dbadx)) ishf=ii1
      if(abs(dbadz).gt.abs(dbadx) .and. & 
     &   abs(dbadz).gt.abs(dbady)) ishf=ii2
        if(ishf.eq.1 .and. dbadx.lt.zero) ishf=-1
        if(ishf.eq.ii1 .and. dbady.lt.zero) ishf=-ii1
        if(ishf.eq.ii2 .and. dbadz.lt.zero) ishf=-ii2
!
      ijkshf=ijk_str2unstr(ijk_str+ishf)
      nob=ijkmob(ijkshf)
!
      if(nob.eq.0) then
        rgh=zero
        goto 500
      else
        rgh=rough_cal(ijkshf,nob,2)
      end if
!
! skip free-slip and porous components
      if(rgh.lt.-ztest .or.  opor(nob).gt.eps) goto 1100
      if(imo(nob).eq.0 .and. ifspin(nob).eq.0) goto 500
      ubdy=zero
      vbdy=zero
      wbdy=zero
      if(iabs(ishf).eq.1) then
! ------ boundary in x direction
        area=one-afr(ijk)
        if(ishf.lt.0) area=one-afr(imjk)
        vcsi=vcsi+area
        wcsi=wcsi+area
        if(imo(nob).eq.2) then
          vbdy=velobs(xi(i),yj(j),zk(k),2,nob)
          wbdy=velobs(xi(i),yj(j),zk(k),3,nob)
        elseif(imo(nob).eq.3 .or. imo(nob).eq.4) then
          call mvbuvw(0,xi(i),yj(j),zk(k),cthetj(j),sthetj(j), & 
     &         nob,ubdy,vbdy,wbdy,cyl)
        endif
        vcs=vcs+area*vbdy
        wcs=wcs+area*wbdy
      elseif(iabs(ishf).eq.ii1) then
! ------ boundary in y direction
        area=one-afb(ijk)
        if(ishf.lt.0) area=one-afb(ijmk)
        ucsi=ucsi+area
        wcsi=wcsi+area
        if(imo(nob).eq.2) then
          ubdy=velobs(xi(i),yj(j),zk(k),1,nob)
          wbdy=velobs(xi(i),yj(j),zk(k),3,nob)
        elseif(imo(nob).eq.3 .or. imo(nob).eq.4) then
          call mvbuvw(0,xi(i),yj(j),zk(k),cthetj(j),sthetj(j), & 
     &         nob,ubdy,vbdy,wbdy,cyl)
        endif
        ucs=ucs+area*ubdy
        wcs=wcs+area*wbdy
      else
! ------ boundary in z direction
        area=one-aft(ijk)
        if(ishf.lt.0) area=one-aft(ijkm)
        ucsi=ucsi+area
        vcsi=vcsi+area
        if(imo(nob).eq.2) then
          ubdy=velobs(xi(i),yj(j),zk(k),1,nob)
          vbdy=velobs(xi(i),yj(j),zk(k),2,nob)
        elseif(imo(nob).eq.3 .or. imo(nob).eq.4) then
          call mvbuvw(0,xi(i),yj(j),zk(k),cthetj(j),sthetj(j), & 
     &         nob,ubdy,vbdy,wbdy,cyl)
        endif
        ucs=ucs+area*ubdy
        vcs=vcs+area*vbdy
      endif
      if(ifspin(nob).eq.1) then
        ucsi=ucsi+area
        vcsi=vcsi+area
        wcsi=wcsi+area
        if(cyl.lt.half) & 
     &  ubdy=ubdy+ospinv(nob)*((zk(k)-cmoz(nob))*amoy(nob) & 
     &                  -(yj(j)-cmoy(nob))*amoz(nob))
        vbdy=vbdy+ospinv(nob)*((xi(i)-cmox(nob))*amoz(nob) & 
     &                   -(zk(k)-cmoz(nob))*amox(nob))
        wbdy=wbdy+ospinv(nob)*((yj(j)-cmoy(nob))*amox(nob) & 
     &                   -(xi(i)-cmox(nob))*amoy(nob))
        ucs=ucs+area*ubdy
        vcs=vcs+area*vbdy
        wcs=wcs+area*wbdy
      endif
!
  500 continue
!
! ------ set relative velocities at sliding walls
!
      if(ucsi.gt.ztest) then
         ucs=ucs/ucsi
         uc=uc-ucs
      endif
      if(vcsi.gt.ztest) then
         vcs=vcs/vcsi
         vc=vc-vcs
      endif
      if(wcsi.gt.ztest) then
         wcs=wcs/wcsi
         wc=wc-wcs
      endif
!
      udotn=uc*dbadx+vc*dbady+wc*dbadz
!
      utan=uc-udotn*dbadx*banrms
      vtan=vc-udotn*dbady*banrms
      wtan=wc-udotn*dbadz*banrms
!
!               estimate distance from wall using
!                d = 0.5* length n dot (delx*i +dely*j +delz*k)
!
      dsdot=(bar-bal)*(bar-bal)+(babk-baf)*(babk-baf)+ & 
     &                       (bat-bab)*(bat-bab)
      if(dsdot.lt.ztest) go to 1000
!
!    evaluate turbulent energy and dissipation with log wall function
!
      delta=half*sqrt(banrms*dsdot)*vf(ijk)
      effnu=mucal(ijk,1)/rho(ijk)
      utmag=sqrt(utan*utan+vtan*vtan+wtan*wtan)
      call logwall(rgh,utmag,delta,effnu,ux,yplus(ijk))
!
      if(ifvis.lt.0) goto 1100 ! in case of LES we get yplus out of this routine
!
      us=ux*utmag
      wstsq=us*us
!
      if(iqsr.ne.0) then
        if(qsr(ijk).gt.zero) wstsq=wstsq+two*utmag* & 
     &               vf(ijk)*delta*qsr(ijk)
      endif
!
      tke(ijk)=wstsq/sqrt(cnu)
      if(ifvis.ge.2) then
        if(ifvis.lt.5) then
          dtke(ijk)=max(ztest,wstsq**1.5/(ckar*delta))
        elseif(ifvis.eq.5) then
          dtke(ijk)=max(ztest,sqrt(wstsq/cnu)/(ckar*delta))
        endif
      endif
      goto 1100
!                  no-normal exit
 1000 continue
!
      if(ifvis.gt.0) then
        tke(ijk)=max(tke(ijk),tkemn)
        if(ifvis.lt.5) then
          dtke(ijk)=max(ztest,cnutq*tke(ijk)**x1p5/tlen)
        elseif(ifvis.eq.5) then
          dtke(ijk)=max(ztest,cnunq*sqrt(tke(ijk))/tlen)
        endif
      endif
!
 1100 continue
!
      return
      end
