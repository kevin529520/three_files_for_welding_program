!> \brief 
!>    Evaluate wall heat flux and heat transfer coefficients
!
!> \param[in] nob         Component number for m>0
!>                        mesh boundary number for m<0
!> \param[in] m           Cell number around component nob (-1 if not relavent)
!> \param[in] ijk 
!> \param[in] twin        Wall temperature
!> \param[in] tfin        Fluid temperature
!> \param[in]  fin        Fluid fraction
!> \param[in]   sa        Wall surface area
!> \param[in]  hi1        Input heat transfer coefficient, to fluid #1
!> \param[in]  hi2        Input heat transfer coefficient, to fluid #2
!> \param[in] scalet      Length scale for correlations
!> \param[in] htcobs      Conduction heat transfer coefficient between obstacle surface and tw node in obstacle
!> \param[in] htcfl       Conduction heat transfer coefficient between obstacle surface and tn node in fluid
!> \param[out]  hasrct    Area-heat transfer coefficient product
!
      subroutine qeval(nob,m,ijk,twin,tfin,fin,sa,hi1,hi2,scalet, & 
     &                 htcobs,htcfl,hasrct)
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
      use fsidata_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/dparam.f90'
      include '../comdeck/acctrm.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/phiou.f90'
      include '../comdeck/state.f90'
!
      real(kind(zzzz)) :: nulam,nuturb,nunat
!
      hasrct=zero
!
      if(nmat.eq.1 .and. fin.lt.emf) return
!
      h1=hi1
!
      if(nmat.eq.1) then
        h2=zero
      else
        h2=hi2
      endif
      if_calc=0
!
      if(h1.gt.ztest) then
!
! treat a given h's as contact resistance
!
          denom=htcobs+htcfl+htcobs*htcfl/h1
          h1=htcobs*htcfl/denom
!
      elseif(h1.lt.-ztest .and. scalet.gt.ztest .and. fin.gt.emf) then
!
!       evaluate htc based on fluid and obstacle conduction
!
          denom=htcfl+htcobs
          hcond=htcfl*htcobs/denom
!
!    use correlations to evaluate htc for fluid #1
!
          if(ifvelp.eq.2) then
              h1=hcond
          else
              if(ifsld.gt.0) then
                  if(fliq(ijk).lt.omfscr) then
                    h1=hcond
                    goto 200
                  endif
              endif
!
              if_calc=1
              tfilm=half*(tfin+twin)
              dtabs=abs(twin-tfin)
              include '../comdeck/ijk2str.f90'
              include '../comdeck/mijk.f90'
! calculate cell-centered fluid velocity
              ufl=(afr(imjk)*u(imjk)+afr(ijk)*u(ijk))/ & 
     &             (afr(imjk)+afr(ijk)+ztest)
              vfl=(afb(ijmk)*v(ijmk)+afb(ijk)*v(ijk))/ & 
     &             (afb(ijmk)+afb(ijk)+ztest)
              wfl=(aft(ijkm)*w(ijkm)+aft(ijk)*w(ijk))/ & 
     &             (aft(ijkm)+aft(ijk)+ztest)
! calculate wall velocity
              ubdy=zero
              vbdy=zero
              wbdy=zero
              i=i_str(ijk_str)
              j=j_str(ijk_str)
              k=k_str(ijk_str)
              if(m.lt.0) then
! at a mesh boundary
                if(nob.ne.1 .and. nob.ne.2) ubdy=ubc(nob,nbl)
                if(nob.ne.3 .and. nob.ne.4) vbdy=vbc(nob,nbl)
                if(nob.ne.5 .and. nob.ne.6) wbdy=wbc(nob,nbl)
              elseif(m.gt.0) then
! at a component surface
                if(imo(nob).eq.2) then
                  ubdy=velobs(xi(i),yj(j),zk(k),1,nob)
                  vbdy=velobs(xi(i),yj(j),zk(k),2,nob)
                  wbdy=velobs(xi(i),yj(j),zk(k),3,nob)
                elseif(imo(nob).ge.3) then
                  call mvbuvw(0,xi(i),yj(j),zk(k),cthetj(j),sthetj(j), & 
     &         nob,ubdy,vbdy,zbdy,cyl)
                end if
! ----- get velocity of rotating components
                if(ifspin(nob).eq.1 .and. cyl.lt.half) then
                  ubdy=ubdy+ospinv(nob)*((zk(k)-cmoz(nob))*amoy(nob) & 
     &                                  -(yj(j)-cmoy(nob))*amoz(nob))
                  vbdy=vbdy+ospinv(nob)*((xi(i)-cmox(nob))*amoz(nob) & 
     &                                  -(zk(k)-cmoz(nob))*amox(nob))
                  wbdy=wbdy+ospinv(nob)*((yj(j)-cmoy(nob))*amox(nob) & 
     &                                  -(xi(i)-cmox(nob))*amoy(nob))
                endif
              endif
!
              vbar=sqrt((ufl-ubdy)*(ufl-ubdy)+(vfl-vbdy)*(vfl-vbdy)+ & 
     &                  (wfl-wbdy)*(wfl-wbdy))
!
              if(iaccf.gt.0) then
                  call accxcl(i,j,k,u,v,w,leneax,icorx,accx)
                  call accycl(i,j,k,u,v,w,leneax,icory,accy)
                  call acczcl(i,j,k,u,v,leneax,icorz,accz)
              else
                  accx=zero
                  accy=zero
                  accz=zero
              endif
              grav=sqrt((gx-accx)**2+(gy-accy)**2+(gz-accz)**2)
!
              den1=rho1cl(p(ijk),tn(ijk))
              ren1=den1*vbar*scalet/mu1ht
              sqre=sqrt(ren1)
              rep8=ren1**0.8d0
              prp33=prn1**ob3
              nulam=0.664d0*prp33*sqre
              nuturb=prp33*(0.036d0*rep8-836.d0)
              nunat=zero
!
              if(grav.gt.ztest .and. ifvelp.ne.1) then
                  den1f=rho1cl(p(ijk),tfilm)
                  beta1=-dr1dt(p(ijk),tfilm)/den1f
                  grsh1=grav*beta1*dtabs*scalet**3*(den1f/mu1ht)**2
                  grpr=abs(grsh1*prn1)
                  if(grpr.le.1.0d+9) nunat=0.57d0*grpr**quarter
                  if(grpr.gt.1.0d+9) nunat=0.14d0*grpr**ob3
              endif
!
              if(ithc1.eq.1) then
                thc1t=alookup(tfin,thc1_var1,thc1_var2,thc1_delta, & 
     &                        thc1_deltar,thc1_t,ntabmx)
              else
                thc1t=thc1
              endif
              h1=max(nulam,nuturb,nunat)*thc1t/scalet
              h1=htcobs*h1/(h1+htcobs)
              h1=max(h1,hcond)
          endif
      else
          h1=zero
      endif
!
  200 continue
!
! Adjust heat transfer coefficient for gap formation if used in the Thermal Stress
! Analysis model
! m here is the component boundary index. It's -1 when qeval is called
! for domain boundaries.
      if(m.gt.0 .and. ithrmst.gt.0 .and. thcgap.gt.-ztest) then
! If the thermal stress model is active, modify the heat transfer coefficient
! due to gap formation
        if(im2face(m).gt.0) then
         gap_length=xm2face(m)
         if(gap_length.gt.ztest) then
          htc_tot=thcgap/gap_length
          if(sbconst.gt.ztest .and. & 
     &            efe1gap .gt.ztest .and. eobegap(nob).gt.ztest) then
            if(one/efe1gap+one/eobegap(nob)-one.gt.em6) then
              t_mean=half*(twin+tfin)
              t_diff=tfin-twin
              htc_rad=sbconst*t_mean*(four*(t_mean**2)+t_diff**2) & 
     &                 /(one/efe1gap+one/eobegap(nob)-one)
              denom=htc_tot+htc_rad
              if(denom.gt.ztest) htc_tot=htc_tot*htc_rad/denom
            endif
          endif
          denom=h1+htc_tot
          if(denom.gt.ztest) h1=h1*htc_tot/denom
         endif
        endif
      endif
!
      if(nmat.eq.1) goto 300
!
      if(h2.gt.ztest) then
!
! treat a given h's as contact resistance
!
          denom=htcobs+htcfl+htcobs*htcfl/h2
          h2=htcobs*htcfl/denom
!
      elseif(h2.lt.-ztest .and. scalet.gt.ztest .and. fin.lt.emf1) then
!
!       evaluate htc based on fluid and obstacle conduction
!
          denom=htcfl+htcobs
          hcond=htcfl*htcobs/denom
!
!    use correlations to evaluate htc for fluid #2
!
          if(ifvelp.eq.2) then
              h2=hcond
          else
! do the calculations if they have not been done for fluid #1 above
              if(if_calc.eq.0) then
                tfilm=half*(tfin+twin)
                dtabs=abs(twin-tfin)
                include '../comdeck/ijk2str.f90'
                include '../comdeck/mijk.f90'
! calculate cell-centered fluid velocity
                ufl=(afr(imjk)*u(imjk)+afr(ijk)*u(ijk))/ & 
     &               (afr(imjk)+afr(ijk)+ztest)
                vfl=(afb(ijmk)*v(ijmk)+afb(ijk)*v(ijk))/ & 
     &               (afb(ijmk)+afb(ijk)+ztest)
                wfl=(aft(ijkm)*w(ijkm)+aft(ijk)*w(ijk))/ & 
     &               (aft(ijkm)+aft(ijk)+ztest)
! calculate wall velocity
                ubdy=zero
                vbdy=zero
                wbdy=zero
                i=i_str(ijk_str)
                j=j_str(ijk_str)
                k=k_str(ijk_str)
                if(m.lt.0) then
! at a mesh boundary
                  if(nob.ne.1 .and. nob.ne.2) ubdy=ubc(nob,nbl)
                  if(nob.ne.3 .and. nob.ne.4) vbdy=vbc(nob,nbl)
                  if(nob.ne.5 .and. nob.ne.6) wbdy=wbc(nob,nbl)
                elseif(m.gt.0) then
! at a component surface
                  if(imo(nob).eq.2) then
                    ubdy=velobs(xi(i),yj(j),zk(k),1,nob)
                    vbdy=velobs(xi(i),yj(j),zk(k),2,nob)
                    wbdy=velobs(xi(i),yj(j),zk(k),3,nob)
                  elseif(imo(nob).ge.3) then
                    call mvbuvw(0,xi(i),yj(j),zk(k),cthetj(j),sthetj(j), & 
     &           nob,ubdy,vbdy,zbdy,cyl)
                  end if
! ----- get velocity of rotating components
                  if(ifspin(nob).eq.1 .and. cyl.lt.half) then
                    ubdy=ubdy+ospinv(nob)*((zk(k)-cmoz(nob))*amoy(nob) & 
     &                                    -(yj(j)-cmoy(nob))*amoz(nob))
                    vbdy=vbdy+ospinv(nob)*((xi(i)-cmox(nob))*amoz(nob) & 
     &                                    -(zk(k)-cmoz(nob))*amox(nob))
                    wbdy=wbdy+ospinv(nob)*((yj(j)-cmoy(nob))*amox(nob) & 
     &                                    -(xi(i)-cmox(nob))*amoy(nob))
                  endif
                endif
!
                vbar=sqrt((ufl-ubdy)*(ufl-ubdy)+(vfl-vbdy)*(vfl-vbdy)+ & 
     &                    (wfl-wbdy)*(wfl-wbdy))
!
                if(iaccf.gt.0) then
                    call accxcl(i,j,k,u,v,w,leneax,icorx,accx)
                    call accycl(i,j,k,u,v,w,leneax,icory,accy)
                    call acczcl(i,j,k,u,v,leneax,icorz,accz)
                else
                    accx=zero
                    accy=zero
                    accz=zero
                endif
                grav=sqrt((gx-accx)**2+(gy-accy)**2+(gz-accz)**2)
              endif
!
              if(incg.gt.0) then
                den2=rho2cl(p(ijk),tn(ijk),rf2t(ijk))
              else
                den2=rho2cl(p(ijk),tn(ijk),rf2)
              endif
              ren2=den2*vbar*scalet/mu2ht
              sqre=sqrt(ren2)
              rep8=ren2**0.8d0
              prp33=prn2**ob3
              nulam=0.664d0*prp33*sqre
              nuturb=prp33*(0.036d0*rep8-836.d0)
              nunat=zero
!
              if(grav.gt.ztest .and. ifvelp.ne.1) then
                  if(incg.gt.0) then
                    den2f=rho2cl(p(ijk),tfilm,rf2t(ijk))
                  else
                    den2f=rho2cl(p(ijk),tfilm,rf2)
                  endif
                  beta2=-dr2dt(p(ijk),tfilm)/den2f
                  grsh2=grav*beta2*dtabs*scalet**3*(den2f/mu2ht)**2
                  grpr=abs(grsh2*prn2)
                  if(grpr.le.1.0d+9) nunat=0.57d0*grpr**quarter
                  if(grpr.gt.1.0d+9) nunat=0.14d0*grpr**ob3
              endif
!
              if(ithc2.eq.1) then
                thc2t=alookup(tfin,thc2_var1,thc2_var2,thc2_delta, & 
     &                        thc2_deltar,thc2_t,ntabmx)
              else
                thc2t=thc2
              endif
              h2=max(nulam,nuturb,nunat)*thc2t/scalet
              h2=htcobs*h2/(h2+htcobs)
              h2=max(h2,hcond)
          endif
      else
          h2=zero
      endif
!
  300 continue
!
      hasrct=sa*((one-fin)*h2+fin*h1)
!
      return
      end
