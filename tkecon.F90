!> \brief 
!>    Convert buoyant and production terms to appropriate energy and dissipation sources
!>          9/86,  rng extension 6/93
!
      subroutine tkecon(i,j,k,avgspr,avgboy)
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
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/dparam.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/edit.f90'
      include '../comdeck/trbcom.f90'
!
      include '../comdeck/ijk.f90'
!
! Limit production term so that it can never go below zero. (could only
! occur in surface cells)
!
      if(nf(ijk).gt.0.and.nf(ijk).lt.8) tkep(ijk)=max(zero,tkep(ijk))
!
      if(tken(ijk).gt.tiny(one) .and. dtken(ijk).gt.tiny(one) .and. & 
     &     mu(ijk).gt.tiny(one)) then
!
        vcell=vf(ijk)*delx(i)*dely(j)/rri(i)*delz(k)
!
        if(ifvis.le.1) then
!
!         mixing length model (production=decay)
!
! Note: tkep(ijk) may become negative at free surface because of the
!       way the shear stress components are computed and summed up.
!
          ptke=cspro*max(zero,tkep(ijk))/(vf(ijk)*mu(ijk))
          gtke=crho*tkeb(ijk)/(delt*delt)
          tke(ijk)=(ptke+gtke)*tlen*tlen/sqrt(cnu)
          if(icstat(ijk_str).gt.0) then
            vcnu=vcell*mu(ijk)/rho(ijk)
            avgspr=avgspr+vcnu*ptke
            avgboy=avgboy+vcnu*gtke
          endif
!
        else
!
!                 turbulent energy transport model
!
! Note: tkep(ijk) may become negative at free surface because of the
!       way the shear stress components are computed and summed up.
!
          ptke=cspro*max(zero,tkep(ijk))/(rho(ijk)*vf(ijk))
          if(icmprs.eq.1) then
            ptke=ptke-cspro*(two/three)*tken(ijk)*div(ijk)
! compressible term can make production negative for expanding gas
            ptke=max(zero,ptke)
          endif
          gtke=crho*mu(ijk)*tkeb(ijk)/(rho(ijk)*delt*delt)
!
          if(ifvis.lt.5) then
!           k-e based (k-e and RNG)
            tke(ijk)=tken(ijk)+delt*(ptke+gtke-dtken(ijk))
          elseif(ifvis.eq.5) then
!           k-w
!           Wilcox (1998)
!
! inverse of dtken is being taken to the power of 3 in chik
! and then chik is taken to the power of 2 in fbtak
!
            if(dtken(ijk).gt.ztest) then
              call slgrad(tken,i,j,k,leneax,0, & 
     &                    gradkx,gradky,gradkz,agradk)
              call slgrad(dtken,i,j,k,leneax,0, & 
     &                    gradwx,gradwy,gradwz,agradw)
              chik=(gradkx*gradwx+gradky*gradwy+ & 
     &                            gradkz*gradwz)
              if(chik.lt.ztest) then
                fbtak(ijk)=one
              else
                chik=chik/(dtken(ijk)**3)
                fbtak(ijk)=(one+680.0d0*chik**2)/(one+400.0d0*chik**2)
              endif
            else
              fbtak(ijk)=680.0d0/400.0d0
            endif
            tke(ijk)=tken(ijk)+delt*(ptke+gtke- & 
     &                         cnu*fbtak(ijk)*dtken(ijk)*tken(ijk))
          endif
!
! do not subtract dissipation if it makes turbulent energy drop below tkemn
!
          if(tke(ijk).lt.tkemn) tke(ijk)=tken(ijk)+delt*(ptke+gtke)
!
          if(icstat(ijk_str).gt.0) then
              avgspr=avgspr+rho(ijk)*vcell*ptke
              avgboy=avgboy+rho(ijk)*vcell*gtke
          endif
!
!
          if(ifvis.ge.3) then
!
!             dissipation terms
!
            drq=one/ttime(ijk)
! ------- check for rng extension to k-e turbulence model
            if(ifvis.lt.5) then
! ###         drq=dtken(ijk)/tken(ijk)
              if(ifvis.eq.3) then
                crng=cdis2
              else
                eta=two*esqr(ijk)/drq
                if(eta.lt.1000.0d0) then
                    etac=eta**3
                    crng=cdis2+cnu*etac*(one-0.2283d0*eta)/ & 
     &                                  (one+0.012d0*etac)
                else
                    crng=cdis2+cnu*(one-0.2283d0*eta)*83.3333d0
                endif
              endif
              if(crng.lt.ztest) then
                dtke(ijk)=dtken(ijk)+delt*drq*(cdis1*(ptke+cdis3*gtke) & 
     &                    -crng*dtken(ijk))
              else
                dtke(ijk)=(dtken(ijk)+delt*cdis1*drq*(ptke+cdis3*gtke))/ & 
     &                    (one+delt*crng*drq)
              endif
            elseif(ifvis.eq.5) then
!
!           k-w
!
! the buoyancy term is derived from buoyancy in k and e equations of k-e
! the effects of cdis3 is transferred from our own k-e derivation which refers
! to the work of Rodi (1980)
!
! ###         drq=cnu*fbtak(ijk)*dtken(ijk)
!
!             Wilcox (1998)
              if(dtken(ijk).gt.ztest) then
                call dvcalcc(i,j,k, & 
     &                   dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
!
! we are calculating the expression:
!
!                chiw=abs((ob8*(((dudy-dvdx)*(dvdz-dwdy)*(dwdx+dudz))+
!     &                         ((dvdz-dwdy)*(dwdx-dudz)*(dudy+dvdx))+
!     &                         ((dwdx-dudz)*(dudy-dvdx)*(dvdz+dwdy))))
!     &                         /(cnu*dtken(ijk))**3)
!                fbtaw=(one+70.0d0*chiw)/(one+80.0d0*chiw)
!
! if one works out the algebra, this reduces to:
                fbtaw=(one+(70.0d0*ob8/(cnu*dtken(ijk))**3)*abs( & 
     &                six*(dudy*dvdz*dwdx+dudz*dvdx*dwdy)-two*( & 
     &                     dudy*dwdx*dwdy+dudz*dvdx*dvdz+ & 
     &                     dvdx*dwdx*dwdy+dvdz*dudz*dudy+ & 
     &                     dwdx*dvdx*dvdz+dwdy*dudy*dudz)))/ & 
     &                (one+(80.0d0*ob8/(cnu*dtken(ijk))**3)*abs( & 
     &                six*(dudy*dvdz*dwdx+dudz*dvdx*dwdy)-two*( & 
     &                     dudy*dwdx*dwdy+dudz*dvdx*dvdz+ & 
     &                     dvdx*dwdx*dwdy+dvdz*dudz*dudy+ & 
     &                     dwdx*dvdx*dvdz+dwdy*dudy*dudz)))
              else
                fbtaw=70.0d0/80.0d0
              endif
              dtke(ijk)=dtken(ijk)+delt*drq/cnu/fbtak(ijk)* & 
     &                  (alphkw/tken(ijk)*(ptke+cdis3*gtke)- & 
     &                  btaw*fbtaw*drq/cnu/fbtak(ijk))
            endif
!
          endif
        endif
      else
        tke(ijk)=tkemn
        dtke(ijk)=tkemn
      endif
!
      return
      end
