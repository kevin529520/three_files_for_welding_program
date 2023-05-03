!> \brief 
!>    Evaluate mass source term QSRM, 
!>    bulk volume source term QSR,
!>    fluid fraction source term QSRF,
!>    thermal energy source term QSRME.
!
      subroutine srcal
!
!         qsrm has unit of density per unit time,
!         qsr and qsrf have unit of 1/time,
!         qsrme has the unit of power per unit volume.
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
      use arrayp_module
!
      use arraymm_module
!
      use meshcb_module
!
      use voids_module
!
      use mblock_module
!
      use obsijk_module
!
      use procs_module
!
      use probe_event_module
!
#ifdef MPI
      use parallel_module
#endif
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/obsd.f90'
      include '../comdeck/diag.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/state.f90'
      include '../comdeck/pardat.f90'
      include '../comdeck/scala.f90'
      include '../comdeck/phiou.f90'
      include '../comdeck/valves.f90'
      include '../comdeck/trbcom.f90'
#ifdef MPI
      include '../comdeck/dparam.f90'
      integer :: nnbl, block_owner
      LOGICAL,save::initflagsrcal   ! controls dynamic array allocation
      DATA          initflagsrcal/.TRUE./
      REAL(kind(zzzz)),save,allocatable :: vsrobs_mpi(:)
      REAL(kind(zzzz)),save,allocatable :: msrobs_mpi(:)
      REAL(kind(zzzz)),save,allocatable :: sendbuf1(:),recvbuf1(:)
      REAL(kind(zzzz)) :: sendbuf(nvorx)
!
! Allocate local MPI variables in the first call
!
      if(initflagsrcal) then
         allocate(vsrobs_mpi(nobs))
         allocate(msrobs_mpi(nobs))
         allocate(sendbuf1(2*nobs))
         allocate(recvbuf1(2*nobs))
         sendbuf1=zero
         recvbuf1=zero
         initflagsrcal=.FALSE.
      endif
!
! Initialize to zero only for the first block on each MPI rank
! (nmybl==1)
!
      if(nmybl.eq.1)then
         vsrobs_mpi(1:nobs)=zero
         msrobs_mpi(1:nobs)=zero
      endif
!
      sendbuf=zero
      recvbuf=zero
#endif
!
      rdeltc=one/deltc
#ifdef MPI
      if(nmybl.eq.1) then
        do l=8,nreg+7
           dvrp_mpi(l)=zero
        enddo
      endif
#endif
!
!     loop over obstacles
!
      do 1000 nob=1,nobs
!
        if(iaqsrb(nob).ne.1) cycle
        if(saobs(nob).lt.ztest .or. & 
     & (abs(qmob(nob)).lt.ztest .and. imo(nob).ne.7)) cycle
!
            mincel=kajk(nob,nbl,1)
            maxcel=kajk(nob,nbl,2)
!
            qbsa=qmob(nob)/saobs(nob)
!
                do 350 m=mincel,maxcel
!
                   ijk=ijkobs(m)
                   if(ijk.lt.1) cycle
                   if(vf(ijk).lt.em6) cycle
!
                   if(imo(nob).eq.7) then
                     if(vf0(ijk).lt.omeps) cycle ! surface area in mixed cells is not accurate
                   endif
!
                   include '../comdeck/ijk2str.f90'
                   i=i_str(ijk_str)
                   j=j_str(ijk_str)
                   k=k_str(ijk_str)
                   vcell=vf(ijk)*delx(i)*dely(j)/rri(i)*delz(k)
                   sa=waobs(m)
!
                   if(imo(nob).eq.7) then
! combustible solid
                     sa=sa*(one-f(ijk))
                     burn_area(ijk)=burn_area(ijk)+sa
                     qbsa=zero
                     tqsrb_l=tn(ijk)
                     rhoqsrb_l=rho2cl(p(ijk),tqsrb_l,rf2)
!
                     if(p(ijk).gt.ztest .and. tn(ijk).gt.ztest) then
! calculate reaction mass rate qbsa
                       call combust_sld(i,j,k,ijk,nob,qbsa,tqsrb_l,rhoqsrb_l)
                     endif
!
                   else
                     rhoqsrb_l=rhoqsrb(nob)
                     tqsrb_l=tqsrb(nob)
                   endif
!
                   qsrt=sa*qbsa
!
!  calculate qsr(ijk), qsrm(ijk), qsrf(ijk)
!
                   dqsr=zero    ! increase of volume per unit open volume in mesh cell per unit time
                   dqsrm=zero   ! increase of mass per unit open volume in mesh cell per unit time
                   dqsrf=zero   ! increase of fluid fraction per unit open volume in mesh cell per unit time
!
                   if(qsrt.gt.ztest) then ! for source
                       if(nmat.eq.1 .and. ifqsrb(nob).eq.2) then ! if the source release bubble
                         dqsrm=zero
                         dqsrf=zero
                         if(ivqsrb(nob).eq.1) then ! volume flow rate is defined
                           dqsr=qsrt/vcell
                         else if(ivqsrb(nob).eq.2) then ! mass flow rate is defined
                           dqsr=qsrt/(vcell*rhoqsrb_l)
                         end if
                       else
                         if(ivqsrb(nob).eq.1) then ! volume flow rate is defined
                           dqsr=qsrt/vcell
                           dqsrm=dqsr*rhoqsrb_l
                           dqsrf=dqsr
                         else if(ivqsrb(nob).eq.2) then ! mass flow rate is defined
                           dqsrm=qsrt/vcell
                           dqsr=dqsrm/rhoqsrb_l
                           dqsrf=dqsr
                         end if
                       end if
                       if(nmat.eq.2 .and. ifqsrb(nob).eq.2) dqsrf=zero
                   elseif(qsrt.lt.-ztest) then ! for sink
                       if(ivqsrb(nob).eq.1) then ! volume flow rate is defined
                         dqsr=qsrt/vcell
                         if(ifrho.ge.2) then
                           dqsrm=dqsr*rho(ijk)
                         else
                           dqsrm=dqsr*rhof
                         endif
                         dqsrf=dqsr
                       else if(ivqsrb(nob).eq.2) then ! mass flow rate is defined
                         dqsrm=qsrt/vcell
                         if(nmat.eq.1 .and. ifsnkb(nob).eq.2) then  ! sink of void for one-fluid problem
                           dqsr=dqsrm/rhosnkb(nob)
                         else if(nmat.eq.1 .and. ifsnkb(nob).eq.3) then  ! sink of mixture of fluid
                           rhomix=f(ijk)*rho(ijk)+(one-f(ijk))* & 
     &                            rhosnkb(nob)
                           dqsr=dqsrm/rhomix
                         else
                           dqsr=dqsrm/rho(ijk)           ! all other cases
                         end if
                         dqsrf=dqsr
                       end if
!
                       if(nmat.eq.1) then  ! one fluid
                         if(ifvof.gt.0 .and. giswat(nbl).eq.0) then
                             if(ifsnkb(nob).eq.1) then
                               dqsrf=f(ijk)*dqsrf
                               if(dcor.lt.half) dqsrm=f(ijk)*dqsrm
                             else if(ifsnkb(nob).eq.2) then
                               dqsr=(one-f(ijk))*dqsr
                               dqsrm=zero
                               dqsrf=zero
                             else if(ifsnkb(nob).eq.3) then  ! mixture
                               dqsrm=f(ijk)*dqsrm
                               dqsrf=f(ijk)*dqsrf
                             end if
                         end if
!
                       else   ! two fluids
                         if(icmprs.eq.0) then
                           if(ifsnkb(nob).eq.3) then
                             dqsrf=f(ijk)*dqsrf
                           else if(ifsnkb(nob).eq.1) then
                             dqsrm=f(ijk)*dqsrm ! fluid mass change rate for fluid #1
                             dqsr=dqsrm/rhof    ! fluid volume change rate for fluid #1
                             dqsrf=dqsr
                           else if(ifsnkb(nob).eq.2) then
                             dqsrm=(one-f(ijk))*dqsrm  ! fluid mass change rate for fluid #2
                             dqsr=dqsrm/rhof2          ! fluid volume change rate for fluid #2
                             dqsrf=zero
                           end if
                         else if(icmprs.eq.1) then
                           if(ivqsrb(nob).eq.2) then   ! mass flow rate is defined
                             if(ifsnkb(nob).eq.3) then
                               dqsrf=f(ijk)*dqsrf
                             else if(ifsnkb(nob).eq.1) then
                               dqsrm=f(ijk)*dqsrm      ! fluid mass change rate for fluid #1
                               dqsr=dqsrm/rhof         ! fluid volume change rate for fluid #1
                               dqsrf=dqsr
                             else if(ifsnkb(nob).eq.2) then
                               if(f(ijk).lt.omepsf) then
                                 dqsrm=(one-f(ijk))*dqsrm    ! fluid mass change rate for fluid #2
                                 rho2=(rho(ijk)-f(ijk)*rhof)/ & 
     &                              (one-f(ijk))
                                 dqsr=dqsrm/rho2  ! fluid volume change rate for fluid #2
                               end if
                               dqsrf=zero
                             end if
                           else if(ivqsrb(nob).eq.1) then  ! volume flow rate is defined
                             if(ifsnkb(nob).eq.3) then
                               dqsrf=f(ijk)*dqsrf
                             else if(ifsnkb(nob).eq.1) then
                               dqsr=f(ijk)*dqsr   ! fluid mass change rate for fluid #1
                               dqsrm=rhof*dqsr    ! fluid volume change rate for fluid #1
                               dqsrf=dqsr
                             else if(ifsnkb(nob).eq.2) then
                               dqsrm=(rho(ijk)-f(ijk)*rhof)*dqsr !mass change rate for fluid #2
                               dqsr=(one-f(ijk))*dqsr
                               dqsrf=zero
                             end if
                           end if
                         end if
                       end if
                   end if
!
                   qsr(ijk)=qsr(ijk)+dqsr
                   qsrm(ijk)=qsrm(ijk)+dqsrm
                   if(ifvof.gt.0) qsrf(ijk)=qsrf(ijk)+dqsrf
!
                   ifqsrsnk=0
                   if((ifqsrb(nob).eq.1.and.dqsr.gt.ztest) & 
     &                   .or.(ifsnkb(nob).eq.1.and.dqsr.lt.-ztest)) then
                      ifqsrsnk=1 ! source or sink is for fluid #1
                   else if((ifqsrb(nob).eq.2.and.dqsr.gt.ztest) & 
     &                   .or.(ifsnkb(nob).eq.2.and.dqsr.lt.-ztest)) then
                      ifqsrsnk=2 ! source or sink is for fluid #2
                   else if(dqsr.lt.-ztest.and.ifsnkb(nob).eq.3) then
                      ifqsrsnk=3 ! sink is for mixture of fluids #1 and #2
                   end if
!
                   if(ifenrg.ge.2) then
                     dqsrme=zero
                     if(dqsr.lt.zero) then ! sink
                       temp=tn(ijk)
                     else                  ! source
                       temp=tqsrb_l
                     end if
                     if(ifqsrsnk.eq.3) then  ! sink is for mixture
                       dqsrme=dqsr*rhoe(ijk)
                       if(icmprs.eq.1) & 
     &                   dqsrme=dqsrme+(one-f(ijk))*dqsr*p(ijk)  ! for compressible flow
                       if(decor.gt.half) dqsrme=dqsrme- & 
     &                      dqsrm*e1cal(p(ijk),tn(ijk),zero,-one)
                     else if(ifqsrsnk.eq.1) then  ! source or sink is for fluid #1
                       dqsrme=dqsrm*e1cal(p(ijk),temp,zero,-one)
                       if(nmat.eq.1 .and. qsrt.lt.-ztest) & 
     &                        dqsrme=dqsrme*f(ijk)
                       if(icmprs.eq.1) then
                         dqsr0=qsrt/(vcell*rhof)
                         dqsrme=dqsrme+(one-f(ijk))*dqsr0*p(ijk)
                       end if
                       if(decor.gt.half) dqsrme=dqsrme- & 
     &                      dqsrm*e1cal(p(ijk),tn(ijk),zero,-one)
                     else if(ifqsrsnk.eq.2 .and. nmat.eq.2) then  ! source or sink is for fluid #2
                       if(icmprs.ne.0) then
                         if(f(ijk).lt.omepsf .and. incg.gt.0) then
                           rf2tt=rf2t(ijk)
                           e2=e2cal_nc(p(ijk),temp,sclr(ijk,incg), & 
     &                                rho(ijk),f(ijk))
                         else
                           rf2tt=rf2
                           e2=e2cal(p(ijk),temp)
                         end if
                       else
                         rf2tt=zero
                         e2=e2cal(p(ijk),temp)
                       end if
                       dqsrme=dqsrm*(e2+rf2tt*temp) ! cp is used instead cv for compressible fluid
                     end if
                     qsrme(ijk)=qsrme(ijk)+dqsrme
                   end if
!
                   if((ifqsrsnk.eq.2 .or. ifqsrsnk.eq.3) & 
     &                                .and. nmat.eq.1) then
                     if(nf(ijk).gt.7) then
#ifdef MPI
                       dvrp_mpi(nf(ijk))=dvrp_mpi(nf(ijk))+gamma* & 
     &                         pr(nf(ijk))*dqsr*vcell*delt
#else
                       dvrp(nf(ijk))=dvrp(nf(ijk))+gamma*pr(nf(ijk))* & 
     &                         dqsr*vcell*delt
#endif
                     elseif(nf(ijk).ge.1 .and. nf(ijk).le.6) then
                       call nvcell(ijk,nvoid,istat)
                       if(nvoid.gt.7) then
! sink or source of pure gas (void)
                         if((qsrt.lt.-ztest.and.ifsnkb(nob).eq.2) .or. & 
     &                              qsrt.gt.ztest) then
#ifdef MPI
                           dvrp_mpi(nvoid)=dvrp_mpi(nvoid)+gamma* & 
     &                       pr(nvoid)*dqsr*vcell*delt
#else
                           dvrp(nvoid)=dvrp(nvoid)+gamma*pr(nvoid)*dqsr* & 
     &                       vcell*delt
#endif
                         else
! sink of gas/liquid mixture
#ifdef MPI
                          dvrp_mpi(nvoid)=dvrp_mpi(nvoid)+gamma* & 
     &                       pr(nvoid)*dqsr*vcell*delt*(one-f(ijk))
#else
                          dvrp(nvoid)=dvrp(nvoid)+gamma*pr(nvoid)*dqsr* & 
     &                       vcell*delt*(one-f(ijk))
#endif
                         end if
                       end if
                     endif
                   endif
!
! Add force due to moving source contribution to fluid momentum:
! F=rhoijk*Vol*qsr/rhot*U_source, rhot=f*rhoijk for nmat=1 and =rhoijk for nmat=2.
! rhoijk cancels in tilde when fsigx,y,z is added. This is consistent with what
! is in drgcl for drg(ijk) for sources
!
                   if(qsrt.gt.ztest) then
!
! calculate fluid velocity Us at source, relative to the source component
!
                     if(ioqsrt(nob).eq.0) then
                       include '../comdeck/mijk.f90'
                       include '../comdeck/pijk.f90'
                       bal=afr(imjk)
                       bar=afr(ijk)
                       if(i.eq.2 .and. wl.le.2) then
                         if(bar.gt.em6) then
                           bal=bar
                         else
                           bal=one
                         endif
                       endif
                       if(i.eq.im1 .and. wr.le.2) then
                         if(bal.gt.em6) then
                           bar=bal
                         else
                           bar=one
                         endif
                       endif
                       dadx=rdx(i)*(bar-bal)
!
                       baf=afb(ijmk)
                       babk=afb(ijk)
                       if(j.eq.2 .and. wf.le.2) then
                         if(babk.gt.em6) then
                           baf=babk
                         else
                           baf=one
                         endif
                       endif
                       if(j.eq.jm1 .and. wbk.le.2) then
                         if(baf.gt.em6) then
                           babk=baf
                         else
                           babk=one
                         endif
                       endif
                       dady=rdy(j)*rri(i)*(babk-baf)
!
                       bab=aft(ijkm)
                       bat=aft(ijk)
                       if(k.eq.2 .and. wb.le.2) then
                         if(bat.gt.em6) then
                           bab=bat
                         else
                           bab=one
                         endif
                       endif
                       if(k.eq.km1 .and. wt.le.2) then
                         if(bab.gt.em6) then
                           bat=bab
                         else
                           bat=one
                         endif
                       endif
                       dadz=rdz(k)*(bat-bab)
                       snorm=dadx*dadx+dady*dady+dadz*dadz
                       if(snorm.gt.ztest) then
                         rsnorm=one/sqrt(snorm)
                         dadx=rsnorm*dadx
                         dady=rsnorm*dady
                         dadz=rsnorm*dadz
                       endif
                       if(ivqsrb(nob).eq.2) then
                         us=dadx*qbsa/rhoqsrb_l
                         vs=dady*qbsa/rhoqsrb_l
                         ws=dadz*qbsa/rhoqsrb_l
                       else
                         us=dadx*qbsa
                         vs=dady*qbsa
                         ws=dadz*qbsa
                       endif
                     else
                       us=zero
                       vs=zero
                       ws=zero
                     endif
!
! calculate source component velocity Us
!
                     if(imo(nob).gt.2) then
                       call mvbuvw(0,xi(i),yj(j),zk(k), & 
     &                      cthetj(j),sthetj(j), & 
     &                      nob,uo,vo,wo,cyl)
                     else
                       uo=zero
                       vo=zero
                       wo=zero
                     endif
!
                     if(nmat.eq.1 .and. f(ijk).gt.emf) then
                       denom=one/f(ijk)
                     else
                       denom=one
                     endif
!
                     if(ivqsrb(nob).eq.1) denom=denom*rhoqsrb_l
!
! fluid absolute velocity at source is Uo+Us
!
                     fsigx(ijk)=fsigx(ijk)+qsrt*(uo+us)*denom
                     fsigy(ijk)=fsigy(ijk)+qsrt*(vo+vs)*denom
                     fsigz(ijk)=fsigz(ijk)+qsrt*(wo+ws)*denom
!
                   endif
!
!          diagnositcs
!
                   if(icstat(ijk_str).gt.0) then
                       delf=dqsr*vcell*delt
                       if(delf.lt.-ztest .and. nmat.eq.1) then
                         delft=delf
                         delf=max(-f(ijk)*vcell,delf)
                         dqsrm=dqsrm*delf/delft
                       endif
                       smass=smass+dqsrm*vcell*delt
#ifdef MPI
                       vsrobs_mpi(nob)=vsrobs_mpi(nob)+delf*rdelt
                       msrobs_mpi(nob)=msrobs_mpi(nob)+dqsrm*vcell
#else
                       vsrobs(nob)=vsrobs(nob)+delf*rdelt
                       msrobs(nob)=msrobs(nob)+dqsrm*vcell
#endif
                       if(ifqsrsnk.eq.1) then
                         sf1=sf1+delf
#ifdef MPI
                         vl1_out_mpi=vl1_out_mpi-min(delf,zero)
                         vl1_in_mpi=vl1_in_mpi+max(delf,zero)
#else
                         vl1_out=vl1_out-min(delf,zero)
                         vl1_in=vl1_in+max(delf,zero)
#endif
                       elseif(ifqsrsnk.eq.2) then
                         sf2=sf2+delf
#ifdef MPI
                         vl2_out_mpi=vl2_out_mpi-min(delf,zero)
#else
                         vl2_out=vl2_out-min(delf,zero)
#endif
                       elseif(ifqsrsnk.eq.3 .and. dqsr.lt.zero) then
                         sf1=sf1+delf*f(ijk)
                         sf2=sf2+delf*(one-f(ijk))
                       endif
                   endif
!
! add sources to scalars
!
                   delf1=dqsrf*delt        ! volume chamge for fluid #1 scalars
                   if(ifqsrsnk.eq.2) then
                       delf2=dqsr*delt     ! volume change for fluid #2 scalars, pure fluid #2 sink
                   else
                       delf2=(one-f(ijk))*dqsr*delt  ! volume change for fluid #2 scalars, mixture sink
                   endif
                   delf3=dqsr*delt
!
                   do ns=1,nsc
                     ff=f(ijk)
                     if((ifqsrsnk.eq.1 .or. ifqsrsnk.eq.3) & 
     &                             .and. isclr(ns).eq.3) then
                       if(ff.gt.emf) then
                         if(delf1.gt.zero) then
                           sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                              sclqsr(ns,nob)*delf1/ff
                           if(dscorsc(ns).gt.half) sclr(ijk,ns)= & 
     &                         sclr(ijk,ns)-sclrn(ijk,ns)*delf1/ff
                         else
                           if(dscorsc(ns).lt.half) sclr(ijk,ns)=   & ! do not apply sink source
     &                         sclr(ijk,ns)+sclrn(ijk,ns)*delf1/ff
                         endif
                       else
! special logic for newly filled cells: goes together with the change in vfsources
                         if(delf1.gt.zero) sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                                sclqsr(ns,nob)*delf1
                       endif
                     elseif((ifqsrsnk.eq.2 .or. ifqsrsnk.eq.3) & 
     &                                   .and.isclr(ns).eq.5) then
!
                       if(ff.lt.emf1) then
                         if(delf2.gt.zero) then
                           if(ns.eq.iburn) then
! update mass fraction of combustion gas
                             if(ff.lt.omepsf) then
                               rho2=(rho(ijk)-f(ijk)*rhof)/(one-f(ijk))
                               sclr(ijk,ns)= & 
     &    (sclr(ijk,ns)*rho2+dqsrm*delt)/(rho2+dqsrm*delt)
                             endif
                           else
                             sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                         sclqsr(ns,nob)*delf2/(one-ff)
                             if(dscorsc(ns).gt.half) sclr(ijk,ns)= & 
     &                         sclr(ijk,ns)-sclrn(ijk,ns)*delf2/(one-ff)
                           endif
                         else
                           if(dscorsc(ns).lt.half) sclr(ijk,ns)=   & ! do not apply sink source
     &                        sclr(ijk,ns)+sclrn(ijk,ns)*delf2/(one-ff)
                         endif
                       else
! special logic for newly filled cells: goes together with the change in vfsources
                         if(delf2.gt.zero) sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                           sclqsr(ns,nob)*delf2
                       endif
                     elseif(isclr(ns).eq.1) then
                       if(delf3.gt.zero) then
                         sclr(ijk,ns)=sclr(ijk,ns)+sclqsr(ns,nob)*delf3
                         if(dscorsc(ns).gt.half) sclr(ijk,ns)= & 
     &                         sclr(ijk,ns)-sclrn(ijk,ns)*delf3
                       else
                         if(dscorsc(ns).lt.half) sclr(ijk,ns)=   & ! do not apply sink source fo
     &                         sclr(ijk,ns)+sclrn(ijk,ns)*delf3
                       endif
                     endif
                   enddo
!
  350           continue
!
 1000 continue
#ifdef MPI
      if(nmybl.eq.myblcks .and. nobs.ge.1) then
         sendbuf1(1:nobs)=vsrobs_mpi(1:nobs)
         sendbuf1(nobs+1:2*nobs)=msrobs_mpi(1:nobs)
         call parallel_sum(sendbuf1(1:2*nobs),recvbuf1)
         vsrobs_mpi(1:nobs)=recvbuf1(1:nobs)
         msrobs_mpi(1:nobs)=recvbuf1(nobs+1:2*nobs)
!
         do nob=1,nobs
            vsrobs(nob)=vsrobs_mpi(nob)+vsrobs(nob)
            msrobs(nob)=msrobs_mpi(nob)+msrobs(nob)
         enddo
      endif
#endif
!
! ------ add mass sources from "particle" defined sources
!
      if(np_source.gt.0) then
!
        do 3200 k1=1,num_particles(np_source)
!
            if(source_particle(k1)%nbl_part.ne.nbl) cycle
!
            mv=source_particle(k1)%istag_part      ! source id
            if(iaqsrp(mv).ne.1) cycle
!
            ijk_str=source_particle(k1)%ijk_str_part
            ijk=ijk_str2unstr(ijk_str)
            if(ijk.ge.ijklim_bc) cycle
            if(vf(ijk).lt.em6) cycle
!
            i=i_str(ijk_str)
            j=j_str(ijk_str)
            k=k_str(ijk_str)
!
! No VF in VCELL since we are distributing the source between the ijk and ijkx cells
! with VF acting as the weighting factor. See this line for ijkx cell below:
!               qsrt=source_particle(k1)%flow_rate*(one-vf(ijk))
!
            vcell=delx(i)*dely(j)/rri(i)*delz(k)
            qsrt=source_particle(k1)%flow_rate
!
!  calculate qsr(ijk), qsrm(ijk), qsrf(ijk)
!
            dqsr=zero
            dqsrm=zero
            dqsrf=zero
!
            ifqsrsnk=0
            if((ifqsrp(mv).eq.1.and.qsrt.gt.ztest) & 
     &            .or.(ifsnkp(mv).eq.1.and.qsrt.lt.-ztest)) then
              ifqsrsnk=1 ! source or sink is for fluid #1
            else if((ifqsrp(mv).eq.2.and.qsrt.gt.ztest) & 
     &            .or.(ifsnkp(mv).eq.2.and.qsrt.lt.-ztest)) then
              ifqsrsnk=2 ! source or sink is for fluid #2
            else if(qsrt.lt.-ztest.and.ifsnkp(mv).eq.3) then
              ifqsrsnk=3 ! sink is for mixture of fluids #1 and #2
            end if
!
            if(qsrt.gt.ztest) then ! for source
                if(nmat.eq.1 .and. ifqsrp(mv).eq.2) then ! if the source release bubble
                  dqsrm=zero
                  dqsrf=zero
                  if(ivqsrp(mv).eq.1) then ! volume flow rate is defined
                    dqsr=qsrt/vcell
                  else if(ivqsrp(mv).eq.2 .or. ivqsrp(mv).eq.3) then ! mass flow rate is defined
                    dqsr=qsrt/(vcell*rhoqsrp(mv))
                  end if
                else
                  if(ivqsrp(mv).eq.1) then ! volume flow rate is defined
                    dqsr=qsrt/vcell
                    dqsrm=dqsr*rhoqsrp(mv)
                    dqsrf=dqsr
                  else if(ivqsrp(mv).eq.2 .or. ivqsrp(mv).eq.3) then ! mass flow rate is defined
                    dqsrm=qsrt/vcell
                    dqsr=dqsrm/rhoqsrp(mv)
                    dqsrf=dqsr
                  end if
                end if
                if(nmat.eq.2 .and. ifqsrp(mv).eq.2) dqsrf=zero
            else if(qsrt.lt.-ztest) then ! for sink
                if(ivqsrp(mv).eq.1) then ! volume flow rate is defined
                  dqsr=qsrt/vcell
                   if(ifrho.ge.2) then
                    dqsrm=dqsr*rho(ijk)
                  else
                    dqsrm=dqsr*rhof
                  endif
                  dqsrf=dqsr
                else if(ivqsrp(mv).eq.2 .or. ivqsrp(mv).eq.3) then ! mass flow rate is defined
                  dqsrm=qsrt/vcell
                  if(nmat.eq.1 .and. ifsnkp(mv).eq.2) then  ! sink of void for one-fluid problems
                    dqsr=dqsrm/rhosnkp(mv)
                  else if(nmat.eq.1 .and. ifsnkp(mv).eq.3) then
                    rhomix=f(ijk)*rho(ijk)+(one-f(ijk))* & 
     &                     rhosnkp(mv)
                    dqsr=dqsrm/rhomix
                  else
                    dqsr=dqsrm/rho(ijk)
                  end if
                  dqsrf=dqsr
                end if
!
                if(nmat.eq.1) then  ! one fluid
                  if(ifvof.gt.0 .and. giswat(nbl).eq.0) then
                    if(ifsnkp(mv).eq.1) then
                      dqsrf=f(ijk)*dqsrf
                      if(dcor.lt.half) dqsrm=f(ijk)*dqsrm
                    else if(ifsnkp(mv).eq.2) then
                      dqsr=(one-f(ijk))*dqsr
                      dqsrm=zero
                      dqsrf=zero
                    else if(ifsnkp(mv).eq.3) then  ! mixture, need test!!!
                      dqsrm=f(ijk)*dqsrm
                      dqsrf=f(ijk)*dqsrf
                    end if
                  end if
!
                else  ! two fluids
                  if(icmprs.eq.0) then
                    if(ifsnkp(mv).eq.3) then
                      dqsrf=f(ijk)*dqsrf
                    else if(ifsnkp(mv).eq.1) then
                      dqsrm=f(ijk)*dqsrm ! fluid mass change rate for fluid #1
                      dqsr=dqsrm/rhof   ! fluid volume change rate for fluid #1
                      dqsrf=dqsr
                    else if(ifsnkp(mv).eq.2) then
                      dqsrm=(one-f(ijk))*dqsrm   ! fluid mass change rate for fluid #2
                      dqsr=dqsrm/rhof2           ! fluid volume change rate for fluid #2
                      dqsrf=zero
                    end if
                  else if(icmprs.eq.1) then
                    if(ivqsrp(mv).eq.2 .or. ivqsrp(mv).eq.3) then  ! mass flow rate is defined
                      if(ifsnkp(mv).eq.3) then
                        dqsrf=f(ijk)*dqsrf
                      else if(ifsnkp(mv).eq.1) then
                        dqsrm=f(ijk)*dqsrm   ! fluid mass change rate for fluid #1
                        dqsr=dqsrm/rhof   ! fluid volume change rate for fluid #1
                        dqsrf=dqsr
                      else if(ifsnkp(mv).eq.2) then
                        if(f(ijk).lt.omepsf) then
                          dqsrm=(one-f(ijk))*dqsrm    ! fluid mass change rate for fluid #2
                          rho2=(rho(ijk)-f(ijk)*rhof)/(one-f(ijk))
                          dqsr=dqsrm/rho2  ! fluid volume change rate for fluid #2
                        end if
                        dqsrf=zero
                      end if
                    else if(ivqsrp(mv).eq.1) then  ! volume flow rate is defined
                      if(ifsnkp(mv).eq.3) then
                        dqsrf=f(ijk)*dqsrf
                      else if(ifsnkp(mv).eq.1) then
                        dqsr=f(ijk)*dqsr   ! fluid mass change rate for fluid #1
                        dqsrm=rhof*dqsr    ! fluid volume change rate for fluid #1
                        dqsrf=dqsr
                      else if(ifsnkp(mv).eq.2) then
                        dqsrm=(rho(ijk)-f(ijk)*rhof)*dqsr !mass change rate for fluid #2
                        dqsr=(one-f(ijk))*dqsr
                       dqsrf=zero
                      end if
                    end if
                  end if
                end if
            end if
!
! NOTE (MRB, 3-24-2014):
!
! F_source=1/wv_ratio is the fluid fraction in the fluid emitted at the source.
!
! If source emits dispersed fluid #1, i.e., F_source < 1, then
! assume that it instantly compacts to F of fluid already in cell. The bulk
! volume source term QSR is then scaled down accordingly.
! (see also PSZOF in PARVEL).
!
            if(nmat.eq.1 .and. f(ijk).gt.emf .and. nf(ijk).eq.0) then
              f_source=one/wv_ratio_t(mv)
              f_weight=min(one,f_source/f(ijk))
            else
              f_weight=one
            endif
!
            qsr(ijk)=qsr(ijk)+dqsr*wv_ratio_t(mv)*f_weight
            qsrm(ijk)=qsrm(ijk)+dqsrm
            if(ifvis.gt.0) then
              tke(ijk)= tkemn
              if(ifvis.lt.5) then
                dtke(ijk)=max(ztest,cnutq*tke(ijk)**x1p5/tlen)
                if(ifvis.ge.3) then
                  tlength(ijk)=cnutq*tke(ijk)**x1p5/dtke(ijk)
                  ttime(ijk)=tke(ijk)/dtke(ijk)
                endif
              endif
              tken(ijk)=tke(ijk)
              dtken(ijk)=dtke(ijk)
            endif
            if(ifvof.gt.0) qsrf(ijk)=qsrf(ijk)+dqsrf
!
            if(ifenrg.ge.2) then
              dqsrme=zero
              if(dqsr.lt.zero) then !sink
                temp=tn(ijk)
              else  ! source
                temp=tqsrp(mv)
              end if
              if(ifqsrsnk.eq.3) then  ! sink is for mixture
                dqsrme=dqsr*rhoe(ijk)
                if(icmprs.eq.1) & 
     &            dqsrme=dqsrme+(one-f(ijk))*dqsr*p(ijk)  ! for compressible flow
                if(decor.gt.half) dqsrme=dqsrme- & 
     &               dqsrm*e1cal(p(ijk),tn(ijk),zero,-one)
              else if(ifqsrsnk.eq.1) then  ! source or sink is for fluid #1
                dqsrme=dqsrm*e1cal(p(ijk),temp,zero,-one)
                if(nmat.eq.1 .and. ifvof.gt.0 .and. giswat(nbl).eq.0 & 
     &               .and. qsrt.lt.-ztest) dqsrme=dqsrme*f(ijk)
                if(icmprs.eq.1) then
                  dqsr0=qsrt/(vcell*rhof)
                  dqsrme=dqsrme+(one-f(ijk))*dqsr0*p(ijk)
                end if
                if(decor.gt.half) dqsrme=dqsrme- & 
     &               dqsrm*e1cal(p(ijk),tn(ijk),zero,-one)
              else if(ifqsrsnk.eq.2 .and. nmat.eq.2) then  ! source or sink is for fluid #2
                if(icmprs.ne.0) then
                  if(f(ijk).lt.omepsf .and. incg.gt.0) then
                    rf2tt=rf2t(ijk)
                    e2=e2cal_nc(p(ijk),temp,sclr(ijk,incg), & 
     &                         rho(ijk),f(ijk))
                  else
                    rf2tt=rf2
                    e2=e2cal(p(ijk),temp)
                  end if
                else
                  rf2tt=zero
                  e2=e2cal(p(ijk),temp)
                end if
                dqsrme=dqsrm*(e2+rf2tt*temp) ! cp is used instead cv for compressible fluid
              end if
              qsrme(ijk)=qsrme(ijk)+dqsrme
            end if
!
! diagnostics
!
            if(icstat(ijk_str).gt.0) then
                delf=dqsr*vcell*delt
                if(delf.lt.-ztest .and. nmat.eq.1) then
                  delft=delf
                  delf=max(-f(ijk)*vcell,delf)
                  dqsrm=dqsrm*delf/delft
                endif
                smass=smass+dqsrm*vcell*delt
                if(ifqsrsnk.eq.1) then
                   sf1=sf1+delf
#ifdef MPI
                   vl1_out_mpi=vl1_out_mpi-min(delf,zero)
                   vl1_in_mpi=vl1_in_mpi+max(delf,zero)
#else
                   vl1_out=vl1_out-min(delf,zero)
                   vl1_in=vl1_in+max(delf,zero)
#endif
                else if(ifqsrsnk.eq.2) then
                   sf2=sf2+delf
#ifdef MPI
                   vl2_out_mpi=vl2_out_mpi-min(delf,zero)
#else
                   vl2_out=vl2_out-min(delf,zero)
#endif
                else if(ifqsrsnk.eq.3.and.dqsr.lt.zero) then
                   sf1=sf1+delf*f(ijk)
                   sf2=sf2+delf*(one-f(ijk))
                end if
            end if
!
            if((ifqsrsnk.eq.2 .or. ifqsrsnk.eq.3) & 
     &                          .and. nmat.eq.1) then
               if(nf(ijk).gt.7) then
#ifdef MPI
                 dvrp_mpi(nf(ijk))=dvrp_mpi(nf(ijk))+gamma*pr(nf(ijk))* & 
     &                         dqsr*vcell*delt
#else
                 dvrp(nf(ijk))=dvrp(nf(ijk))+gamma*pr(nf(ijk))* & 
     &                         dqsr*vcell*delt
#endif
               elseif(nf(ijk).ge.1 .and. nf(ijk).le.6) then
                 call nvcell(ijk,nvoid,istat)
                 if(nvoid.gt.7) then
! sink or source of pure gas (void)
                   if((qsrt.lt.-ztest.and.ifsnkp(mv).eq.2) .or. & 
     &                        qsrt.gt.ztest) then
#ifdef MPI
                     dvrp_mpi(nvoid)=dvrp_mpi(nvoid)+gamma*pr(nvoid)* & 
     &                 dqsr*vcell*delt
#else
                     dvrp(nvoid)=dvrp(nvoid)+gamma*pr(nvoid)*dqsr* & 
     &                 vcell*delt
#endif
                   else
! sink of gas/liquid mixture
#ifdef MPI
                     dvrp_mpi(nvoid)=dvrp_mpi(nvoid)+gamma*pr(nvoid)* & 
     &                           dqsr*vcell*delt*(one-f(ijk))
#else
                     dvrp(nvoid)=dvrp(nvoid)+gamma*pr(nvoid)*dqsr* & 
     &                           vcell*delt*(one-f(ijk))
#endif
                   end if
                 end if
               end if
            end if
!
! scalar source
!
            delf1=dqsrf*delt        ! volume change for fluid #1 scalars
            if(ifqsrsnk.eq.2) then
              delf2=dqsr*delt     ! volume change for fluid #2 scalars, pure fluid #2 sink or source
            else
              delf2=(one-f(ijk))*dqsr*delt  ! volume change for fluid #2 scalars, mixture sink
            endif
            delf3=dqsr*delt
!
            do ns=1,nsc
              if(iscpqsrc(mv,ns).gt.1) then
! external table
                call readtimedata(ifile_unit_par(nvar_par+ns,mv), & 
     &                            fscpqsr(mv,ns), & 
     &                            if_open_par(nvar_par+ns,mv), & 
     &                            num_vals_par(nvar_par+ns,mv), & 
     &                            iposition_par(nvar_par+ns,mv), & 
     &                            times_par(:,nvar_par+ns,mv), & 
     &                            values_par(:,nvar_par+ns,mv), & 
     &                            iscpqsrc(mv,ns),scqsr,0)
!
! convert NC gas volume fraction at boundaries (from external table) to density
!
                if(incg.gt.0.and.ns.eq.incg) then
                  if(ifqsrp(mv).eq.2) then
                    if(iscpqsrc(mv,incg).gt.1) then
                      if(irsat.eq.0) then
                        vfnc1=max(zero,min(scqsr,one))
                        scqsr=vfnc1*rhoqsrp(mv)*rf2/ & 
     &                  (rfnc-vfnc1*(rfnc-rf2))
                      else
                        psattmp=psat(tqsrp(mv),1)
                        bval=-rhoqsrp(mv)*rfnc*tqsrp(mv)
                        cval=rhoqsrp(mv)*scqsr*psattmp*tqsrp(mv) & 
     &                                                  *(rfnc-rf2)
                        discrim=max(zero,bval**2-four*cval)
                        ptot=(-bval+sqrt(discrim))/two
                        vfvap=scqsr*psattmp/ptot
                        vfnc1=max(zero,min(one-vfvap,one))
                        scqsr=vfnc1*rhoqsrp(mv)*rf2/ & 
     &                      (rfnc-vfnc1*(rfnc-rf2))
                      endif
                    endif
                  endif
                endif
              else
! table in prepin
                scqsr=alinint(tref,tms,scpqsr(1,mv,ns),ntms)
              endif
!
              ff=f(ijk)
              if((ifqsrsnk.eq.1 .or. ifqsrsnk.eq.3) & 
     &                      .and. isclr(ns).eq.3) then
                if(ff.gt.emf) then
                  if(delf1.gt.zero) then
                    sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                       scqsr*delf1/ff
                    if(dscorsc(ns).gt.half) sclr(ijk,ns)=sclr(ijk,ns)- & 
     &                       sclrn(ijk,ns)*delf1/ff
                  else
                    if(dscorsc(ns).lt.half) sclr(ijk,ns)=sclr(ijk,ns)+   & ! do not apply sink source f
     &                         sclrn(ijk,ns)*delf1/ff
                  endif
                else
! special logic for newly filled cells: goes together with the change in vfsources
                  if(delf1.gt.zero) sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                        scqsr*delf1
                endif
              elseif((ifqsrsnk.eq.2 .or. ifqsrsnk.eq.3) & 
     &                                   .and.isclr(ns).eq.5) then
!
                if(ff.lt.emf1) then
                  if(delf2.gt.zero) then
                    sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                        scqsr*delf2/(one-ff)
                    if(dscorsc(ns).gt.half) sclr(ijk,ns)=sclr(ijk,ns)- & 
     &                        sclrn(ijk,ns)*delf2/(one-ff)
                  else
                    if(dscorsc(ns).lt.half) sclr(ijk,ns)=sclr(ijk,ns)+   & ! do not apply sink source f
     &                          sclrn(ijk,ns)*delf2/(one-ff)
                  endif
                else
! special logic for newly filled cells: goes together with the change in vfsources
                  if(delf2.gt.zero) sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                    scqsr*delf2
                endif
              elseif(isclr(ns).eq.1) then
                if(delf3.gt.zero) then
                  sclr(ijk,ns)=sclr(ijk,ns)+scqsr*delf3
                  if(dscorsc(ns).gt.half) sclr(ijk,ns)=sclr(ijk,ns)- & 
     &                         sclrn(ijk,ns)*delf3
                else
                  if(dscorsc(ns).lt.half) sclr(ijk,ns)=sclr(ijk,ns)+   & ! do not apply sink source for
     &                       sclrn(ijk,ns)*delf3
                endif
              endif
            enddo
!
! Apportion (1-vf) part of the source in a neighbor, if ijk cell is partially blocked
!
            ijkx=0
!
            if(vf(ijk).lt.omeps) then
              if(cyl.gt.half) then
! Convert Cartesian source velocities to cylindrical velocities.
! This is necessary because they are always Cartesian.
                  call vel_trans(1,nbl,uv_t(mv),vv_t(mv), & 
     &                source_particle(k1)%x_part, & 
     &                source_particle(k1)%y_part,uflo,vflo)
              else
                  uflo=uv_t(mv)
                  vflo=vv_t(mv)
              endif
              wflo=wv_t(mv)
!
              include '../comdeck/mijk.f90'
              include '../comdeck/pijk.f90'
!
              if(abs(uflo).ge.abs(vflo) .and. & 
     &           abs(uflo).ge.abs(wflo)) then
                if(uflo.gt.zero .and. afr(ijk).gt.em6) then
                  ijkx=ipjk
                elseif(afr(imjk).gt.em6) then
                  ijkx=imjk
                else
                  ijkx=ijk
                endif
              elseif(abs(wflo).ge.abs(vflo)) then
                if(wflo.gt.zero .and. aft(ijk).gt.em6) then
                  ijkx=ijkp
                elseif(aft(ijkm).gt.em6) then
                  ijkx=ijkm
                else
                  ijkx=ijk
                endif
              else
                if(vflo.gt.zero .and. afb(ijk).gt.em6) then
                  ijkx=ijpk
                elseif(afb(ijmk).gt.em6) then
                  ijkx=ijmk
                else
                  ijkx=ijk
                endif
              endif
              if(vf(ijkx).lt.em6) ijkx=ijk
              ijkx_str=ijk_unstr2str(ijkx)
              ix=i_str(ijkx_str)
              jx=j_str(ijkx_str)
              kx=k_str(ijkx_str)
              vcellx=vf(ijkx)*delx(ix)*dely(jx)/rri(ix)*delz(kx)
              qsrt=source_particle(k1)%flow_rate*(one-vf(ijk))
!
!  calculate qsr(ijkx), qsrm(ijkx), qsrf(ijkx)
!
              dqsr=zero
              dqsrm=zero
              dqsrf=zero
              if(qsrt.gt.ztest) then ! for source
                if(ivqsrp(mv).eq.1) then ! volume flow rate is defined
                  dqsr=qsrt/vcellx
                  dqsrm=dqsr*rhoqsrp(mv)
                  dqsrf=dqsr
                else if(ivqsrp(mv).eq.2 .or. ivqsrp(mv).eq.3) then ! mass flow rate is defined
                  dqsrm=qsrt/vcellx
                  dqsr=dqsrm/rhoqsrp(mv)
                  dqsrf=dqsr
                end if
                if(nmat.eq.1 .and. ifqsrp(mv).eq.2) then ! if the source release bubble
                  dqsrm=zero
                  dqsrf=zero
                end if
                if(nmat.eq.2 .and. ifqsrp(mv).eq.2) dqsrf=zero
!
              else if(qsrt.lt.-ztest) then ! for sink
                if(ivqsrp(mv).eq.1) then ! volume flow rate is defined
                  dqsr=qsrt/vcellx
                  if(ifrho.ge.2) then
                    dqsrm=dqsr*rho(ijkx)
                  else
                    dqsrm=dqsr*rhof
                  endif
                  dqsrf=dqsr
                else if(ivqsrp(mv).eq.2 .or. ivqsrp(mv).eq.3) then ! mass flow rate is defined
                  dqsrm=qsrt/vcellx
                  if(nmat.eq.1 .and. ifsnkp(mv).eq.2) then  ! sink of void for one-fluid problems
                    dqsr=dqsrm/rhosnkp(mv)
                  else if(nmat.eq.1 .and. ifsnkp(mv).eq.3) then
                    rhomix=f(ijkx)*rho(ijkx)+(one-f(ijkx))* & 
     &                     rhosnkp(mv)
                    dqsr=dqsrm/rhomix
                  else
                    dqsr=dqsrm/rho(ijkx)
                  end if
                  dqsrf=dqsr
                end if
!
                if(nmat.eq.1) then  ! one fluid
                  if(ifvof.gt.0 .and. giswat(nbl).eq.0) then
                    if(ifsnkp(mv).eq.1) then
                      dqsrf=f(ijkx)*dqsrf
                      if(dcor.lt.half) dqsrm=f(ijkx)*dqsrm
                    else if(ifsnkp(mv).eq.2) then
                      dqsr=(one-f(ijkx))*dqsr
                      dqsrm=zero
                      dqsrf=zero
                    else if(ifsnkp(mv).eq.3) then  ! mixture, need test!!!
                      dqsrm=f(ijkx)*dqsrm
                      dqsrf=f(ijkx)*dqsrf
                    end if
                  end if
!
                else  ! two fluids
                  if(icmprs.eq.0) then
                    if(ifsnkp(mv).eq.3) then
                      dqsrf=f(ijkx)*dqsrf
                    else if(ifsnkp(mv).eq.1) then
                      dqsrm=f(ijkx)*dqsrm ! fluid mass change rate for fluid #1
                      dqsr=dqsrm/rhof   ! fluid volume change rate for fluid #1
                      dqsrf=dqsr
                    else if(ifsnkp(mv).eq.2) then
                      dqsrm=(one-f(ijkx))*dqsrm   ! fluid mass change rate for fluid #2
                      dqsr=dqsrm/rhof2            ! fluid volume change rate for fluid #2
                      dqsrf=zero
                    end if
                  else if(icmprs.eq.1) then
                    if(ivqsrp(mv).eq.2 .or. ivqsrp(mv).eq.3) then  ! mass flow rate is defined
                      if(ifsnkp(mv).eq.3) then
                        dqsrf=f(ijkx)*dqsrf
                      else if(ifsnkp(mv).eq.1) then
                        dqsrm=f(ijkx)*dqsrm   ! fluid mass change rate for fluid #1
                        dqsr=dqsrm/rhof   ! fluid volume change rate for fluid #1
                        dqsrf=dqsr
                      else if(ifsnkp(mv).eq.2) then
                        if(f(ijkx).lt.omepsf) then
                          dqsrm=(one-f(ijkx))*dqsrm    ! fluid mass change rate for fluid #2
                          rho2=(rho(ijkx)-f(ijkx)*rhof)/(one-f(ijkx))
                          dqsr=dqsrm/rho2  ! fluid volume change rate for fluid #2
                        end if
                        dqsrf=zero
                      end if
                    else if(ivqsrp(mv).eq.1) then  ! volume flow rate is defined
                      if(ifsnkp(mv).eq.3) then
                        dqsrf=f(ijkx)*dqsrf
                      else if(ifsnkp(mv).eq.1) then
                        dqsr=f(ijkx)*dqsr   ! fluid mass change rate for fluid #1
                        dqsrm=rhof*dqsr    ! fluid volume change rate for fluid #1
                        dqsrf=dqsr
                      else if(ifsnkp(mv).eq.2) then
                        dqsrm=(rho(ijkx)-f(ijkx)*rhof)*dqsr
                        dqsr=(one-f(ijkx))*dqsr   ! fluid mass change rate for fluid #2
                        dqsrf=zero
                      end if
                    end if
                  end if
                end if
              end if
!
! NOTE (MRB, 3-24-2014):
!
! F_source=1/wv_ratio is the fluid fraction in the fluid emitted at the source.
!
! If source emits dispersed fluid #1, i.e., F_source < 1, then
! assume that it instantly compacts to F of fluid already in cell. The bulk
! volume source term QSR is then scaled down accordingly.
! (see also PSZOF in PARVEL).
!
              if(nmat.eq.1 .and. f(ijkx).gt.emf .and. nf(ijkx).eq.0) & 
     &                                   then
                f_source=one/wv_ratio_t(mv)
                f_weight=min(one,f_source/f(ijkx))
              else
                f_weight=one
              endif
!
              qsr(ijkx)=qsr(ijkx)+dqsr*wv_ratio_t(mv)*f_weight
              qsrm(ijkx)=qsrm(ijkx)+dqsrm
              if(ifvis.gt.0) then
                tke(ijk)= tkemn
                if(ifvis.lt.5) then
                  dtke(ijk)=max(ztest,cnutq*tke(ijk)**x1p5/tlen)
                  if(ifvis.ge.3) then
                    tlength(ijk)=cnutq*tke(ijk)**x1p5/dtke(ijk)
                    ttime(ijk)=tke(ijk)/dtke(ijk)
                  endif
                endif
                tken(ijk)=tke(ijk)
                dtken(ijk)=dtke(ijk)
              endif
              if(ifvof.gt.0) qsrf(ijkx)=qsrf(ijkx)+dqsrf
!
              if(ifenrg.ge.2) then
                dqsrme=zero
                if(dqsr.lt.zero) then !sink
                  temp=tn(ijkx)
                else  ! source
                  temp=tqsrp(mv)
                end if
                if(ifqsrsnk.eq.3) then  ! sink is for mixture
                  dqsrme=dqsr*rhoe(ijkx)
                  if(icmprs.eq.1) & 
     &              dqsrme=dqsrme+(one-f(ijkx))*dqsr*p(ijkx)  ! for compressible flow
                  if(decor.gt.half) dqsrme=dqsrme- & 
     &               dqsrm*e1cal(p(ijkx),tn(ijkx),zero,-one)
                else if(ifqsrsnk.eq.1) then  ! source or sink is for fluid #1
                  dqsrme=dqsrm*e1cal(p(ijkx),temp,zero,-one)
                  if(nmat.eq.1 .and. ifvof.gt.0 .and. giswat(nbl).eq.0 & 
     &                   .and. qsrt.lt.-ztest) dqsrme=dqsrme*f(ijkx)
                  if(icmprs.eq.1) then
                    dqsr0=qsrt/(vcellx*rhof)
                    dqsrme=dqsrme+(one-f(ijkx))*dqsr0*p(ijkx)
                  end if
                  if(decor.gt.half) dqsrme=dqsrme- & 
     &               dqsrm*e1cal(p(ijkx),tn(ijkx),zero,-one)
                else if(ifqsrsnk.eq.2 .and. nmat.eq.2) then  ! source or sink is for fluid #2
                  if(icmprs.ne.0) then
                    if(f(ijkx).lt.omepsf .and. incg.gt.0) then
                      rf2tt=rf2t(ijkx)
                      e2=e2cal_nc(p(ijkx),temp,sclr(ijkx,incg), & 
     &                                rho(ijkx),f(ijkx))
                    else
                      rf2tt=rf2
                      e2=e2cal(p(ijkx),temp)
                    end if
                  else
                    rf2tt=zero
                    e2=e2cal(p(ijk),temp)
                  end if
                  dqsrme=dqsrm*(e2+rf2tt*temp) ! cp is used instead cv for compressible fluid
                end if
                qsrme(ijkx)=qsrme(ijkx)+dqsrme
              end if
!
!             diagnositcs
!
              ijkx_str=ijk_unstr2str(ijkx)
              if(icstat(ijkx_str).gt.0) then
                   delf=dqsr*vcellx*delt
                   if(delf.lt.-ztest .and. nmat.eq.1) then
                     delft=delf
                     delf=max(-f(ijk)*vcell,delf)
                     dqsrm=dqsrm*delf/delft
                   endif
                   smass=smass+dqsrm*vcellx*delt
                   if(ifqsrsnk.eq.1) then
                     sf1=sf1+delf
#ifdef MPI
                     vl1_out_mpi=vl1_out_mpi-min(delf,zero)
                     vl1_in_mpi=vl1_in_mpi+max(delf,zero)
#else
                     vl1_out=vl1_out-min(delf,zero)
                     vl1_in=vl1_in+max(delf,zero)
#endif
                   else if(ifqsrsnk.eq.2) then
                     sf2=sf2+delf
#ifdef MPI
                     vl2_out_mpi=vl2_out_mpi-min(delf,zero)
#else
                     vl2_out=vl2_out-min(delf,zero)
#endif
                   else if(ifqsrsnk.eq.3.and.dqsr.lt.zero) then
                     sf1=sf1+delf*f(ijkx)
                     sf2=sf2+delf*(one-f(ijkx))
                   end if
              end if
!
              if((ifqsrsnk.eq.2 .or. ifqsrsnk.eq.3) & 
     &                          .and. nmat.eq.1) then
                if(nf(ijkx).gt.7) then
#ifdef MPI
                  dvrp_mpi(nf(ijkx))=dvrp_mpi(nf(ijkx))+gamma* & 
     &                         pr(nf(ijkx))*dqsr*vcellx*delt
#else
                  dvrp(nf(ijkx))=dvrp(nf(ijkx))+gamma*pr(nf(ijkx))* & 
     &                         dqsr*vcellx*delt
#endif
                elseif(nf(ijkx).ge.1 .and. nf(ijkx).le.6) then
                  call nvcell(ijkx,nvoid,istat)
                  if(nvoid.gt.7) then
! sink or source of pure gas (void)
                    if((qsrt.lt.-ztest.and.ifsnkp(mv).eq.2) .or. & 
     &                        qsrt.gt.ztest) then
#ifdef MPI
                      dvrp_mpi(nvoid)=dvrp_mpi(nvoid)+gamma*pr(nvoid)* & 
     &                 dqsr*vcellx*delt
#else
                      dvrp(nvoid)=dvrp(nvoid)+gamma*pr(nvoid)*dqsr* & 
     &                 vcellx*delt
#endif
                    else
! sink of gas/liquid mixture
#ifdef MPI
                      dvrp_mpi(nvoid)=dvrp_mpi(nvoid)+gamma*pr(nvoid)* & 
     &                           dqsr*vcellx*delt*(one-f(ijkx))
#else
                      dvrp(nvoid)=dvrp(nvoid)+gamma*pr(nvoid)*dqsr* & 
     &                           vcellx*delt*(one-f(ijkx))
#endif
                    end if
                  end if
                end if
              end if
!
            delf1=dqsrf*delt        ! volume change for fluid #1 scalars
            if(ifqsrsnk.eq.2) then
              delf2=dqsr*delt     ! volume change for fluid #2 scalars, pure fluid #2 sink or source
            else
              delf2=(one-f(ijkx))*dqsr*delt  ! volume change for fluid #2 scalars, mixture sink
            endif
            delf3=dqsr*delt
!
            do ns=1,nsc
              if(iscpqsrc(mv,ns).gt.1) then
! external table
                call readtimedata(ifile_unit_par(nvar_par+ns,mv), & 
     &                            fscpqsr(mv,ns), & 
     &                            if_open_par(nvar_par+ns,mv), & 
     &                            num_vals_par(nvar_par+ns,mv), & 
     &                            iposition_par(nvar_par+ns,mv), & 
     &                            times_par(:,nvar_par+ns,mv), & 
     &                            values_par(:,nvar_par+ns,mv), & 
     &                            iscpqsrc(mv,ns),scqsr,0)
!
! convert NC gas volume fraction at boundaries (from external table) to density
!
                if(incg.gt.0.and.ns.eq.incg) then
                  if(ifqsrp(mv).eq.2) then
                    if(iscpqsrc(mv,incg).gt.1) then
                      if(irsat.eq.0) then
                        vfnc1=max(zero,min(scqsr,one))
                        scqsr=vfnc1*rhoqsrp(mv)*rf2/ & 
     &                  (rfnc-vfnc1*(rfnc-rf2))
                      else
                        psattmp=psat(tqsrp(mv),1)
                        bval=-rhoqsrp(mv)*rfnc*tqsrp(mv)
                        cval=rhoqsrp(mv)*scqsr*psattmp*tqsrp(mv) & 
     &                                                  *(rfnc-rf2)
                        discrim=max(zero,bval**2-four*cval)
                        ptot=(-bval+sqrt(discrim))/two
                        vfvap=scqsr*psattmp/ptot
                        vfnc1=max(zero,min(one-vfvap,one))
                        scqsr=vfnc1*rhoqsrp(mv)*rf2/ & 
     &                      (rfnc-vfnc1*(rfnc-rf2))
                      endif
                    endif
                  endif
                endif
              else
! table in prepin
                scqsr=alinint(tref,tms,scpqsr(1,mv,ns),ntms)
              endif
              ff=f(ijkx)
              if((ifqsrsnk.eq.1 .or. ifqsrsnk.eq.3) & 
     &                      .and. isclr(ns).eq.3) then
                if(ff.gt.emf) then
                  if(delf1.gt.zero) then
                    sclr(ijkx,ns)=sclr(ijkx,ns)+ & 
     &                       scqsr*delf1/ff
                    if(dscorsc(ns).gt.half) sclr(ijkx,ns)=sclr(ijkx,ns)- & 
     &                       sclrn(ijkx,ns)*delf1/ff
                  else
                    if(dscorsc(ns).lt.half) sclr(ijkx,ns)=sclr(ijkx,ns)+   & ! do not apply sink source
     &                         sclrn(ijkx,ns)*delf1/ff
                  endif
                else
! special logic for newly filled cells: goes together with the change in vfsources
                  if(delf1.gt.zero) sclr(ijkx,ns)=sclr(ijkx,ns)+ & 
     &                             scqsr*delf1
                endif
              elseif((ifqsrsnk.eq.2 .or. ifqsrsnk.eq.3) & 
     &                                   .and.isclr(ns).eq.5) then
!
                if(ff.lt.emf1) then
                  if(delf2.gt.zero) then
                    sclr(ijkx,ns)=sclr(ijkx,ns)+ & 
     &                        scqsr*delf2/(one-ff)
                    if(dscorsc(ns).gt.half) sclr(ijkx,ns)=sclr(ijkx,ns)- & 
     &                        sclrn(ijkx,ns)*delf2/(one-ff)
                  else
                    if(dscorsc(ns).lt.half) sclr(ijkx,ns)=sclr(ijkx,ns)+   & ! do not apply sink source
     &                          sclrn(ijkx,ns)*delf2/(one-ff)
                  endif
                else
! special logic for newly filled cells: goes together with the change in vfsources
                  if(delf2.gt.zero) sclr(ijkx,ns)=sclr(ijkx,ns)+ & 
     &                    scqsr*delf2
                endif
              elseif(isclr(ns).eq.1) then
                if(delf3.gt.zero) then
                  sclr(ijkx,ns)=sclr(ijkx,ns)+scqsr*delf3
                  if(dscorsc(ns).gt.half) sclr(ijkx,ns)=sclr(ijkx,ns)- & 
     &                         sclrn(ijkx,ns)*delf3
                else
                  if(dscorsc(ns).lt.half) sclr(ijkx,ns)=sclr(ijkx,ns)+   & ! do not apply sink source f
     &                       sclrn(ijkx,ns)*delf3
                endif
              endif
            enddo
!
            end if
!
 3200   continue
!
      endif
#ifdef MPI
      nvoids=nreg+7
      if(nmybl.eq.myblcks .and. nvoids.gt.7) then
        sendbuf(8:nreg+7)=dvrp_mpi(8:nreg+7)
        call parallel_sum(sendbuf(8:nreg+7),dvrp_mpi)
        do l=8,nreg+7
          dvrp(l)=dvrp(l)+dvrp_mpi(l)
        enddo
      endif
#endif
!
! ------ Add Valve/Vents for Fluid #2 to 2-Fluid Cases with Compressible Fluid 2
!
      if(icmprs.gt.0) then
!
        do 4000 l=1,nvals
!
! active simulation control of valve events
! check if the valve is turned off by events
! since this subroutine is called after CICSTP, where the valve event status has already been checked
! thus only need to check the OFF status here
         if(ievent_val_off(l).gt.0) cycle
!
          do 4010 nn=ijkval(1,l),ijkval(2,l)
            ijk=ijkvnt(nn)
#ifdef MPI
            call ijknbl(ijk,nnbl)
            if(rankmpi.ne.block_owner(nnbl)) cycle
#endif
            if(ijk.eq.0) cycle
            if(vf(ijk).lt.emf) cycle
            if(f(ijk).ge.valf_max(l)) cycle
!
            pdiff=p(ijk)-valp_cnt(l)
! cycle if outflow is not allowed
            if(pdiff.gt. ztest .and. if_2way(l).eq.2) cycle
! cycle if inflow is not allowed
            if(pdiff.lt.-ztest .and. if_2way(l).eq.1) cycle
!
            apdiff=abs(pdiff)
            if(apdiff.lt.ztest) cycle
!
! * * * VOF source term
!
            dqsrf=zero
!
            include '../comdeck/ijk2str.f90'
            i=i_str(ijk_str)
            j=j_str(ijk_str)
            k=k_str(ijk_str)
!
            vcell=vf(ijk)*delx(i)*dely(j)/rri(i)*delz(k)
!
! * * * gas velocity in vent based on the loss coefficient and pressure differential
            vent_vel_valc=valc(l)*sqrt(apdiff)/valarea(l)
!
! * * * gas density at source
            if(incg.gt.0) then
                rho2=rho2cl_nc(p(ijk),tn(ijk),sclr(ijk,incg))
            else
                rho2=rho2cl(p(ijk),tn(ijk),rf2)
            endif
!
            if(pdiff.gt.zero) then
!
! * * * gas is leaving through the vent
!
              if(incg.gt.0) then
! gamm2=Cp/Cv of gas in cell
                  gamm2=one+p(ijk)/(rhocv(ijk)*tn(ijk))
! gas constant in cell
                  rf_gas=rf2t(ijk)
              else
! gamm2=Cp/Cv of gas in cell
                  gamm2=one+rf2/cv2
! gas constant in cell
                  rf_gas=rf2
              endif
!
              if(gamm2.le.one .or. rf_gas.lt.ztest) cycle
 
! * * * gas velocity in cell
              include '../comdeck/mijk.f90'
! use max value for each component since flow may be in opposite directions at opposing faces
              u_velocity=max(abs(u(imjk)),abs(u(ijk)))
              v_velocity=max(abs(v(ijmk)),abs(v(ijk)))
              w_velocity=max(abs(w(ijkm)),abs(w(ijk)))
              vel_magnitude2=u_velocity*u_velocity+ & 
     &                       v_velocity*v_velocity+ & 
     &                       w_velocity*w_velocity
! * * * gas critical velocity squared
              vel_critical2=(gamm2-one)/(gamm2+one)* & 
     &         (vel_magnitude2+(two/(gamm2-one))*gamm2*rf_gas*tn(ijk))
!
              if(vel_critical2.lt.ztest) cycle
!
              vel_ratio2=vel_magnitude2/vel_critical2
!
! stagnation temperature, pressure and density, from Bernoulli
!
              temp_stag=tn(ijk)/ & 
     &  (one-((gamm2-one)/(gamm2+one))*vel_ratio2)
              press_stag=p(ijk)/ & 
     & ((one-((gamm2-one)/(gamm2+one))*vel_ratio2)**(gamm2/(gamm2-one)))
              rho_stag=press_stag/(rf_gas*temp_stag)
!
! Assume pressure at vent exit is the same as the external (L.I.Sedov, Vol II, page 49)
!
              press_ratio=valp_cnt(l)/press_stag
!
            else
!
! * * * gas is entering through the vent
!
              if(valt_cnt(l).lt.ztest) cycle
!
              if(incg.gt.0) then
! gamm2=Cp/Cv of gas in cell
                  gamm2=one+rfnc/cvnc
! gas constant in cell
                  rf_gas=rfnc
              else
! gamm2=Cp/Cv of gas in cell
                  gamm2=one+rf2/cv2
! gas constant in cell
                  rf_gas=rf2
              endif
!
              if(gamm2.le.one .or. rf_gas.lt.ztest) cycle
! * * * gas critical velocity squared
              vel_critical2=two/(gamm2+one)*gamm2*rf_gas*valt_cnt(l)
!
              if(vel_critical2.lt.ztest) cycle
!
! assume external gas pressure, temperature and density are at stagnation
!
              press_stag=valp_cnt(l)
              rho_stag=press_stag/(rf_gas*valt_cnt(l))
!
! Assume pressure at vent entrance is the same as the cell pressure (L.I.Sedov, Vol II, page 49)
!
              press_ratio=p(ijk)/press_stag
!
            endif
!
! velocity through the vent, from Bernoulli
!
            vent_vel_bernoulli=sqrt((two*gamm2/(gamm2-one))* & 
     &                       (press_stag/rho_stag)* & 
     &             (one-press_ratio**((gamm2-one)/gamm2)))
!
            if(vent_vel_bernoulli.lt.ztest) cycle
!
! gas mass flux through the vent, rho*vel, per Bernoulli
            vent_mass_flux_bernoulli= & 
     &        sqrt((gamm2+one)/(gamm2-one))*rho_stag* & 
     &        (press_ratio**(one/gamm2))*sqrt(vel_critical2)* & 
     &        sqrt(one-press_ratio**((gamm2-one)/gamm2))
! gas density in vent, from Bernoulli
            rho_vent=vent_mass_flux_bernoulli/vent_vel_bernoulli
!
            if(rho_vent.lt.ztest) cycle
!
! gas mass flow rate computed from VALC and pressure differential, vent_vel_valc*rho_vent,
! cannot exceed vent_mass_flux_bernoulli.
            vent_mass_flux_valc=min(vent_vel_valc*rho_vent, & 
     &                                vent_mass_flux_bernoulli)
!
! * * * Mass source
!
! gas mass flow rate is the conserved quantity along a stream line
! scale by 1-f assuming the vent is partially covered with fluid #1
!
! MRB: this mass flow rate, if used without dqsrm2 below, requires a time step
!      size limit for stability. Have not figured out the limit.
!
            dqsrm1=(one-f(ijk))*vent_mass_flux_valc*valarea(l)/vcell
!
! * * * combine with implicit pressure treatment in the cell
!
! ------ Estimate volume inflow across cell boundaries
            include '../comdeck/mijk.f90'
            dvf0=rri(i)*(rdx(i)*(afr(ijk)*u(ijk)/rr(i)- & 
     &                   afr(imjk)*u(imjk)/rr(i-1))+ & 
     &           rdy(j)*(afb(ijk)*v(ijk)-afb(ijmk)*v(ijmk)))+ & 
     &           rdz(k)*(aft(ijk)*w(ijk)-aft(ijkm)*w(ijkm))
!
            if(abs(dvf0).gt.ztest) then
! ------ Compute Implicit Flow Rate for Valve/Vent
              if(pdiff.gt.zero) then
                dvf0=min(dvf0,zero)
              else
                dvf0=max(dvf0,zero)
              endif
              dvf02=-delt*gamm2*p(ijk)*min(abs(dvf0), & 
     &           valc(l)*sqrt(apdiff)/vcell)
              vgold=(one-f(ijk))*vcell
              dpa=vgold/(delt*gamm2*p(ijk)*valc(l))
              dpc=dpa*(apdiff-dvf02)    ! combine out and in volume flows
              dpsol=qsolv(dpa,-dpc)     ! pnew is resulting adiabatic pressure
              pnew=valp_cnt(l)+dpsol**2
            else
              pnew=p(ijk)
            endif
!
! ------ Compute Valve/Vent Source Terms
!
! * * * Volume Source at Valve/Vent
!
            dqsrm2=rho2*valc(l)*sqrt(abs(pnew-valp_cnt(l)))/vcell
!
! * * * combine the two mass flow rates by using the smallest of the two
!
            dqsrm=-sign(one,pdiff)*min(dqsrm1,dqsrm2)
!
! * * * Volume source: convert using cell gas density
!
            dqsr=dqsrm/rho2
!
! * * * Internal energy changes by flow and volume change
!
! These are approximations since it should be using p, T and rho inside the vent,
! which can be different than those outside the vent on either side of it, used
! below, if the flow is sufficiently fast.
!
            if(pdiff.gt.zero) then
              if(incg.gt.0) then
                e2= & 
     &   e2cal_nc(p(ijk),tn(ijk),sclr(ijk,incg),rho(ijk),f(ijk))
              else
                e2=e2cal(p(ijk),tn(ijk))
              endif
              dqsrme=dqsrm*(e2+rf_gas*tn(ijk))
            else
              if(incg.gt.0) then
! As long as the mixture density and scalar are positive and the same, and
! fluid fraction is zero, the e2 energy comes out as that of pure NC gas.
                e2=e2cal_nc(valp_cnt(l),valt_cnt(l),one,one,zero)
              else
                e2=e2cal(valp_cnt(l),valt_cnt(l))
              endif
              dqsrme=dqsrm*(e2+rf_gas*valt_cnt(l))
            endif
!
!       Add fluid #1 to the source
!
            if(pdiff.gt.zero .and. f(ijk).lt.valf_min(l)) then
              dqsr1=dqsr*f(ijk)/(one-f(ijk)) ! assume both fluids are removed at the same rate
              dqsrf=dqsr1
              dqsrm1=dqsr1*rhof
              dqsrme1=dqsrm1*e1cal(p(ijk),tn(ijk),zero,-one)          
            else
              dqsr1=zero
              dqsrf=zero
              dqsrm1=zero
              dqsrme1=zero
            endif
!
! Update Source Arrays
            qsr(ijk)=qsr(ijk)+dqsr+dqsr1
            qsrf(ijk)=qsrf(ijk)+dqsrf            
            qsrm(ijk)=qsrm(ijk)+dqsrm+dqsrm1
            qsrme(ijk)=qsrme(ijk)+dqsrme+dqsrme1
!
! Diagnostics
            if(icstat(ijk_str).gt.0) then
              smass=smass+(dqsrm+dqsrm1)*vcell*delt
              sf1=sf1+dqsr1*vcell*delt
              sf2=sf2+dqsr*vcell*delt  
              vl1_out=vl1_out-min(dqsr1*vcell*delt,zero)              
              vl2_out=vl2_out-min(dqsr*vcell*delt,zero)
              valq(l)=valq(l)-(dqsr+dqsr1)*vcell ! valve volume flow rate
            endif
!
! Scalar Sources
            delf1=dqsrf*delt
            delf2=dqsr*delt
            delf=delf1+delf2
            ff=f(ijk)
!            
            do ns=1,nsc
              if(isclr(ns).eq.1) then
! delf1 cannot be positive here       
                if(ff.gt.emf) then
                  if(dscorsc(ns).lt.half) sclr(ijk,ns)=   & ! do not apply sink source
     &                sclr(ijk,ns)+sclrn(ijk,ns)*delf1/ff     
                endif
              elseif(isclr(ns).eq.5) then
                if(ff.lt.emf1) then
                  if(delf2.gt.zero) then
                    if(dscorsc(ns).gt.half) sclr(ijk,ns)=sclr(ijk,ns)- & 
     &                             sclrn(ijk,ns)*delf2/(one-ff)
                  else
                    if(dscorsc(ns).lt.half) sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                             sclrn(ijk,ns)*delf2/(one-ff)
                  endif  
                endif
              elseif(isclr(ns).eq.1) then
                if(delf.gt.zero) then
                  if(dscorsc(ns).gt.half) sclr(ijk,ns)=sclr(ijk,ns)- & 
     &                           sclrn(ijk,ns)*delf
                else
                  if(dscorsc(ns).lt.half) sclr(ijk,ns)=sclr(ijk,ns)+ & 
     &                           sclrn(ijk,ns)*delf
                endif
              endif
            end do
!
 4010     continue
 4000   continue
!
      endif
!      
! ------ End of Valve/Vents for Fluid #2
!
! ------ add volume change for thermal expansion in ifrho=1 model
!
      if(ifthexp.eq.1 .and. cycle.gt.1) then
!
        density_ratio=rhofs/rhof
!
        if(ifomp.eq.0) then
          ich=1
        else
          ich=nompch
        endif
!
        call load_bal(ich*numthrds,ijklim2(nbl,1),ijklim2(nbl,2))
!
        call et_on(45)
!
!$omp parallel do schedule(static,1) & 
!$omp& private(n,nn,n2,nnn,ijk) & 
!$omp& private(rhoe_sav,f_sav,tnold,rho1,rho2,rho1_2,rho2_2,rcsqfv) & 
!$omp& private(tnnew,denom) & 
!$omp& private(tn_sav,fliqold,qsr_liq_shr,frslv,frslvn,qsr_sld_shr) & 
!$omp& firstprivate(density_ratio)
!
      do 2502 n=1,numthrds
!
        do 2501 n2=1,ich
!
          nnn=(n2-1)*numthrds+n
!
          do 2500 nn=klo(nnn),khi(nnn)
!
            ijk=ijk2m1(nn)
!
            if(vf(ijk).lt.em6) cycle
            if((nmat.eq.1 .or. icmprs.eq.1) .and. & 
     &         (f(ijk).lt.emf .or. f_aux(ijk).lt.emf)) cycle
!
! skip qsrm > 0.0 cells since they may have temperature changes
! because of the introduction of colder or hotter fluid and not
! because of heating or cooling
!
            if(iqsr.gt.0) then
              if(qsrm(ijk).gt.zero) cycle
            endif
!
            rhoe_sav=rhoe(ijk)
            f_sav=f(ijk)
            tn_sav=tn(ijk)
!
            rhoe(ijk)=rhoe_aux(ijk)
            f(ijk)=f_aux(ijk)
!
! compute temperature and liquid fraction at the previous time-step
!
            if(ishrnk.eq.1) then
              tn(ijk)=teval(ijk)
              fliqold=fliqcl(ijk)
              tnold=tn(ijk)
            else
              tnold=teval(ijk)
            endif
!
            rhoe(ijk)=rhoe_sav
            f(ijk)=f_sav
            tn(ijk)=tn_sav
!
! rhoe could be incremented with heat sources in cdsoli, therefore,
! must call teval to get the latest temperature
!
            tnnew=teval(ijk)
!
! ### shrinkage computed from the change in temperature
!
            if(ishrnk.eq.1) then
!
! compute volume source for the dynamic shrinkage model
!
              if(ifrho.eq.1) then
!
                  rho1=rho1cl(p(ijk),tnold)
                  rho2=rho1cl(p(ijk),tnnew)
!
                  qsr_liq_shr=one-rho2/rho1
!
! use liquid phase density at current temperature TN to compute
! the solidification shrinkage
!
                  density_ratio=rhofs/rho2
!
              else
!
                  qsr_liq_shr=zero
!
              endif
!
              if(irhof.eq.0 .and. fliq(ijk).lt.one) then
!
! ### shrinkage computed from the change in solid fraction
!
! convert mass fraction of solid to volume fraction of solid
!
                  frslv=(one-fliq(ijk))/ & 
     &                  (one-fliq(ijk)+fliq(ijk)*density_ratio)
                  frslvn=(one-fliqold)/ & 
     &                   (one-fliqold+fliqold*density_ratio)
!
!    calculate qsr_sld_shr - cell volume change per unit open cell volume
!    caused by shrinkage: negative for shrinkage, positive for expansion.
!
                  denom=frslv*(density_ratio-one)+one
                  if(denom.gt.ztest) then
                    qsr_sld_shr=(one-density_ratio)*(frslv-frslvn)/denom
                  else
                    qsr_sld_shr=(one-density_ratio)*(frslv-frslvn)
                  endif
!
! non-tabular density variation with temperature applies only to the liquid phase, so
! we multiply it here by the volume fraction of liquid,
! while tabular density variation with temperature applies to the solid-liquid mixture
! across the whole solidification range and, therefore, does not need to be scaled by
! liquid fraction.
!
                  qsr_liq_shr=(one-frslv)*qsr_liq_shr
!
              else
!
                  qsr_sld_shr=zero
!
              endif
!
              qsrthexp(ijk)=(qsr_liq_shr+qsr_sld_shr)*rdeltc
!
            else
!
! compute volume source for thermal expansion model
!
              rho2=rho1cl(p(ijk),tnnew)
              rho1=rho1cl(p(ijk),tnold)
!
              rcsqfv=one-rho2/rho1
!
              if(nmat.eq.2) then
                if(icmprs.eq.0) then
                  rho2_2=rho2cl(p(ijk),tnnew,zero)
                  rho1_2=rho2cl(p(ijk),tnold,zero)
                  rcsqfv=f(ijk)*rcsqfv+(one-f(ijk))*(one-rho2_2/rho1_2)
                else
                  rcsqfv=f(ijk)*rcsqfv
                endif
              endif
!
              qsrthexp(ijk)=qsrthexp(ijk)+rcsqfv*rdeltc
!
            endif
!
 2500     continue
 2501   continue
 2502 continue
!
      call et_off(45)
!
      endif
!
      return
      end
