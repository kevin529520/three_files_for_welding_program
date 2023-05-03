!> \brief 
!>    Add particle of class np_gas to replace collapsed voids
!
      subroutine void_2particle
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
      use mblock_module
!
      use arrays_module
!
      use arrayp_module
!
      use meshcb_module
!
      use avgco_module
!
      use voids_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/phiou.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/diag.f90'
      include '../comdeck/edit.f90'
      include '../comdeck/grfdat.f90'
      include '../comdeck/state.f90'
      include '../comdeck/acctrm.f90'
      include '../comdeck/pardat.f90'
      include '../comdeck/scala.f90'
      include '../comdeck/dumn.f90'
!
      if(num_particles(np_void).ge.ip_nmax(np_void)) return
!
      block_loop: do nbl=1,nblcks
          call lblock(0)
!
          do ijk=ijklim(nbl,1),ijklim(nbl,2)
!
              if(vf(ijk).lt.em6) cycle
              include '../comdeck/ijk2str.f90'
              if(icstat(ijk_str).le.0) cycle
!
! * * * determine if a void collapsed
!
! skip cells with void now or no void before
!
              if(nf(ijk).gt.0 .or. nfn(ijk).lt.8) cycle
!
              include '../comdeck/mijk.f90'
              include '../comdeck/pijk.f90'
!
! check to see if the void moved from cell ijk to a neighbor cell
!
              if(afr(imjk).gt.em6 .and. nf(imjk).gt.0) then
                call nvcell(imjk,nvoid,istat)
                if(istat.eq.1) then
                  if(nfc(nvoid).eq.nfn(ijk) .or. & 
     &           nfnc(nfn(ijk)).eq.nvoid) cycle ! check if void nvoid is the same as nfn
                endif
              endif
              if(afr(ijk).gt.em6 .and. nf(ipjk).gt.0) then
                call nvcell(ipjk,nvoid,istat)
                if(istat.eq.1) then
                  if(nfc(nvoid).eq.nfn(ijk) .or. & 
     &           nfnc(nfn(ijk)).eq.nvoid) cycle ! check if void nvoid is the same as nfn
                endif
              endif
              if(afb(ijmk).gt.em6 .and. nf(ijmk).gt.0) then
                call nvcell(ijmk,nvoid,istat)
                if(istat.eq.1) then
                  if(nfc(nvoid).eq.nfn(ijk) .or. & 
     &           nfnc(nfn(ijk)).eq.nvoid) cycle ! check if void nvoid is the same as nfn
                endif
              endif
              if(afb(ijk).gt.em6 .and. nf(ijpk).gt.0) then
                call nvcell(ijpk,nvoid,istat)
                if(istat.eq.1) then
                  if(nfc(nvoid).eq.nfn(ijk) .or. & 
     &           nfnc(nfn(ijk)).eq.nvoid) cycle ! check if void nvoid is the same as nfn
                endif
              endif
              if(aft(ijkm).gt.em6 .and. nf(ijkm).gt.0) then
                call nvcell(ijkm,nvoid,istat)
                if(istat.eq.1) then
                  if(nfc(nvoid).eq.nfn(ijk) .or. & 
     &           nfnc(nfn(ijk)).eq.nvoid) cycle ! check if void nvoid is the same as nfn
                endif
              endif
              if(aft(ijk).gt.em6 .and. nf(ijkp).gt.0) then
                call nvcell(ijkp,nvoid,istat)
                if(istat.eq.1) then
                  if(nfc(nvoid).eq.nfn(ijk) .or. & 
     &           nfnc(nfn(ijk)).eq.nvoid) cycle ! check if void nvoid is the same as nfn
                endif
              endif
!
! * * * void did collapse in cell ijk
!
              i=i_str(ijk_str)
              j=j_str(ijk_str)
              k=k_str(ijk_str)
!
              nvoid=nfn(ijk) ! not presently used
!
!        void collapsed: add void particle
!
              kp_void=num_particles(np_void)+1
!
              num_particles(np_void)=kp_void
! place particle in the middle of the cell
              void_particle(kp_void)%x_part=xi(i)
              void_particle(kp_void)%y_part=yj(j)
              void_particle(kp_void)%z_part=zk(k)
              void_particle(kp_void)%t_part=zero
              void_particle(kp_void)%ijk_str_part=ijk_str
              void_particle(kp_void)%nbl_part=nbl
              void_particle(kp_void)%motn_part=1
              void_particle(kp_void)%istag_part=0
              void_particle(kp_void)%ibtag_part=0
              void_particle(kp_void)%u_part= & 
     & (afr(imjk)*u(imjk)+afr(ijk)*u(ijk))/(afr(imjk)+afr(ijk)+ztest)
              void_particle(kp_void)%v_part= & 
     & (afb(ijmk)*v(ijmk)+afb(ijk)*v(ijk))/(afb(ijmk)+afb(ijk)+ztest)
              void_particle(kp_void)%w_part= & 
     & (aft(ijkm)*w(ijkm)+aft(ijk)*w(ijk))/(aft(ijkm)+aft(ijk)+ztest)
! void volume in cell
              if(p_diameter(np_void).lt.ztest) then
                void_volume=(one-f(ijk))*vf(ijk)* & 
     &                      delx(i)*dely(j)/rri(i)*delz(k)
                if(void_volume.gt.ztest) then
                  void_diameter=(six*void_volume/cpi)**ob3
                else
                  void_diameter=0.001d0*delmnn
                endif
                void_particle(kp_void)%size_part=void_diameter
              else
                void_particle(kp_void)%size_part=p_diameter(np_void)
              endif
              void_particle(kp_void)%charge_part=zero
              if(ifenrg.gt.0) then
                  void_particle(kp_void)%temp_part=tn(ijk) 
              else
                  void_particle(kp_void)%temp_part=tv2m
              endif
              void_particle(kp_void)%pressure_part=p(ijk)  
! calculate particle density              
              if(ipvoid_type.eq.1 .or.  & 
     &  void_particle(kp_void)%pressure_part.lt.ztest .or. & 
     &  void_particle(kp_void)%temp_part+t2abs.lt.ztest) then
                  if(p_density(np_void).gt.ztest) then
                    void_particle(kp_void)%density_part= & 
     &                  p_density(np_void)
                  else
                    void_particle(kp_void)%density_part=0.001d0*rhof
                  endif
              else
                  void_particle(kp_void)%density_part= & 
     &               void_particle(kp_void)%pressure_part/ & 
     &    (pvoid_rf*(void_particle(kp_void)%temp_part+t2abs))
              endif
!              
              id_max(np_void)=id_max(np_void)+1
              void_particle(kp_void)%id_species=1
              void_particle(kp_void)%id_part=id_max(np_void)
!
              num_particles_g(np_void)=num_particles_g(np_void)+1
!
              if(num_particles(np_void).ge.ip_nmax(np_void)) & 
     &                 exit block_loop
!
          enddo
      enddo block_loop
!
      return
      end
