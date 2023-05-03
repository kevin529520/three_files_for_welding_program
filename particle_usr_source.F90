!> \brief 
!>    Place user particlesin the domain
!
      subroutine particle_usr_source
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
      use obsijk_module
!
      use avgco_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/dparam.f90'
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
      include '../comdeck/bafdef.f90'
      include '../comdeck/obsd.f90'
      include '../comdeck/dumn.f90'
!
! *** input:
!
! ip_class - particle class:
!
!          <= NSC_PART, - user particle class
!
! *** ADD USER CODE HERE
!
      RETURN ! comment out this line if making chanes below
!
! EXAMPLE. ALL VARIABLES WITH 'placement' MUST BE DEFINED BY USER:
!
      do ip_class=1,nsc_part
!
        if(num_particles(ip_class).ge.ip_nmax(ip_class)) cycle ! particle count in class ip_class reached maximum
!
! x_placement, y_placement, z_placement - Cartesian particle placement coordinates
! u_placement, v_placement, w_placement - initial particle velocities, mesh coordinates.
!
        x_placement=zero
        y_placement=zero
        z_placement=zero
        u_placement=zero
        v_placement=zero
        w_placement=zero
!
! find cell index for particle CARTESIAN location x_placement, y_placement, z_placement
!
        call gpcijk(x_placement,y_placement,z_placement, & 
     &                     xout,yout,ijk_str,nblp,1,istat)
!
! check if cell found and if it's legit
!
        if(ijk_str.gt.0) then   ! a structured mesh cell found
          ijk=ijk_str2unstr(ijk_str)
          if(vf(ijk).gt.em6 .and. ijk.lt.ijklim_bc) then ! active fluid cell
!
! ADD A PARTICLE FOR USER CLASS #IP_CLASS
!
! store the user-class particle
              kp_user=num_particles(ip_class)+1                       ! increment particle count
              num_particles(ip_class)=kp_user
!
              user_particle(kp_user,ip_class)%x_part=xout             ! x-coordinate (r in cylindrical)
              user_particle(kp_user,ip_class)%y_part=yout             ! y-coordiante (rad*x(im1) in cylindrical)
              user_particle(kp_user,ip_class)%z_part=z_placement      ! z-coordinate
              user_particle(kp_user,ip_class)%t_part=zero             ! life time
              user_particle(kp_user,ip_class)%ijk_str_part=ijk_str    ! structured index of containing cell
              user_particle(kp_user,ip_class)%nbl_part=nblp           ! containing mesh block
              user_particle(kp_user,ip_class)%motn_part= & 
     &                                     ip_motion(ip_class)        ! movement flag (0 - stationary, 1 - moving)
              user_particle(kp_user,ip_class)%istag_part=0            ! source tag
              user_particle(kp_user,ip_class)%ibtag_part=0            ! initial block tag
              user_particle(kp_user,ip_class)%u_part=u_placement      ! initial x-velocity
              user_particle(kp_user,ip_class)%v_part=v_placement      ! initial y-velocity
              user_particle(kp_user,ip_class)%w_part=w_placement      ! initial z-velocity
              user_particle(kp_user,ip_class)%size_part= & 
     &                                              p_diameter(ip_class)  ! diameter
              user_particle(kp_user,ip_class)%density_part= & 
     &                                               p_density(ip_class)  ! density
              user_particle(kp_user,ip_class)%charge_part= & 
     &                                                p_charge(ip_class)  ! electric charge
              user_particle(kp_user,ip_class)%temp_part= & 
     &                                          p_temperature(ip_class)   ! initial temperature
              id_max(ip_class)=id_max(ip_class)+1                         ! particle id counter
              user_particle(kp_user,ip_class)%id_part=id_max(ip_class)    ! particle id
!
! free paticle attributes
!
              do nn=1,n_usr_attr(ip_class)
                  user_particle(kp_user,ip_class)%attribute(nn)=zero
              enddo
!
              num_particles_g(ip_class)=num_particles_g(ip_class)+1   ! number of paticles generated
          endif
!
        endif
!
      enddo
!
      return
      end
