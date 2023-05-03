!> \brief 
!>    Specify additional force on a material particle
!
      subroutine particle_usr_force(ip_class, & 
     &                   x_particle,y_particle,z_particle, & 
     &                   u_particle,v_particle,w_particle, & 
     &                   ijk_str_particle,motn_particle, & 
     &                   size_particle,density_particle,charge_particle, & 
     &                   id_particle,ufrcx,ufrcy,ufrcz)
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
!  x_particle - particle x-coordiante (radius in cylindrical coordiantes)
!  y_particle - particle y-coordinate (rad*x(im1) in cylindrical coordiantes)
!  z_particle - particle z-coordinate
!
!  u_particle - particle x-velocity (radial in cylindrical coordiantes)
!  v_particle - particle y-velocity (linear azimuthal in cylindrical coordiantes)
!  w_particle - particle z-velocity
!
!  ijk_str_particle - structured index of the cell containing the particle
!  motn_particle    - motion flag:
!          = 0 - stationary particle (not called here)
!          = 1 - moving particle
!  size_particle    - particle diameter
!  density_particle - particle density
!  charge_particle  - particle electric charge
!  id_particle      - particle unique id ithing the class
!
! *** output:
!
!  ufrcx,ufrcy,ufrcz - user force components
!
! *** ADD USER CODE HERE
!
!   ijk=ijk_str2unstr(ijk_str_particle) - fluid cell index
!   i=i_str(ijk_str_particle)         - cell x-direction index
!   j=j_str(ijk_str_particle)         - cell y-direction index
!   k=k_str(ijk_str_particle)         - cell z-direction index
!
      ufrcx=zero
      ufrcy=zero
      ufrcz=zero
!
      return
      end
