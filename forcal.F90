                                   subroutine forcal
!c-----------------------------------------------------------------------
!c   evaluate additional forces on fluid -- user customizable
!c-----------------------------------------------------------------------
!c
!c    ***************************************************************
!c    **                          notice                           **
!c    **  this subprogram contains flow science, inc. proprietary  **
!c    **       trade secret and confidential information.          **
!c    **                                                           **
!c    **                unauthorized use prohibited                **
!c    **           copyright 1985-2005 flow science, inc.          **
!c    ***************************************************************
!c
      use arrays_module
!c
      use meshcb_module
!c
      use avgco_module
!c
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/dparam.f90'
      include '../comdeck/const.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/edit.f90'
      include '../comdeck/phiou.f90'
      include '../comdeck/scala.f90'
       include '../comdeck/dumn.f90'
!c      include '../comdeck/multifile.f90'
!c
!c      fsigx(ijk)= x component of added force on fluid in cell ijk
!c      fsigy(ijk)= y component of added force on fluid in cell ijk
!c      fsigz(ijk)= z component of added force on fluid in cell ijk
!c
!c      these forces are either zero or contain non-zero contributions
!c      from thermo-capillary surface forces if ifsten<>0
!c
!c      the forces are located at grid cell centers
!c
!c ------- definition of quantities for customizing --------
!c      variable        description
!c      --------        --------------------------------
!c       ijk            current cell index
!c       ipjk             cell to right
!c       imjk             cell to left
!c       ijpk             cell to back
!c       ijmk             cell to front
!c       ijkp             cell to top
!c       ijkm             cell to bottom
!c       i                current x index
!c       j                current y index
!c       k                current z index
!c
!c       t                time
!c       delt             time step size
!c       x(i)             mesh coordinate at right of cell ijk
!c       xi(i)            cell ijk center
!c       y(j)             mesh coordinate at back of cell ijk
!c       yj(j)            cell ijk center
!c       z(k)             mesh coordinate at top of cell ijk
!c       zk(k)            cell ijk center
!c       delx(i)          cell size in x direction
!c       dely(j)          cell size in y direction
!c       delz(k)          cell size in z direction
!c       rri(i)           correction factor for cylindrical coordinates
!c                        i.e., delta y at x(i) is dely(j)/rri(i)
!c       vf(ijk)        open volume fraction in cell
!c       afr(ijk)       open area fraction at right face
!c       afb(ijk)       open area fraction at back face
!c       aft(ijk)       open area fraction at top face
!c
!c       u(ijk)         x velocity at right face
!c       v(ijk)         y velocity at back face
!c       w(ijk)         z velocity at top face
!c       fn(ijk)        fluid fraction in cell at beginning of cycle
!c       p(ijk)         pressure in cell
!c       tn(ijk)        temperature in cell
!c       rhoe(ijk)      density*specific energy in cell
!c       arint(ijk)     free surface area in cell
!c       rho(ijk)       density in cell (only for variable density)
!c
!c       nf(ijk)        free surface indicator in cell
!c          =0           interior fluid cell
!c          =1           surface cell - fluid at left
!c          =2           surface cell - fluid at right
!c          =3           surface cell - fluid at front
!c          =4           surface cell - fluid at back
!c          =5           surface cell - fluid at bottom
!c          =6           surface cell - fluid at top
!c          =7           surface cell - cavitating cell
!c          >=8          void cell -- void id number
!c
!c       nsc            number of scalars
!c       sclr(ijk,ns)   concentration of scalar ns at cell ijk
!c                      after advection and diffusion
!c                      (update this variable to change scalar
!c                      concentration)
!c       sclrn(ijk,ns)  concentration of scalar ns at cell ijk
!c                      at beginning of time step
!c
!c -------  initialize force arrays if not done in thermocapillary model
!c
      if(ifsten.eq.0) then
        do 100 k=1,km1
        do 100 j=1,jm1
        do 100 i=1,im1
        include '../comdeck/ijk.f90'
          fsigx(ijk)=zero
          fsigy(ijk)=zero
          fsigz(ijk)=zero
  100   continue
      endif
!c
!c         ... enter changes here ...
!c
!c
!c ============ SPECIAL CUSTOMIZATION FOR APPLIED SURFACE PRESSURE=======
!c
!c ------ Add applied pressure to fluid surface.  This routine
!c        computes a surface force arising from an applied pressure,
!c        e.g., by a gas jet or plasma arc.
!c
!c ---  loop over real cells (set boundary cells in subroutine bc)

!c ============= User Input ======================================================
!c         U=400.0
!c         I=30.0
!c         M=0.65
        apmax=900
        af=0.001
        ar=0.001
         bh=0.001
         apdia=0.0021
         fapdia=0.0019
         Fxx=29.23
         Fyy=29.23	
         Fzz=0.0265
         H=0.003
         apx=dum10
         apy=dum11
         apz=dum12
         apnrx=0.0
         apnry=0.0
         apnrz=-1.0
!c ============= End User Input ==================================================
!c
       if(apdia.lt.ztest) goto 200
!c
       do 1000 k=kprb,kprt
       do 1000 j=jprf,jprbk
       do 1000 i=iprl,iprr
!c -----  calculate current cell index
         include '../comdeck/ijk.f90'
!c
!c ------   skip calculation for completely blocked cells
         if(vf(ijk).lt.em6) cycle
!c ------ skip passive cells
         if(icstat(ijk).le.0) cycle
!c ------   skip non-surface cells
!c
!c ------  Other tests could be added to restrict application of surface pressures
!c ------ Here we assume the applied pressure is from a jet with a center location
!c        (apx,apy,apz) and having an axial normal vector (apnrx,apnry,apnrz) in 
!c        the direction of the jet.
!c        The pressure distribution is assumed to be Gaussian about the axis of the
!c        jet with amplitude apmax and standard deviation diameter apdia.       

!c        Get outward normal to free surface
!c
         call srfnorm(f,ijk,snx,sny,snz,sar,icount_mpi,0)
!c
         dotnorm=snx*apnrx+sny*apnry+snz*apnrz
         if(dotnorm.ge.zero) goto 1000
!c
!c ------ normal distance from center of distribution is disnorm.
!c        this distance is used to impose a shape factor on the pressure,
!c        i.e., a Gaussian with radius from the center:
         rsq=(xi(i)-apx-dum13*t)*(xi(i)-apx-dum13*t)+(zk(k)-apz)*(zk(k)-apz)
!c ------ add y distance if not cylindrical coordinates
         if(cyl.lt.half) rsq=rsq+(yj(j)-apy)*(yj(j)-apy)
         rdotz=(xi(i)-apx-dum13*t)*apnrx+(yj(j)-apy)*apnry+(zk(k)-apz)*apnrz
         disnorm=sqrt(max(zero,rsq-rdotz*rdotz))
!c
!c ------ evaluate Guassian distribution
         if((xi(i)-dum13*t).ge.apx)  then
         argexp=(xi(i)-apx-dum13*t)**2/af**2+(yj(j)-apy)**2/bh**2
         else
         argexp=(xi(i)-apx-dum13*t)**2/ar**2+(yj(j)-apy)**2/bh**2    
         endif
!c          argexp=(disnorm/apdia)**2
          fargexp=(disnorm/fapdia)**2
         if(nf(ijk).eq.0 .or. nf(ijk).ge.8) goto 1100          
         papp=apmax*exp(-3*argexp)
1100  continue
!c ------ compute surface force magnitude as pressure*area
!c        (factor of two needed because of force averaging)
         pforce=two*arint(ijk)*papp*dotnorm
!c ------ evaluate electromagnetic force distribution
         Fx=-(Fxx/disnorm)*exp(-0.5*fargexp)*(1-exp(-0.5*fargexp))*((zk(k)/H)**2)*(xi(i)-apx-dum13*t)/disnorm
     
         Fy=-(Fyy/disnorm)*exp(-0.5*fargexp)*(1-exp(-0.5*fargexp))*((zk(k)/H)**2)*(yj(j)-apy)/disnorm
         
         Fz=-(Fzz/disnorm**2)*((1-exp(-0.5*fargexp))**2)*(zk(k)/H)
         
!c ------ compute electromagnetic force magnitude as pressure*volume         
         eforcex=delx(i)*dely(j)*delz(k)*Fx
         eforcey=delx(i)*dely(j)*delz(k)*Fy
         eforcez=delx(i)*dely(j)*delz(k)*Fz
!c ------ update cell centered forces in surface cells 
         fsigx(ijk)=fsigx(ijk)+pforce*snx+eforcex*f(ijk)
         fsigy(ijk)=fsigy(ijk)+pforce*sny+eforcey*f(ijk)
         fsigz(ijk)=fsigz(ijk)+pforce*snz+eforcez*f(ijk)
!c
!c ------ There is a possible instability here because of forces applied
!c        in directions normal to the axis of pressure jet.  This happens
!c        because the applied pressure is constant and does not adjust
!c        to surface elevation changes.  The instability will not occur
!c        if surface tension or gravity forces are large enough to smooth
!c        any horizontal variations in elevation caused by the applied pressure.  
!c
 1000 continue
!c
  200 continue
!c
      return
      end