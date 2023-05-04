subroutine qsadd

!c
!c     this subroutine is called when nsc>0.  the call is near the   (NSC	0	Number of user-defined scalar species)

!c     end of the cycle, after the pressure/velocity update and fluid

!c     advection and diffusion,

!c     but before new cells are initialized, nf's are set,

!c     and the chemistry routine is called.

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
      use mblock_module

!c
      use arrays_module

!c
      use arrayp_module

!c
      use meshcb_module

!c
      use voids_module
!c
      use obsijk_module
!c
#ifdef SINGLE
      
      include '../comdeck/precis4.f90'
#else
     
      include '../comdeck/precis.f90'
#endif
     
      include '../comdeck/params.f90'
      include '../comdeck/dparam.f90'
      include '../comdeck/bafdef.f90'
      include '../comdeck/obsd.f90'
      
      include '../comdeck/cntrl.f90'
     
      include '../comdeck/const.f90'
      
      include '../comdeck/dumn.f90'
     
      include '../comdeck/phiou.f90'
      
      include '../comdeck/scala.f90'
     
      include '../comdeck/state.f90'
      
      include '../comdeck/pardat.f90'

!c
!c             scalar species sources and sinks

!c
!c     (not currently implemented)

!c
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

!c
!c       nbl              current mesh block number

!c
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

!c
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
!c
!c      skip over if no scalars exist and this subroutine is used for

!c         scalar sources

!c      if(nsc.eq.0) return

!c
!c ---  loop over real cells (set boundary cells in subroutine bc)

!c       do 100 k=kprb,kprt

!c       do 100 j=jprf,jprbk

!c       do 100 i=iprl,iprr

!c -----  calculate current cell index

!c        include '../comdeck/ijk.f'

!c ------  skip calculation for completely blocked cells

!c         if(vf(ijk).lt.em6) goto 100

!c --------  calculate "neighbor indices"

!c           include '../comdeck/mijk.f'

!c           include '../comdeck/pijk.f'

!c ----------  skip empty (void) cells

!c             if(fn(ijk).lt.emf .and. nmat.eq.1) go to 100

!c
!c ====================SPECIAL CUSTOMIZATION ===========================
!c
!c ------ Add thermal energy beam source to fluid surfaces.  A phantom 
!c        obstacle is used to define the beam location.  This routine
!c        computes an energy depostion where a fluid surface overlaps
!c        the phantom obstacle.  The obstacle is defined as a type ifob=2
!c        with ospin<>0, oadrg=0, and obdrg=0.  (Note: ospin<>0 is needed.)
!c        If a surface cell is less than half filled (f<0.5) energy is 
!c!        distributed over surface cell plus its neighbor
!c
!c       droplet
	timdrp=0.01   !c  0.01秒滴一滴
	temdrop=dum18
	udrop=0.0
	vdrop=0.0
	wdrop=-0.5
    rdrop=0.001
    !c rdrop=0.001  !c半径
	nsub=10	
    
       !c ydrop=dum11  ! dum11=0.063     dum12=0.001,    dum13=0.005,  dum10=0.0,dum17=0.1,    dum18=1800,
	   !c zdrop=dum12+0.001
       !c xdrop=dum10+dum13*t
   
    if(t.lt.40)then
       
	    ydrop=dum11  ! dum11=0.063     dum12=0.001,    dum13=0.005,  dum10=0.0,dum17=0.1,    dum18=1800,
	    zdrop=dum12+0.001
	    xdrop=dum10+dum13*t
    elseif(t.ge.40 .AND. t.lt.80) then
        ydrop=dum11  ! dum11=0.063     dum12=0.001,    dum13=0.005,  dum10=0.0,dum17=0.1,    dum18=1800,
	    zdrop=dum12+0.007!+0.002  更改
	    xdrop=dum10+dum13*(t-40)   !c 更改 t-40
    elseif(t.ge.80 .AND. t.lt.120) then
        ydrop=dum11  ! dum11=0.063     dum12=0.001,    dum13=0.005,  dum10=0.0,dum17=0.1,    dum18=1800,
	    zdrop=dum12+0.014!+0.002  更改
	    xdrop=dum10+dum13*(t-80)   !c 更改 t-40   
    
    endif
    
    
!c       add droplet
	if(dum17.lt.t-half*delt.or.dum17.gt.t+half*delt) goto 1100
	dum17=dum17+timdrp
	rdrsqr=rdrop*rdrop
	ncube=nsub*nsub*nsub
	kbd=0
	ked=0
	jbd=0
	jed=0
	ibd=0
	ied=0
	do 300 k=kprb,kprt!ckprb	first real cell index in z-direction  kprt	last real cell index in z-direction
		if(kbd.eq.0.and.z(k).gt.zdrop-rdrop) kbd=k  !z(k) mesh coordinate at top of cell ijk
		if(ked.eq.0.and.z(k).gt.zdrop+rdrop) ked=k
  300 continue
	if(kbd.eq.0) goto 1100
	if(ked.eq.0) ked=kprt
	if(cyl.gt.0.5) then
		jbd=jprf
		jed=jprbk
		ibd=iprl
		ied=iprr
		goto 550
	endif
	do 400 j=jprf,jprbk
		if(jbd.eq.0.and.y(j).gt.ydrop-rdrop) jbd=j
		if(jed.eq.0.and.y(j).gt.ydrop+rdrop) jed=j
  400 continue
	if(jbd.eq.0) goto 1100
	if(jed.eq.0) jed=jprbk
	do 500 i=iprl,iprr
		if(ibd.eq.0.and.x(i).gt.xdrop-rdrop) ibd=i
		if(ied.eq.0.and.x(i).gt.xdrop+rdrop) ied=i
  500 continue
	if(ibd.eq.0) goto 1100
	if(ied.eq.0) ied=iprr
  550 continue
	do 2000 k=kbd,ked
	do 2000 j=jbd,jed
	do 2000 i=ibd,ied
	include '../comdeck/ijk.f90'
	if(vf(ijk).lt.em6) goto 2000
	vdv=one/float(ncube)
	xdx=(x(i)-x(i-1))/float(nsub)
	ydy=(y(j)-y(j-1))/float(nsub)
	zdz=(z(k)-z(k-1))/float(nsub)
	xc=x(i-1)-half*xdx
	yc=y(j-1)-half*ydy
	zc=z(k-1)-half*zdz
	xcc=xc
	ycc=yc
	do 600 nk=1,nsub
		zc=zc+zdz
	do 580 nj=1,nsub
		ycc=ycc+ydy
	do 560 ni=1,nsub
		xcc=xcc+xdx
		vdvc=vdv
		if(cyl.gt.0.5) then
			vdvc=vdv*xcc/xi(i)
			th=ycc/x(iprr)
			ycc=xcc*sin(th)
			xcc=xcc*cos(th)
		endif
		dissqr=(xcc-xdrop)**2+(ycc-ydrop)**2+(zc-zdrop)**2
		if(dissqr.gt.rdrsqr) goto 560
		f(ijk)=f(ijk)+vdvc
  560 continue
		xcc=xc
  580 continue
		ycc=yc
  600 continue
	if(f(ijk).gt.one) f(ijk)=one
	if(f(ijk).lt.emf) then
		f(ijk)=zero
		goto 2000
	endif
	u(ijk-1)=udrop
	u(ijk)=udrop
	v(ijk-ii1)=vdrop
	v(ijk)=vdrop
	w(ijk-ii2)=wdrop
	w(ijk)=wdrop
	tn(ijk)=temdrop
	if(ifenrg.gt.0) rhoe(ijk)=f(ijk)*rhof*e1cal(0.0,tn(ijk),0.0,0.0)
 2000 continue
 1100 continue 
      return
      endsubroutine qsadd

!c
!c     this subroutine is called when nsc>0.  the call is near the   (NSC	0	Number of user-defined scalar species)

!c     end of the cycle, after the pressure/velocity update and fluid

!c     advection and diffusion,

!c     but before new cells are initialized, nf's are set,

!c     and the chemistry routine is called.

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
      use mblock_module

!c
      use arrays_module

!c
      use arrayp_module

!c
      use meshcb_module

!c
      use voids_module
!c
      use obsijk_module
!c
#ifdef SINGLE
      
      include '../comdeck/precis4.f90'
#else
     
      include '../comdeck/precis.f90'
#endif
     
      include '../comdeck/params.f90'
      include '../comdeck/dparam.f90'
      include '../comdeck/bafdef.f90'
      include '../comdeck/obsd.f90'
      
      include '../comdeck/cntrl.f90'
     
      include '../comdeck/const.f90'
      
      include '../comdeck/dumn.f90'
     
      include '../comdeck/phiou.f90'
      
      include '../comdeck/scala.f90'
     
      include '../comdeck/state.f90'
      
      include '../comdeck/pardat.f90'

!c
!c             scalar species sources and sinks

!c
!c     (not currently implemented)

!c
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

!c
!c       nbl              current mesh block number

!c
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

!c
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
!c
!c      skip over if no scalars exist and this subroutine is used for

!c         scalar sources

!c      if(nsc.eq.0) return

!c
!c ---  loop over real cells (set boundary cells in subroutine bc)

!c       do 100 k=kprb,kprt

!c       do 100 j=jprf,jprbk

!c       do 100 i=iprl,iprr

!c -----  calculate current cell index

!c        include '../comdeck/ijk.f'

!c ------  skip calculation for completely blocked cells

!c         if(vf(ijk).lt.em6) goto 100

!c --------  calculate "neighbor indices"

!c           include '../comdeck/mijk.f'

!c           include '../comdeck/pijk.f'

!c ----------  skip empty (void) cells

!c             if(fn(ijk).lt.emf .and. nmat.eq.1) go to 100

!c
!c ====================SPECIAL CUSTOMIZATION ===========================
!c
!c ------ Add thermal energy beam source to fluid surfaces.  A phantom 
!c        obstacle is used to define the beam location.  This routine
!c        computes an energy depostion where a fluid surface overlaps
!c        the phantom obstacle.  The obstacle is defined as a type ifob=2
!c        with ospin<>0, oadrg=0, and obdrg=0.  (Note: ospin<>0 is needed.)
!c        If a surface cell is less than half filled (f<0.5) energy is 
!c!        distributed over surface cell plus its neighbor
!c
!c       droplet
	timdrp=0.01   !c  0.01秒滴一滴
	temdrop=dum18
	udrop=0.0
	vdrop=0.0
	wdrop=-0.5
    rdrop=0.001
    !c rdrop=0.001  !c半径
	nsub=10	
    
       !c ydrop=dum11  ! dum11=0.063     dum12=0.001,    dum13=0.005,  dum10=0.0,dum17=0.1,    dum18=1800,
	   !c zdrop=dum12+0.001
       !c xdrop=dum10+dum13*t
   
    if(t.lt.40)then
       
	    ydrop=dum11  ! dum11=0.063     dum12=0.001,    dum13=0.005,  dum10=0.0,dum17=0.1,    dum18=1800,
	    zdrop=dum12+0.001
	    xdrop=dum10+dum13*t
    elseif(t.ge.40 .AND. t.lt.80) then
        ydrop=dum11  ! dum11=0.063     dum12=0.001,    dum13=0.005,  dum10=0.0,dum17=0.1,    dum18=1800,
	    zdrop=dum12+0.007!+0.002  更改
	    xdrop=dum10+dum13*(t-40)   !c 更改 t-40
    elseif(t.ge.80 .AND. t.lt.120) then
        ydrop=dum11  ! dum11=0.063     dum12=0.001,    dum13=0.005,  dum10=0.0,dum17=0.1,    dum18=1800,
	    zdrop=dum12+0.014!+0.002  更改
	    xdrop=dum10+dum13*(t-80)   !c 更改 t-40   
    
    endif
    
    
!c       add droplet
	if(dum17.lt.t-half*delt.or.dum17.gt.t+half*delt) goto 1100
	dum17=dum17+timdrp
	rdrsqr=rdrop*rdrop
	ncube=nsub*nsub*nsub
	kbd=0
	ked=0
	jbd=0
	jed=0
	ibd=0
	ied=0
	do 300 k=kprb,kprt!ckprb	first real cell index in z-direction  kprt	last real cell index in z-direction
		if(kbd.eq.0.and.z(k).gt.zdrop-rdrop) kbd=k  !z(k) mesh coordinate at top of cell ijk
		if(ked.eq.0.and.z(k).gt.zdrop+rdrop) ked=k
  300 continue
	if(kbd.eq.0) goto 1100
	if(ked.eq.0) ked=kprt
	if(cyl.gt.0.5) then
		jbd=jprf
		jed=jprbk
		ibd=iprl
		ied=iprr
		goto 550
	endif
	do 400 j=jprf,jprbk
		if(jbd.eq.0.and.y(j).gt.ydrop-rdrop) jbd=j
		if(jed.eq.0.and.y(j).gt.ydrop+rdrop) jed=j
  400 continue
	if(jbd.eq.0) goto 1100
	if(jed.eq.0) jed=jprbk
	do 500 i=iprl,iprr
		if(ibd.eq.0.and.x(i).gt.xdrop-rdrop) ibd=i
		if(ied.eq.0.and.x(i).gt.xdrop+rdrop) ied=i
  500 continue
	if(ibd.eq.0) goto 1100
	if(ied.eq.0) ied=iprr
  550 continue
	do 2000 k=kbd,ked
	do 2000 j=jbd,jed
	do 2000 i=ibd,ied
	include '../comdeck/ijk.f90'
	if(vf(ijk).lt.em6) goto 2000
	vdv=one/float(ncube)
	xdx=(x(i)-x(i-1))/float(nsub)
	ydy=(y(j)-y(j-1))/float(nsub)
	zdz=(z(k)-z(k-1))/float(nsub)
	xc=x(i-1)-half*xdx
	yc=y(j-1)-half*ydy
	zc=z(k-1)-half*zdz
	xcc=xc
	ycc=yc
	do 600 nk=1,nsub
		zc=zc+zdz
	do 580 nj=1,nsub
		ycc=ycc+ydy
	do 560 ni=1,nsub
		xcc=xcc+xdx
		vdvc=vdv
		if(cyl.gt.0.5) then
			vdvc=vdv*xcc/xi(i)
			th=ycc/x(iprr)
			ycc=xcc*sin(th)
			xcc=xcc*cos(th)
		endif
		dissqr=(xcc-xdrop)**2+(ycc-ydrop)**2+(zc-zdrop)**2
		if(dissqr.gt.rdrsqr) goto 560
		f(ijk)=f(ijk)+vdvc
  560 continue
		xcc=xc
  580 continue
		ycc=yc
  600 continue
	if(f(ijk).gt.one) f(ijk)=one
	if(f(ijk).lt.emf) then
		f(ijk)=zero
		goto 2000
	endif
	u(ijk-1)=udrop
	u(ijk)=udrop
	v(ijk-ii1)=vdrop
	v(ijk)=vdrop
	w(ijk-ii2)=wdrop
	w(ijk)=wdrop
	tn(ijk)=temdrop
	if(ifenrg.gt.0) rhoe(ijk)=f(ijk)*rhof*e1cal(0.0,tn(ijk),0.0,0.0)
 2000 continue
 1100 continue 
      return
      end