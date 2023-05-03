      subroutine heatfl(i,j,k,heatfl1,heatfl2)
!c
!c  add heat sources to fluid
!c
!c    ***************************************************************
!c    **                          notice                           **
!c    **  this subprogram contains flow science, inc. proprietary  **
!c    **       trade secret and confidential information.          **
!c    **                                                           **
!c    **                unauthorized use prohibited                **
!c    **           copyright 1985-2012 flow science, inc.          **
!c    ***************************************************************
!c
      use arrays_module
!c
      use meshcb_module
!c
      use obsijk_module
!c
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/const.f90'
      include '../comdeck/diag.f90'
      include '../comdeck/logvar.f90'
      include '../comdeck/obsd.f90'
      include '../comdeck/bcinfo.f90'
      include '../comdeck/acctrm.f90'
      include '../comdeck/scala.f90'
      include '../comdeck/state.f90'
      include '../comdeck/dumn.f90'
      include '../comdeck/ijk.f90'
!!c
!c  i,j,k         = current cell indices
!c  heatfl1       = heat source for fluid 1
!c  heatfl2       = heat source for fluid 2
!c
!c  heatfl1,2 is thermal energy per unit mass of fluid material
!c  per unit time (or power per unit mass)
!c
!c        add heat source           
        apx=dum10  !c 0 0.063 0.001
        apy=dum11
        apz=dum12
        pi=3.141593
	Qm=6*3**(1.0/2.0)*dum1*dum2*dum5/(dum8*dum9*pi**(3.0/2.0))!c������Դ  �����ܶ�
!c skip if energy trnsport is not turned on
!c
      heatfl1=zero
      heatfl2=zero
      if(ifenrg.eq.0) goto 200  !c   Flag for internal energy evaluation
!c
!c
!c
!c ------   skip passive cells
!c
       if(icstat(ijk).le.0) goto 200  !c
!c
!c ------   skip non-surface cells
!c
       if(f(ijk).lt.emf .or. fn(ijk).lt.emf) goto 200
!c ------  look for phantom obstacle of type ifob=2 in cell
!c
!c     
       !c һ������
             if((xi(i)-dum13*t).ge.apx)  then
         	    argexp=(xi(i)-apx-dum13*t)**2/dum7**2+(yj(j)-apy)**2/dum8**2+(zk(k)-apz)**2/dum9**2 !c��Դ��Χ  ��̬λ��  apz=0.001
         	    if(argexp.lt.3.0) then  !c�޶���Դ��Χ
         		    bamp1=Qm*1.4*exp(-3*argexp)/dum7
         		    vph=10**(dum20+6.121-18836/tn(ijk)-0.5*log10(tn(ijk)))
         		    bamp2=dum19*vph
         		    bamp=bamp1-bamp2
         	    endif
             else
        	    argexp=(xi(i)-apx-dum13*t)**2/dum6**2+(yj(j)-apy)**2/dum8**2+(zk(k)-apz)**2/dum9**2
                    if(argexp.lt.3.0) then
         		    bamp1=Qm*0.6*exp(-3*argexp)/dum6
         		    vph=10**(dum20+6.121-18836/tn(ijk)-0.5*log10(tn(ijk)))
         		    bamp2=dum19*vph
         		    bamp=bamp1-bamp2
                    endif
             endif
       !c��������  xi��ʱ��仯 ��һ��Ϊxi(i)-apx-dum13*t �ڶ���Ӧ��Ϊxi(i)-apx-dum13*��t-40��zk(k)-apz-0.06
       !if(t.lt.40)then
       !      if((xi(i)-dum13*t).ge.apx)  then
       !  	    argexp=(xi(i)-apx-dum13*t)**2/dum7**2+(yj(j)-apy)**2/dum8**2+(zk(k)-apz)**2/dum9**2 !c��Դ��Χ  ��̬λ��  apz=0.001
       !  	    if(argexp.lt.3.0) then  !c�޶���Դ��Χ
       !  		    bamp1=Qm*1.4*exp(-3*argexp)/dum7
       !  		    vph=10**(dum20+6.121-18836/tn(ijk)-0.5*log10(tn(ijk)))
       !  		    bamp2=dum19*vph
       !  		    bamp=bamp1-bamp2
       !  	    endif
       !      else
       ! 	    argexp=(xi(i)-apx-dum13*t)**2/dum6**2+(yj(j)-apy)**2/dum8**2+(zk(k)-apz)**2/dum9**2
       !             if(argexp.lt.3.0) then
       !  		    bamp1=Qm*0.6*exp(-3*argexp)/dum6
       !  		    vph=10**(dum20+6.121-18836/tn(ijk)-0.5*log10(tn(ijk)))
       !  		    bamp2=dum19*vph
       !  		    bamp=bamp1-bamp2
	      !  	    endif
       !      endif
       !else
       !    if((xi(i)-dum13*t).ge.apx)  then
       ! 	    argexp=(xi(i)-apx-dum13*t)**2/dum7**2+(yj(j)-apy)**2/dum8**2+(zk(k)-apz-0.06)**2/dum9**2 !c��Դ��Χ  ��̬λ��
       !  	    if(argexp.lt.3.0) then  !c�޶���Դ��Χ
       ! 		    bamp1=Qm*1.4*exp(-3*argexp)/dum7
       !  		    vph=10**(dum20+6.121-18836/tn(ijk)-0.5*log10(tn(ijk)))
       !  		    bamp2=dum19*vph
       !  		    bamp=bamp1-bamp2
       !  	    endif
       !      else
       ! 	    argexp=(xi(i)-apx-dum13*t)**2/dum6**2+(yj(j)-apy)**2/dum8**2+(zk(k)-apz-0.06)**2/dum9**2
       !             if(argexp.lt.3.0) then
       !  		    bamp1=Qm*0.6*exp(-3*argexp)/dum6
       !  		    vph=10**(dum20+6.121-18836/tn(ijk)-0.5*log10(tn(ijk)))
       !  		    bamp2=dum19*vph
       !  		    bamp=bamp1-bamp2
	      !  	    endif
       !      endif
       !endif
!c
!c
      heatfl1=bamp/rho(ijk)  !c       rho(ijk)       density in cell (only for variable density)
      heatfl2=zero
!c
!c
!ccc     include '../comdeck/ijk.f'
!c
  200   continue
      return
      end