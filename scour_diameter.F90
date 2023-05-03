!> \brief 
!>    Calculate average d50 in cell using sclrn if iopt=1 and sclr otherwise
!
      subroutine scour_diameter(ijk,d50,iopt)
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
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/state.f90'
      include '../comdeck/scala.f90'
      include '../comdeck/const.f90'
!
! Initialize the D50 particle size and total packed sediment concentration
!
      d50=zero
!
      scr1t=zero
      if(nsspc.gt.1) then
        do isspc=1,nsspc
          iscr1=isspk(isspc)
          if(iopt.eq.1) then
            sclr_packed=sclrn(ijk,iscr1)
          else
            sclr_packed=sclr(ijk,iscr1)
          endif
          d50=d50+ssdia(isspc)*sclr_packed
          scr1t=scr1t+sclr_packed
        enddo
!
! Compute the D50 particle size; if no packed sediment present, get out to avoid
! divide by zero error
!
        if(scr1t.gt.ztest) then
          d50=d50/scr1t
        else
          return
        endif
      else
        d50=ssdia(1)
      endif
!
      return
      end
!
