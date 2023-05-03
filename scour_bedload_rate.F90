!> \brief 
!>    Calculate bedload transport rate and depth for sediment species isspc
!
!>    References:
!> \verbatim
!>        1. Meyer-Peter, E. and Müller, R., 1948, Formulas for bed-load transport.
!>           Proceedings of the 2nd Meeting of the International Association for
!>           Hydraulic Structures Research. pp. 39–64.
!>        2. Soulsby, R., 1997, Dynamics of Marine Sands, Thomas Telford Publications,
!>           London, pp 160.
!>        3. Van Rijn, L. C., 1984, Sediment Transport, Part I: Bed load transport,
!>           Journal of Hydraulic Engineering 110(10), pp 1431-1456.
!> \endverbatim
!
!> \param[in]     dbed    Bedload transport depth (dimensional)
!> \param[in]     qd      Bed-load volume transport rate per unit width (volume/(time*width)) (dimensional)
!> \param[in]     dia     Dimensional grain size of species isspc
!> \param[in]     rhof    Fluid mass density
!> \param[in]     den     Material mass density of species isspc
!> \param[in]     gtot    Gravity acceleration (absolute value)
!> \param[in]     visc    Dynamic molecular viscosity of fluid
!> \param[in]     shields Shields number
!> \param[in]     crtshields Critical Shields number
!> \param[in]     isspc   Sediment species index
!
      subroutine scour_bedload_rate(dbed,qd,dia,rhof,den,gtot, & 
     &                              visc,shields,crtshields,isspc)
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
!      qdl:   bed-load volume transport rate per unit width (volume/(time*width)) (dimensionaless)
!      dstar: dimensionless grain size of species isspc
!      ibedload: flag for bedload transport equation, defined in scala.f
!               =1, Meyer-Peter & Müller equatin (1948)
!               =2, Nielsen equation (1992), from Soulsby (1997)
!               =3, Van Rijn equation (1984)
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
      include '../comdeck/params.f90'
      include '../comdeck/const.f90'
      include '../comdeck/scala.f90'
!
! calculate bedload transport depth: interpretation of saltation height (van Rijn (1984))
      dstar=dia*(rhof*max(zero,den-rhof)*gtot/visc**2)**ob3
      dbed=dia*0.3d0*dstar**0.7d0*sqrt(shields/crtshields-one)
!
      if(ibedload.eq.1) then
! dimensionless bed-load volume (net) transport rate from Meyer-Peter & Müller (1948)
        qdl=ssbdlf(isspc)*(shields-crtshields)**1.5d0
!
      else if(ibedload.eq.2) then
! dimensionless bed-load volume transport rate from Nielsen (1992)-in Soulsby (1997))
        qdl=ssbdlf2(isspc)*(shields-crtshields)*sqrt(shields)
!
      else if(ibedload.eq.3) then
! dimensionless bed-load volume transport rate from van Rijn (1984)
 
        if(dstar.gt.ztest) then
          qdl=ssbdlf3(isspc)*((shields/crtshields-one)**2.1d0)/(dstar**0.3d0)
        else
          qdl=zero
        endif
      end if
!
! dimensional bed-load volume transport rate per unit width (volume/(time*width))
      qd=qdl*sqrt((dia**3)*gtot*max(zero,den-rhof)/rhof)    ! net volume transport rate per unit width
!
      return
!
      end
