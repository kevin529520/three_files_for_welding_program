      subroutine mvbfrc2_usr(fx,fy,fz,nob,m)
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
!    Customization subroutine for m-th control force on GMO acting at fixed points
!    defined in prepin file: xfpnt, yfpnt and zfpnt.
!
!    It is called if igf(nob)=2 and imo(nob)=4 are both satisfied.
!    (igf(nob)=2 means to provide up to 3 forces with application points on the GMO,
!    namely, m can be 1, 2 or 3. imo(nob)=4 means coupled motion,
!    nob is the index of the geometry component)
!
!    *** Note, this routine is called from outside over all moving objects and m values
!    (m=1,2,3).
!    If the m-th force is not the one for customization or nob is not the
!    object for customization, then use "return" as shown in
!    the example. You can also use "if ... then" for the object index and force
!    index instead of using "if ... return", e.g., if(m.eq.2.and.nob.eq.3) then ...end if"
!    Once this subroutine is called, the force will be automatically transformed
!    into space system if it is orginally defined in body system.
!
!
!                     Arguments in this subroutine
!
!    fx         x component of the m-th control force on the GMO. The force is
!               in the reference system defined by irefgfm(nob,m).
!               Namely, if irefgfm(nob,m)=1, then the force is in space system.
!               If irefgfm(nob,m)=2, then the force is in body system.
!               irefgfm(nob,m) is defined in the "Control Force and Torques" tab
!               for the GMO setup in the GUI. The same is true for fy and fz.
!    fy         y component of the m-th control force on mass center in
!               the reference system defined by irefgfm(nob,m).
!    fz         z component of the m-th control force on mass center in
!               the reference system defined by irefgfm(nob,m).
!
!   nob         object index
!
!    m          index of force. If igf(nob)=2 and nob has coupled 6-DOF motion (imo(nob)=4),
!               then user is allowed to define up to three independent forces with
!               their application points on the GMO. Here m is the index of the force.
!
!
!                     Useful variables available in this subroutine
!
!   t:  current time
!
!   nobs :    total number of geometry components including both GMO and
!             non-GMO objects.
!
!   mvbmx:    total number of moving objects (imo(nob)=3 and 4).
!
!   imo(nob): type of object.
!             =0,  non-GMO object
!             =3,  moving object with prescribed motion
!             =4,  moving object with coupled motion
!
!   imoaxs(nob): flag for motion constraints on an moving object
!                =0 (default) object has 6-DOF motion
!                =1 fixed-axis rotation, axis is parallel to x axis of space system
!                =2 fixed-axis rotation, axis is parallel to y axis of space system
!                =3 fixed-axis rotation, axis is parallel to z axis of space system
!                =4 fixed-point rotation
!
!   xaxs(nob): x coordinate of fixed rotation-axis/point in space system (imo(nob)=3 or 4)
!   yaxs(nob): y coordinate of fixed rotation-axis/point in space system (imo(nob)=3 or 4)
!   zaxs(nob): z coordinate of fixed rotation-axis/point in space system (imo(nob)=3 or 4)
!
!   massmvb(nob): mass of moving object
!
!   xmcmvb0(nob)
!   ymcmvb0(nob)
!   zmcmvb0(nob): mass center's initial coordinates in space system (at t=0).
!
!   xmcmvb(nob)
!   ymcmvb(nob)
!   zmcmvb(nob): mass center's coordinates in space system at current time t.
!
!   tjmvb(n,i,j): moment of inertia tensor [J] in body system.
!                 i is row index, j is collum index.
!                 See FLOW-3D manual for mathematical definition.
!                 1) For motion with 6-DOF, it is about the mass center.
!                 2) For motion about a fixed point, [J] is about mass center in input,
!                 but here it is already tramsformed into moment of inertia tensor about
!                 the fixed point.
!
!   tjaxs(nob): Moment of inertia about fixed axis (for imoaxis(nob)=1 or 2 or 3 only)
!
!   umc(nob),
!   vmc(nob),
!   wmc(nob): x, y and z components of mass center velocity in space system.
!
!   omxmvb(nob)
!   omymvb(nob)
!   omzmvb(nob): angular velocity components in x, y and z of body system
!
!   omxmvb_s(nob)
!   omymvb_s(nob)
!   omzmvb_s(nob): angular velocity components in x, y and z of space system,
!
!   anglmvb(nob): rotation angle measured from the initial orientation. For
!                 rotation with fixed axis only.
!                 E.g., if the rotation axis is parallel to x, the the positive
!                 rotation direction is in +x based on right hand convention.
!
!   xmcmvbn(nob)
!   ymcmvbn(nob)
!   zmcmvbn(nob)
!   umcn(nob)
!   vmcn(nob)
!   wmcn(nob)
!   omxmvbn(nob)
!   omymvbn(nob)
!   omzmvbn(nob): variables at previous time step
!
!   xref0(nob),
!   yref0(nob),
!   zref0(nob): initial coordinates in space system for the reference point fixed on
!              moving object nob. For coupled motion:
!              1) 6-DOF coupled motion, it is mass center.
!              2) fixed axis motion, it is a point on the axis.
!              3) fixed point motion, it is the fixed point.
!
!   xref(nob),
!   yref(nob),
!   zref(nob): location of the reference point in space system at time t.
!
!   trfmvb(nob,i,j):  coordinate transformation tensor for moving object nob.
!                   i is row index, j is collum index. It is the matrix [R] in tech note.
!                   [R] is an orthogonal matrix: its inverse and its transpose are equal.
!
!                   1) To transform a vector a(i) (collum vector, i=1,2,3) from body-fixed system to
!                      space system, use [as]=[R]*[ab], where [ab] and [as] are in
!                      body-fixed system and space system respectively.
!
!                        do i=1,3
!                           as(i)=0.0
!                           do j=1,3
!                             as(i)=as(i)+trfmvb(nob,i,j)*ab(j)
!                           end do
!                        end do
!
!                   2) To transform a vector a(i) from space system to body-fixed system,
!                      [ab]=[R]'*[as], [R]' is transpose/inverse of [R]
!
!                        do i=1,3
!                           ab(i)=0.0
!                           do j=1,3
!                             ab(i)=ab(i)+trfmvb(nob,j,i)*as(j)
!                           end do
!                        end do
!
!                   3) To find location of a point fixed on GMO at time t.
!                      Assume its initial location in space system is (x0,y0,z0).
!                      Its current location (xpt,ypt,zpt) in space system can be
!                      obtained in the following way.
!
!                     xr0=x0-xref0(nob)
!                     yr0=y0-yref0(nob)
!                     zr0=z0-zref0(nob)
!                     xpt=trfmvb(nob,1,1)*xr0+trfmvb(nob,1,2)*yr0+
!          &               trfmvb(nob,1,3)*zr0+xref(nob)
!                     ypt=trfmvb(nob,2,1)*xr0+trfmvb(nob,2,2)*yr0+
!          &                trfmvb(nob,2,3)*zr0+yref(nob)
!                     zpt=trfmvb(nob,3,1)*xr0+trfmvb(nob,3,2)*yr0+
!          &                trfmvb(nob,3,3)*zr0+zref(nob)
!
!
      use mblock_module
!
      use arrays_module
!
      use meshcb_module
!
      use mblock_module
!
      use obsijk_module
!
#ifdef SINGLE
      include '../comdeck/precis4.f90'
#else
      include '../comdeck/precis.f90'
#endif
!
      include '../comdeck/params.f90'
      include '../comdeck/const.f90'
      include '../comdeck/acctrm.f90'
      include '../comdeck/mparam.f90'
      include '../comdeck/cntrl.f90'
      include '../comdeck/dumn.f90'
!
!
!cccccccccccccccccccc Example cccccccccccccc
!
!     Example: to customize the 2nd control force in x direction for object #3.
!     Assume it decreases exponentially with time. The customization is
!
!     nob_usr=3
!     m_usr=2
!     if(nob.ne.nob_usr .or. m.ne.m_usr) return
!
!     frc0=1.342e5      ! initial value of force
!     c0=0.1            ! exponential coefficient
!     fx=frc0*exp(-c0*t)  ! force at current time t
!
!     Note: in this example, the y and z components of the m-th control force on the GMO
!     are not affected by the customization. All other control forces/torques are not
!     affect either. They are calculated in other subroutines based on the input data.
!cccccccccccccccccccc end of example  cccccccccccccccccc
!
      return
      end
