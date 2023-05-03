      subroutine mvbvel_usr(nob)
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
!    Customization subroutine for velocity of object nob. The object must be a GMO under
!    prescribed motion (imo(nob)=3). The motion can be prescribed 6-DOF motion, or it is
!    a fixed-axis or fixed-point rotation. Argument nob is object (geometry component) index.
!
!    *** Note, this routine is called from outside over all moving objects. If nob is
!    not the object for customization, use "return" as in the example.
!
!    The following existing variables can be defined by user in this subroutine:
!
!    uref(nob): x component of velocity of the reference point for object nob.
!
!    vref(nob): y component of velocity of the reference point for object nob.
!
!    wref(nob): z component of velocity of the reference point for object nob.
!
!    omxmvb(nob): x component of angular velocity (rad/s) in body-fixed system for object nob.
!
!    omymvb(nob): y component of angular velocity (rad/s) in body-fixed system for object nob.
!
!    omzmvb(nob): z component of angular velocity (rad/s) in body-fixed system for object nob.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!               useful variables available in this subroutine
!
!   t:  current time
!
!   xref0(nob),
!   yref0(nob),
!   zref0(nob): reference point's initial coordinates in space system for object nob.
!               The reference point is fixed on the GMO.
!              1) For prescribed 6-DOF motion: they are input variables.
!                 By default, the reference point is the object's mass center
!                 if mass properties (density,or mass center location) are given.
!                 If the mass properties are not given, by default the reference point is
!                 the object's geometry center.
!              2) For fixed-axis motion, the reference point is a point on the axis
!              3) For fixed-point motion, the refernce is the fixed point.
!
!   xref(nob),
!   yref(nob),
!   zref(nob): reference point location at time t. The reference point is fixed on object nob.
!
!   imo(nob): type of object.
!             =0,  non-GMO object
!             =3,  moving object with prescribed motion
!             =4,  moving object with coupled motion
!
!   imoaxs(nob): flag for motion constraints
!                =0 (default) 6-DOF motion
!                =1 fixed-axis rotation, the axis is parallel to x axis of space system
!                =2 fixed-axis rotation, the axis is parallel to y axis of space system
!                =3 fixed-axis rotation, the axis is parallel to z axis of space system
!                =4 fixed-point rotation
!
!   xaxs(nob): x coordinate of fixed rotation-axis/point in space system
!   yaxs(nob): y coordinate of fixed rotation-axis/point in space system
!   zaxs(nob): z coordinate of fixed rotation-axis/point in space system
!
!   umc(nob),
!   vmc(nob),
!   wmc(nob): x, y and z components of mass center velocity in space system at time t.
!
!   anglmvb(nob): rotation angle measured from the initial orientation for
!                 moving object with fixed axis. E.g., if the rotation axis is
!                 parallel to x, the the positive rotation angle is in +x
!                 based on right hand convention.
!
!   anglmin(nob): maximum angle (absolute value, in degree) allowed in the
!                 negative rotational direction, default is infinitive.
!
!   anglmax(nob): maximum angle (absolute value, in degree) allowed in the
!                 positive rotatinal direction, default is infinitive.
!
!   xrefn(nob)
!   yrefn(nob)
!   zrefn(nob)
!   urefn(nob)
!   vrefn(nob)
!   wrefn(nob)
!   omxmvbn(nob)
!   omymvbn(nob)
!   omzmvbn(nob): variables at previous time step
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
!     Example: to define x and y-velocities of reference point for object #3 and #4.
!     Both objects are under prescribed 6-DOF motion. Assume the x-velocity decreases
!     exponentially with t, and y-velocity is a linear function of the object's
!     displacement in y. The customization is:
!
!     if(nob.ne.3 .and. nob.ne.4) return
!
!     vx0=5.0      ! initial x velocity of the reference point
!     c0=0.1            ! exponential coefficient
!     uref(nob)=vx0*exp(-c0*t)  ! x velocity at current time t
!
!     c1=2.5       ! linear coefficient
!     c2=20.0
!     dyref=yref(nob)-yref0(nob)   ! displacement of reference point in y
!     vref(nob)=c1*dyref+c2   ! velocity in y of reference point
!
!cccccccccccccccccccc end of example  cccccccccccccccccc
!
      return
      end
