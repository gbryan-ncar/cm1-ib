  module ib_module
  implicit none

  public

    ! use immersed boundary technique?
!!!    logical, parameter :: do_ib  =  .true.


    ! variables and arrays for immersed boundary system:
    integer :: ibib,ieib,jbib,jeib,kbib,keib
    logical :: do_ib

  CONTAINS

!-----------------------------------------------------------------------

    subroutine ib_init(xh,yh,xf,yf,sigma,sigmaf,zs,zh,bndy,ibdy,jbdy,kbdy,xdamp)

    use input
    use constants

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      logical, intent(inout), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(inout), dimension(ibib:ieib,jbib:jeib) :: ibdy,jbdy,kbdy
      real, intent(inout), dimension(ibib:ieib,jbib:jeib) :: xdamp

      integer :: i,j,k,n,case
      real :: aa,hh,xc,xloc,yloc
      real, dimension(:,:), allocatable :: zsfoo

      ! Set case here:

      case = 2


      allocate( zsfoo(ib:ie,jb:je) )
      zsfoo = 0.0

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

      IF( case.eq.1 )THEN

        !  standard nh mountain wave case:

        hh =      400.0              ! max. height (m)
        aa =     1000.0              ! half width (m)
        xc =        0.0 + 0.5*dx     ! x-location (m)

        do j=jb,je
        do i=ib,ie
          zsfoo(i,j) = hh/( 1.0+( (xh(i)-xc)/aa )**2 )
        enddo
        enddo

      ENDIF

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

      IF( case.eq.2 )THEN

        !  2D block for tests of advection scheme:

        hh =      400.0              ! max. height (m)

        if( nx.eq.1 )then
          do j=jb,je
          do i=ib,ie
            if( abs(yh(j)-centery).le.400.0 ) zsfoo(i,j) = hh
          enddo
          enddo
        else
          do j=jb,je
          do i=ib,ie
            if( abs(xh(i)-centerx).le.400.0 ) zsfoo(i,j) = hh
          enddo
          enddo
        endif

      ENDIF

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

      IF( case.eq.3 )THEN

        !  standard test case:  cube in a wind tunnel

      do n=1,3
        if( n.eq.1 ) yloc = 0.5*maxy
        if( n.eq.2 ) yloc = 0.5*maxy + 3.0*100.0
        if( n.eq.3 ) yloc = 0.5*maxy - 3.0*100.0
        print *,'  n,yloc = ',n,yloc
        do j=jb,je
        do i=ib,ie
          if( abs(xh(i)-0.25*maxx).le.50.0 .and. abs(yh(j)-yloc).le.50.0 ) zsfoo(i,j) = 100.0
        enddo
        enddo
      enddo

      ENDIF

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

      do j=jb,je
      do i=ib,ie
        if( abs(xh(i)-0.95*maxx).le.0.05*maxx )then
          xdamp(i,j) = sin( 0.5*pi*(xh(i)-0.90*maxx)/(0.05*maxx) )
        endif
      enddo
      enddo

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

    !---  DO NOT CHANGE ANYTHING BELOW HERE  ---!

      do k=1,nk
      do j=jb,je
      do i=ib,ie
        if( sigma(k).le.zsfoo(i,j) ) bndy(i,j,k) = .true.
      enddo
      enddo
      enddo

      ! last step:  get ibdy,jbdy,kbdy

      do k=1,nk
      do j=jb+1,je-1
      do i=ib+1,ie-1
        if( bndy(i,j,k) )then
          kbdy(i,j) = k+1
          if( .not. bndy(i-1,j,k) ) ibdy(i,j) = 1  ! west facing wall
          if( .not. bndy(i+1,j,k) ) ibdy(i,j) = 2  ! east facing wall
          if( ( .not. bndy(i-1,j,k) ) .and. ( .not. bndy(i+1,j,k) ) ) ibdy(i,j) = 3  ! west and east facing wall
          if( .not. bndy(i,j-1,k) ) jbdy(i,j) = 1  ! south facing wall
          if( .not. bndy(i,j+1,k) ) jbdy(i,j) = 2  ! north facing wall
          if( ( .not. bndy(i,j-1,k) ) .and. ( .not. bndy(i,j+1,k) ) ) jbdy(i,j) = 3  ! south and north facing wall
        endif
      enddo
      enddo
      enddo

      deallocate( zsfoo )

    end subroutine ib_init

!-----------------------------------------------------------------------

  end module ib_module
