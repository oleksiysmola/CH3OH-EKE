   recursive subroutine ML_symmetry_transformation_XY3_IV(nsym,src,dst,ndeg)
    implicit none 
    !
    !integer, parameter :: ark         = selected_real_kind(12,25)  ! "Accurate" reals
    integer,intent(in)    :: nsym  ! n of irreps
    integer,intent(in)    :: ndeg  ! n degrees of freedom 
    double precision,intent(in)      :: src(1:ndeg)
    double precision,intent(out)     :: dst(1:ndeg,nsym)
    !
    integer :: ioper
    !
    double precision         :: repres(nsym,ndeg,ndeg),a,b,e,o,pi
    !
    a = 0.5d0 ; b = 0.5d0*sqrt(3.0d0) ; e = 1.0d0 ; o = 0.0d0
    pi = 4.0d0 * datan2(1.0d0,1.0d0)
    !
    if (verbose>=6) print*,"ML_symmetry_transformation_XY3_IV"
    !
    if (nsym>6) then
      write (6,"('symmetry_transformation_local: illegal nsym = ',i8)") nsym
      stop 'symmetry_transformation_local: illegal nsym'
    endif
    !
    repres = 0 
    !
    ! repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    !
    repres(:,12,12) = 1.0_ark
    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    !d9
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    !
    !C3+/(123)
    repres(2,3,5) = 1.0_ark
    repres(2,4,3) = 1.0_ark
    repres(2,5,4) = 1.0_ark
    !
    repres(2,7,9) = 1.0_ark
    repres(2,8,7) = 1.0_ark
    repres(2,9,8) = 1.0_ark
    !
    repres(2,10,10) = -a
    repres(2,10,11) = -b
    repres(2,11,10) =  b
    repres(2,11,11) = -a
    !
    !C3-/(132)
    !
    repres(3,3,4) = 1.0_ark
    repres(3,4,5) = 1.0_ark
    repres(3,5,3) = 1.0_ark
    !
    repres(3,7,8) = 1.0_ark
    repres(3,8,9) = 1.0_ark
    repres(3,9,7) = 1.0_ark
    !
    repres(3,10,10) = -a
    repres(3,10,11) =  b
    repres(3,11,10) = -b
    repres(3,11,11) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,7,7) = 1.0_ark
    repres(4,8,9) = 1.0_ark
    repres(4,9,8) = 1.0_ark
    !
    repres(4,10,10) =  1.0_ark
    repres(4,11,11) = -1.0_ark
    !
    !C2'/(9)->(12)*
    repres(6,3,4) = 1.0_ark
    repres(6,4,3) = 1.0_ark
    repres(6,5,5) = 1.0_ark
    !
    repres(6,7,8)  = 1.0_ark
    repres(6,8,7)  = 1.0_ark
    repres(6,9,9)  = 1.0_ark
    !
    repres(6,10,10) = -a
    repres(6,10,11) =  b
    repres(6,11,10) =  b
    repres(6,11,11) =  a
    !
    !(13)->(35)
    repres(5,3,5) = 1.0_ark
    repres(5,4,4) = 1.0_ark
    repres(5,5,3) = 1.0_ark
    !
    repres(5,7,9) = 1.0_ark
    repres(5,8,8) = 1.0_ark
    repres(5,9,7) = 1.0_ark
    !
    repres(5,10,10) = -a
    repres(5,10,11) = -b
    repres(5,11,10) = -b
    repres(5,11,11) =  a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo
    !
    dst(12,1) = src(12)
    dst(12,2) = src(12)+2.0d0*pi/3.0d0
    dst(12,3) = src(12)-2.0d0*pi/3.0d0
    dst(12,4) =-src(12)
    dst(12,5) =-src(12)-2.0d0*pi/3.0d0
    dst(12,6) =-src(12)+2.0d0*pi/3.0d0
    !
  end subroutine ML_symmetry_transformation_XY3_IV 