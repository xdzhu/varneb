! Diagnostic wrapper around a user's already licensed VASP lattlib.o.
! Contains no VASP lattice-classification implementation. Never runs DFT,
! modifies VASP, changes a structure, or certifies an electronic result.
! Intel classic object/module ABI is required by this particular harness.
module diagnostic_tolerance
  use iso_c_binding, only: c_double
  implicit none
  real(c_double), bind(C, name="sym_prec_mp_tiny_") :: tolerance
end module

program lattice_probe
  use diagnostic_tolerance
  use iso_c_binding, only: c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(c_double) :: cell(3,3), reciprocal(3,3), temporary(3,3), dimensions(6), volume
  integer :: axis, real_type, reciprocal_type, expected, ibrava2b
  external :: lattyp, ibrava2b
  read(*,*) tolerance
  do axis=1,3
    read(*,*) cell(:,axis)
  end do
  if (.not. all(ieee_is_finite(cell)) .or. .not. ieee_is_finite(tolerance)) error stop 2
  if (tolerance <= 0) error stop 2
  reciprocal(:,1)=cross_product(cell(:,2),cell(:,3))
  reciprocal(:,2)=cross_product(cell(:,3),cell(:,1))
  reciprocal(:,3)=cross_product(cell(:,1),cell(:,2))
  volume=dot_product(cell(:,1),reciprocal(:,1))
  if (volume <= 0) error stop 2
  reciprocal=reciprocal/volume
  temporary=cell
  call lattyp(temporary(:,1),temporary(:,2),temporary(:,3),real_type,dimensions,-1)
  temporary=reciprocal
  call lattyp(temporary(:,1),temporary(:,2),temporary(:,3),reciprocal_type,dimensions,-1)
  expected=ibrava2b(real_type)
  write(*,'(3(I0,1X))') real_type,reciprocal_type,expected
contains
  function cross_product(a,b) result(c)
    real(c_double), intent(in) :: a(3),b(3)
    real(c_double) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function
end program

! Fail-closed dependencies for unused/error paths in the monolithic object.
! LATTYP uses its object's CELVOL; it must never call the k-point HNFORM.
subroutine diagnostic_unused_hnform() bind(C,name="mkpoints_mp_hnform_")
  error stop "Unexpected HNFORM call: diagnostic harness is not valid for this path"
end subroutine

subroutine error(routine,message,index)
  implicit none
  character(*), intent(in) :: routine,message
  integer, intent(in) :: index
  write(*,*) "Licensed lattice routine rejected input: ",routine,message,index
  error stop 3
end subroutine
