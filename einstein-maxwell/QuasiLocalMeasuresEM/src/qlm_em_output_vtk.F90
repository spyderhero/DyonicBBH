#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_output_vtk (CCTK_ARGUMENTS, hn, file_name)
  use cctk
  use constants
  use qlm_em_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer,      intent(in) :: hn
  character(*), intent(in) :: file_name

  integer, parameter :: unit = 194

  integer   :: nptheta, npphi
  integer   :: i, j
  CCTK_REAL :: xx, yy, zz

  nptheta = qlm_em_ntheta(hn) - 2*qlm_em_nghoststheta(hn) 
  npphi = qlm_em_nphi(hn) - 2*qlm_em_nghostsphi(hn)

  open (unit=unit, file=file_name, action='write')

  write (unit, '(A)') '# vtk DataFile Version 2.0' 
  write (unit, '(A)') 'Horizon data' 
  write (unit, '(A)') 'ASCII' 
  write (unit, '(A)') 'DATASET POLYDATA'
  write (unit, '(A,X,I10,X,A)') 'POINTS', nptheta*npphi, 'float'

  do i = 1+qlm_em_nghoststheta(hn), qlm_em_ntheta(hn)-qlm_em_nghoststheta(hn)
     do j = 1+qlm_em_nghostsphi(hn), qlm_em_nphi(hn)-qlm_em_nghostsphi(hn)
        xx = qlm_em_x(i, j, hn) 
        yy = qlm_em_y(i, j, hn) 
        zz = qlm_em_z(i, j, hn) 
        write (unit, *) xx, yy, zz 
     end do
  end do

  write (unit, '()')
  write (unit, '(A,X,I10,X,I10)') &
       'POLYGONS', npphi*(nptheta-1), 5*npphi*(nptheta-1)

  do i = 0, nptheta-2
     do j = 0, npphi-2
        write (unit,'(I1,4(I10))') &
             4, i*npphi+j, (i+1)*npphi+j, (i+1)*npphi+j+1, i*npphi+j+1
     end do
     write (unit,'(I1,4(I10))') &
          4, i*npphi+npphi-1, (i+1)*npphi+npphi-1, (i+1)*npphi, i*npphi
  end do

  write (unit, '()')
  write (unit, '(A,X,I10)') 'POINT_DATA', nptheta*npphi

  call writescalar ('shape', qlm_em_shape(:,:,hn))
  call writescalar ('l0', qlm_em_l0(:,:,hn))
  call writescalar ('l1', qlm_em_l1(:,:,hn))
  call writescalar ('l2', qlm_em_l2(:,:,hn))
  call writescalar ('l3', qlm_em_l3(:,:,hn))
  call writescalar ('n0', qlm_em_n0(:,:,hn))
  call writescalar ('n1', qlm_em_n1(:,:,hn))
  call writescalar ('n2', qlm_em_n2(:,:,hn))
  call writescalar ('n3', qlm_em_n3(:,:,hn))
  call writescalar_complex ('m0', qlm_em_m0(:,:,hn))
  call writescalar_complex ('m1', qlm_em_m1(:,:,hn))
  call writescalar_complex ('m2', qlm_em_m2(:,:,hn))
  call writescalar_complex ('m3', qlm_em_m3(:,:,hn))
  call writescalar_complex ('npkappa', qlm_em_npkappa(:,:,hn))
  call writescalar_complex ('nptau', qlm_em_nptau(:,:,hn))
  call writescalar_complex ('npsigma', qlm_em_npsigma(:,:,hn))
  call writescalar_complex ('nprho', qlm_em_nprho(:,:,hn))
  call writescalar_complex ('npepsilon', qlm_em_npepsilon(:,:,hn))
  call writescalar_complex ('npgamma', qlm_em_npgamma(:,:,hn))
  call writescalar_complex ('npbeta', qlm_em_npbeta(:,:,hn))
  call writescalar_complex ('npalpha', qlm_em_npalpha(:,:,hn))
  call writescalar_complex ('nppi', qlm_em_nppi(:,:,hn))
  call writescalar_complex ('npnu', qlm_em_npnu(:,:,hn))
  call writescalar_complex ('npmu', qlm_em_npmu(:,:,hn))
  call writescalar_complex ('nplambda', qlm_em_nplambda(:,:,hn))
  call writescalar_complex ('psi0', qlm_em_psi0(:,:,hn))
  call writescalar_complex ('psi1', qlm_em_psi1(:,:,hn))
  call writescalar_complex ('psi2', qlm_em_psi2(:,:,hn))
  call writescalar_complex ('psi3', qlm_em_psi3(:,:,hn))
  call writescalar_complex ('psi4', qlm_em_psi4(:,:,hn))
  call writescalar ('xit', qlm_em_xi_t(:,:,hn))
  call writescalar ('xip', qlm_em_xi_p(:,:,hn))
  call writescalar ('chi1', qlm_em_chi(:,:,hn))

  close (unit)

contains

  subroutine writescalar (array_name, array)
    character(*), intent(in) :: array_name
    CCTK_REAL,    intent(in) :: array(:, :)

    integer :: i, j

    write (unit, '(/A,X,A,X,A)') 'SCALARS', array_name, 'float 1'
    write (unit, '(A)') 'LOOKUP_TABLE default'
    do i = 1+qlm_em_nghoststheta(hn), qlm_em_ntheta(hn)-qlm_em_nghoststheta(hn)
       do j = 1+qlm_em_nghostsphi(hn), qlm_em_nphi(hn)-qlm_em_nghostsphi(hn)
          write (unit, *) array(i, j)
       end do
    end do
  end subroutine writescalar

  subroutine writescalar_complex (array_name, array)
    character(*), intent(in) :: array_name
    CCTK_COMPLEX, intent(in) :: array(:, :)

    call writescalar ('re' // array_name, real(array))
    call writescalar ('im' // array_name, aimag(array))
  end subroutine writescalar_complex

end subroutine qlm_em_output_vtk
