#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



! TODO:
! instead of interpolating to the symmetry points, copy them



! A convenient shortcut
#define P(x) CCTK_PointerTo(x)



subroutine qlm_em_interpolate (CCTK_ARGUMENTS, hn)
  use cctk
  use qlm_em_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_INT,  parameter :: izero = 0
  integer,   parameter :: ik = kind(izero)
  integer,   parameter :: sk = kind(interpolator)
  CCTK_REAL, parameter :: one = 1
  
  CCTK_REAL, parameter :: poison_value = -42
  
  integer      :: len_coordsystem
  integer      :: len_interpolator
  integer      :: len_interpolator_options
  
  character    :: fort_coordsystem*100
  character    :: fort_interpolator*100
  character    :: fort_interpolator_options*1000
  
  integer      :: nvars
  
  integer      :: coord_handle
  integer      :: interp_handle
  integer      :: options_table
  
  integer      :: ninputs
  integer      :: noutputs
  
  CCTK_REAL, allocatable :: xcoord(:,:)
  CCTK_REAL, allocatable :: ycoord(:,:)
  CCTK_REAL, allocatable :: zcoord(:,:)
  
  integer      :: ind_gxx, ind_gxy, ind_gxz, ind_gyy, ind_gyz, ind_gzz
  integer      :: ind_kxx, ind_kxy, ind_kxz, ind_kyy, ind_kyz, ind_kzz
  integer      :: ind_alpha
  integer      :: ind_betax, ind_betay, ind_betaz
  integer      :: ind_ttt
  integer      :: ind_ttx, ind_tty, ind_ttz
  integer      :: ind_txx, ind_txy, ind_txz, ind_tyy, ind_tyz, ind_tzz
  integer      :: ind_ex, ind_ey, ind_ez
  integer      :: ind_ax, ind_ay, ind_az
  integer      :: ind_phi1, ind_phi2
  
  integer      :: coord_type
  CCTK_POINTER :: coords(3)
  CCTK_INT     :: inputs(34)
  CCTK_INT     :: output_types(115)
  CCTK_POINTER :: outputs(115)
  CCTK_INT     :: operand_indices(115)
  CCTK_INT     :: operation_codes(115)
  integer      :: npoints
  
  character    :: msg*1000
  
  integer      :: ni, nj
  
  integer      :: ierr
  
  
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Interpolating 3d grid functions")
  end if
  
  
  
  if (shift_state==0) then
     call CCTK_WARN (0, "The shift must have storage")
  end if
  
!!$  if (stress_energy_state==0) then
!!$     call CCTK_WARN (0, "The stress-energy tensor must have storage")
!!$  end if
  
  
  
  ! Get coordinate system
  call CCTK_FortranString &
       (len_coordsystem, int(coordsystem,sk), fort_coordsystem)
  call CCTK_CoordSystemHandle (coord_handle, fort_coordsystem)
  if (coord_handle<0) then
     write (msg, '("The coordinate system """, a, """ does not exist")') &
          trim(fort_coordsystem)
     call CCTK_WARN (0, msg)
  end if
  
  ! Get interpolator
  call CCTK_FortranString &
       (len_interpolator, int(interpolator,sk), fort_interpolator)
  call CCTK_InterpHandle (interp_handle, fort_interpolator)
  if (interp_handle<0) then
     write (msg, '("The interpolator """,a,""" does not exist")') &
          trim(fort_interpolator)
     call CCTK_WARN (0, msg)
  end if
  
  ! Get interpolator options
  call CCTK_FortranString &
       (len_interpolator_options, int(interpolator_options,sk), &
       fort_interpolator_options)
  call Util_TableCreateFromString (options_table, fort_interpolator_options)
  if (options_table<0) then
     write (msg, '("The interpolator_options """,a,""" have a wrong syntax")') &
          trim(fort_interpolator_options)
     call CCTK_WARN (0, msg)
  end if
  
  
  
  if (hn > 0) then
     
     ni = qlm_em_ntheta(hn)
     nj = qlm_em_nphi(hn)
     
     allocate (xcoord(ni,nj))
     allocate (ycoord(ni,nj))
     allocate (zcoord(ni,nj))
     
     xcoord(:,:) = qlm_em_x(:ni,:nj,hn)
     ycoord(:,:) = qlm_em_y(:ni,:nj,hn)
     zcoord(:,:) = qlm_em_z(:ni,:nj,hn)
     
  end if
  
  
  
  ! TODO: check the excision mask
  
  ! Get variable indices
  call CCTK_VarIndex (ind_gxx  , "ADMBase::gxx"   )
  call CCTK_VarIndex (ind_gxy  , "ADMBase::gxy"   )
  call CCTK_VarIndex (ind_gxz  , "ADMBase::gxz"   )
  call CCTK_VarIndex (ind_gyy  , "ADMBase::gyy"   )
  call CCTK_VarIndex (ind_gyz  , "ADMBase::gyz"   )
  call CCTK_VarIndex (ind_gzz  , "ADMBase::gzz"   )
  call CCTK_VarIndex (ind_kxx  , "ADMBase::kxx"   )
  call CCTK_VarIndex (ind_kxy  , "ADMBase::kxy"   )
  call CCTK_VarIndex (ind_kxz  , "ADMBase::kxz"   )
  call CCTK_VarIndex (ind_kyy  , "ADMBase::kyy"   )
  call CCTK_VarIndex (ind_kyz  , "ADMBase::kyz"   )
  call CCTK_VarIndex (ind_kzz  , "ADMBase::kzz"   )
  call CCTK_VarIndex (ind_alpha, "ADMBase::alp"   )
  call CCTK_VarIndex (ind_betax, "ADMBase::betax" )
  call CCTK_VarIndex (ind_betay, "ADMBase::betay" )
  call CCTK_VarIndex (ind_betaz, "ADMBase::betaz" )
  if (stress_energy_state /= 0) then
     call CCTK_VarIndex (ind_ttt  , "TmunuBase::eTtt")
     call CCTK_VarIndex (ind_ttx  , "TmunuBase::eTtx")
     call CCTK_VarIndex (ind_tty  , "TmunuBase::eTty")
     call CCTK_VarIndex (ind_ttz  , "TmunuBase::eTtz")
     call CCTK_VarIndex (ind_txx  , "TmunuBase::eTxx")
     call CCTK_VarIndex (ind_txy  , "TmunuBase::eTxy")
     call CCTK_VarIndex (ind_txz  , "TmunuBase::eTxz")
     call CCTK_VarIndex (ind_tyy  , "TmunuBase::eTyy")
     call CCTK_VarIndex (ind_tyz  , "TmunuBase::eTyz")
     call CCTK_VarIndex (ind_tzz  , "TmunuBase::eTzz")
  else
     ind_ttt = -1
     ind_ttx = -1
     ind_tty = -1
     ind_ttz = -1
     ind_txx = -1
     ind_txy = -1
     ind_txz = -1
     ind_tyy = -1
     ind_tyz = -1
     ind_tzz = -1
  end if
  if (calc_charge /= 0) then
     call CCTK_VarIndex (ind_ex  , "ProcaBase::Ex")
     call CCTK_VarIndex (ind_ey  , "ProcaBase::Ey")
     call CCTK_VarIndex (ind_ez  , "ProcaBase::Ez")
     call CCTK_VarIndex (ind_ax  , "ProcaBase::Ax")
     call CCTK_VarIndex (ind_ay  , "ProcaBase::Ay")
     call CCTK_VarIndex (ind_az  , "ProcaBase::Az")
  else
     ind_ex = -1
     ind_ey = -1
     ind_ez = -1
     ind_ax = -1
     ind_ay = -1
     ind_az = -1
  end if
  call CCTK_VarIndex (ind_phi1, "ScalarBase::phi1")
  call CCTK_VarIndex (ind_phi2, "ScalarBase::phi2")
  
  
  
  ! Set up the interpolator arguments
  coord_type = CCTK_VARIABLE_REAL
  if (hn > 0) then
     npoints = ni * nj
     coords(:) = (/ P(xcoord), P(ycoord), P(zcoord) /)
  else
     npoints = 0
     coords(:) = CCTK_NullPointer()
  end if
  
  inputs = (/ &
       ind_gxx, ind_gxy, ind_gxz, ind_gyy, ind_gyz, ind_gzz, &
       ind_kxx, ind_kxy, ind_kxz, ind_kyy, ind_kyz, ind_kzz, &
       ind_alpha, &
       ind_betax, ind_betay, ind_betaz, &
       ind_ttt, &
       ind_ttx, ind_tty, ind_ttz, &
       ind_txx, ind_txy, ind_txz, ind_tyy, ind_tyz, ind_tzz, &
       ind_ex, ind_ey, ind_ez, &
       ind_ax, ind_ay, ind_az, &
       ind_phi1, ind_phi2 /)
  
  call CCTK_NumVars (nvars)
  if (nvars < 0) call CCTK_WARN (0, "internal error")
  if (any(inputs /= -1 .and. (inputs < 0 .or. inputs >= nvars))) then
     call CCTK_WARN (0, "internal error")
  end if
  
  operand_indices = (/ &
       00, 01, 02, 03, 04, 05, & ! g_ij
       00, 01, 02, 03, 04, 05, & ! g_ij,k
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, & ! g_ij,kl
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       06, 07, 08, 09, 10, 11, & ! K_ij
       06, 07, 08, 09, 10, 11, & ! K_ij,k
       06, 07, 08, 09, 10, 11, &
       06, 07, 08, 09, 10, 11, &
       12, &                     ! alp
       13, 14, 15, &             ! beta^i
       16, &                     ! T_tt
       17, 18, 19, &             ! T_ti
       20, 21, 22, 23, 24, 25, & ! T_ij
       26, 27, 28, &             ! E^i
       29, 30, 31, &             ! A_i
       29, 30, 31, &             ! A_i,J
       29, 30, 31, &
       29, 30, 31, &
       32, 33 /)                 ! phi1, phi2
  
  operation_codes = (/ &
       0, 0, 0, 0, 0, 0, &      ! g_ij
       1, 1, 1, 1, 1, 1, &      ! g_ij,k
       2, 2, 2, 2, 2, 2, &
       3, 3, 3, 3, 3, 3, &
       11, 11, 11, 11, 11, 11, & ! g_ij,kl
       12, 12, 12, 12, 12, 12, &
       13, 13, 13, 13, 13, 13, &
       22, 22, 22, 22, 22, 22, &
       23, 23, 23, 23, 23, 23, &
       33, 33, 33, 33, 33, 33, &
       0, 0, 0, 0, 0, 0, &      ! K_ij
       1, 1, 1, 1, 1, 1, &      ! K_ij,k
       2, 2, 2, 2, 2, 2, &
       3, 3, 3, 3, 3, 3, &
       0, &                     ! alp
       0, 0, 0, &               ! beta^i
       0, &                     ! T_tt
       0, 0, 0, &               ! T_ti
       0, 0, 0, 0, 0, 0, &      ! T_ij
       0, 0, 0, &               ! E^i
       0, 0, 0, &               ! A_i
       1, 1, 1, &               ! A_i,j
       2, 2, 2, &
       3, 3, 3, &
       0, 0 /)                  ! phi1, phi2

  output_types(:) = CCTK_VARIABLE_REAL
  if (hn > 0) then
     outputs = (/ &
          P(qlm_em_gxx), P(qlm_em_gxy), P(qlm_em_gxz), P(qlm_em_gyy), P(qlm_em_gyz), P(qlm_em_gzz), &
          P(qlm_em_dgxxx), P(qlm_em_dgxyx), P(qlm_em_dgxzx), P(qlm_em_dgyyx), P(qlm_em_dgyzx), P(qlm_em_dgzzx), &
          P(qlm_em_dgxxy), P(qlm_em_dgxyy), P(qlm_em_dgxzy), P(qlm_em_dgyyy), P(qlm_em_dgyzy), P(qlm_em_dgzzy), &
          P(qlm_em_dgxxz), P(qlm_em_dgxyz), P(qlm_em_dgxzz), P(qlm_em_dgyyz), P(qlm_em_dgyzz), P(qlm_em_dgzzz), &
          P(qlm_em_ddgxxxx), P(qlm_em_ddgxyxx), P(qlm_em_ddgxzxx), P(qlm_em_ddgyyxx), P(qlm_em_ddgyzxx), P(qlm_em_ddgzzxx), &
          P(qlm_em_ddgxxxy), P(qlm_em_ddgxyxy), P(qlm_em_ddgxzxy), P(qlm_em_ddgyyxy), P(qlm_em_ddgyzxy), P(qlm_em_ddgzzxy), &
          P(qlm_em_ddgxxxz), P(qlm_em_ddgxyxz), P(qlm_em_ddgxzxz), P(qlm_em_ddgyyxz), P(qlm_em_ddgyzxz), P(qlm_em_ddgzzxz), &
          P(qlm_em_ddgxxyy), P(qlm_em_ddgxyyy), P(qlm_em_ddgxzyy), P(qlm_em_ddgyyyy), P(qlm_em_ddgyzyy), P(qlm_em_ddgzzyy), &
          P(qlm_em_ddgxxyz), P(qlm_em_ddgxyyz), P(qlm_em_ddgxzyz), P(qlm_em_ddgyyyz), P(qlm_em_ddgyzyz), P(qlm_em_ddgzzyz), &
          P(qlm_em_ddgxxzz), P(qlm_em_ddgxyzz), P(qlm_em_ddgxzzz), P(qlm_em_ddgyyzz), P(qlm_em_ddgyzzz), P(qlm_em_ddgzzzz), &
          P(qlm_em_kxx), P(qlm_em_kxy), P(qlm_em_kxz), P(qlm_em_kyy), P(qlm_em_kyz), P(qlm_em_kzz), &
          P(qlm_em_dkxxx), P(qlm_em_dkxyx), P(qlm_em_dkxzx), P(qlm_em_dkyyx), P(qlm_em_dkyzx), P(qlm_em_dkzzx), &
          P(qlm_em_dkxxy), P(qlm_em_dkxyy), P(qlm_em_dkxzy), P(qlm_em_dkyyy), P(qlm_em_dkyzy), P(qlm_em_dkzzy), &
          P(qlm_em_dkxxz), P(qlm_em_dkxyz), P(qlm_em_dkxzz), P(qlm_em_dkyyz), P(qlm_em_dkyzz), P(qlm_em_dkzzz), &
          P(qlm_em_alpha), &
          P(qlm_em_betax), P(qlm_em_betay), P(qlm_em_betaz), &
          P(qlm_em_ttt), &
          P(qlm_em_ttx), P(qlm_em_tty), P(qlm_em_ttz), &
          P(qlm_em_txx), P(qlm_em_txy), P(qlm_em_txz), P(qlm_em_tyy), P(qlm_em_tyz), P(qlm_em_tzz), &
          P(qlm_em_ex), P(qlm_em_ey), P(qlm_em_ez), P(qlm_em_ax), P(qlm_em_ay), P(qlm_em_az), &
          P(qlm_em_daxx), P(qlm_em_dayx), P(qlm_em_dazx), P(qlm_em_daxy), P(qlm_em_dayy), P(qlm_em_dazy), &
          P(qlm_em_daxz), P(qlm_em_dayz), P(qlm_em_dazz),&
          P(qlm_em_phi1), P(qlm_em_phi2) /)
  else
     outputs(:) = CCTK_NullPointer()
  end if
  
  
  
  ninputs = size(inputs)
  noutputs = size(outputs)
  
  
  
#if 0
  ! Poison the output variables
  call poison (qlm_em_gxx    )
  call poison (qlm_em_gxy    )
  call poison (qlm_em_gxz    )
  call poison (qlm_em_gyy    )
  call poison (qlm_em_gyz    )
  call poison (qlm_em_gzz    )
  call poison (qlm_em_dgxxx  )
  call poison (qlm_em_dgxyx  )
  call poison (qlm_em_dgxzx  )
  call poison (qlm_em_dgyyx  )
  call poison (qlm_em_dgyzx  )
  call poison (qlm_em_dgzzx  )
  call poison (qlm_em_dgxxy  )
  call poison (qlm_em_dgxyy  )
  call poison (qlm_em_dgxzy  )
  call poison (qlm_em_dgyyy  )
  call poison (qlm_em_dgyzy  )
  call poison (qlm_em_dgzzy  )
  call poison (qlm_em_dgxxz  )
  call poison (qlm_em_dgxyz  )
  call poison (qlm_em_dgxzz  )
  call poison (qlm_em_dgyyz  )
  call poison (qlm_em_dgyzz  )
  call poison (qlm_em_dgzzz  )
  call poison (qlm_em_ddgxxxx)
  call poison (qlm_em_ddgxyxx)
  call poison (qlm_em_ddgxzxx)
  call poison (qlm_em_ddgyyxx)
  call poison (qlm_em_ddgyzxx)
  call poison (qlm_em_ddgzzxx)
  call poison (qlm_em_ddgxxxy)
  call poison (qlm_em_ddgxyxy)
  call poison (qlm_em_ddgxzxy)
  call poison (qlm_em_ddgyyxy)
  call poison (qlm_em_ddgyzxy)
  call poison (qlm_em_ddgzzxy)
  call poison (qlm_em_ddgxxxz)
  call poison (qlm_em_ddgxyxz)
  call poison (qlm_em_ddgxzxz)
  call poison (qlm_em_ddgyyxz)
  call poison (qlm_em_ddgyzxz)
  call poison (qlm_em_ddgzzxz)
  call poison (qlm_em_ddgxxyy)
  call poison (qlm_em_ddgxyyy)
  call poison (qlm_em_ddgxzyy)
  call poison (qlm_em_ddgyyyy)
  call poison (qlm_em_ddgyzyy)
  call poison (qlm_em_ddgzzyy)
  call poison (qlm_em_ddgxxyz)
  call poison (qlm_em_ddgxyyz)
  call poison (qlm_em_ddgxzyz)
  call poison (qlm_em_ddgyyyz)
  call poison (qlm_em_ddgyzyz)
  call poison (qlm_em_ddgzzyz)
  call poison (qlm_em_ddgxxzz)
  call poison (qlm_em_ddgxyzz)
  call poison (qlm_em_ddgxzzz)
  call poison (qlm_em_ddgyyzz)
  call poison (qlm_em_ddgyzzz)
  call poison (qlm_em_ddgzzzz)
  call poison (qlm_em_kxx    )
  call poison (qlm_em_kxy    )
  call poison (qlm_em_kxz    )
  call poison (qlm_em_kyy    )
  call poison (qlm_em_kyz    )
  call poison (qlm_em_kzz    )
  call poison (qlm_em_dkxxx  )
  call poison (qlm_em_dkxyx  )
  call poison (qlm_em_dkxzx  )
  call poison (qlm_em_dkyyx  )
  call poison (qlm_em_dkyzx  )
  call poison (qlm_em_dkzzx  )
  call poison (qlm_em_dkxxy  )
  call poison (qlm_em_dkxyy  )
  call poison (qlm_em_dkxzy  )
  call poison (qlm_em_dkyyy  )
  call poison (qlm_em_dkyzy  )
  call poison (qlm_em_dkzzy  )
  call poison (qlm_em_dkxxz  )
  call poison (qlm_em_dkxyz  )
  call poison (qlm_em_dkxzz  )
  call poison (qlm_em_dkyyz  )
  call poison (qlm_em_dkyzz  )
  call poison (qlm_em_dkzzz  )
  call poison (qlm_em_alpha  )
  call poison (qlm_em_betax  )
  call poison (qlm_em_betay  )
  call poison (qlm_em_betaz  )
  call poison (qlm_em_ttt    )
  call poison (qlm_em_ttx    )
  call poison (qlm_em_tty    )
  call poison (qlm_em_ttz    )
  call poison (qlm_em_txx    )
  call poison (qlm_em_txy    )
  call poison (qlm_em_txz    )
  call poison (qlm_em_tyy    )
  call poison (qlm_em_tyz    )
  call poison (qlm_em_tzz    )
  call poison (qlm_em_ex     )
  call poison (qlm_em_ey     )
  call poison (qlm_em_ez     )
  call poison (qlm_em_ax     )
  call poison (qlm_em_ay     )
  call poison (qlm_em_az     )
  call poison (qlm_em_daxx   )
  call poison (qlm_em_daxy   )
  call poison (qlm_em_daxz   )
  call poison (qlm_em_dayx   )
  call poison (qlm_em_dayy   )
  call poison (qlm_em_dayz   )
  call poison (qlm_em_dazx   )
  call poison (qlm_em_dazy   )
  call poison (qlm_em_dazz   )
  call poison (qlm_em_phi1   )
  call poison (qlm_em_phi2   )
#endif
  

  
  ! Call the interpolator
  call Util_TableSetIntArray &
       (ierr, options_table, noutputs, &
       operand_indices, "operand_indices")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  call Util_TableSetIntArray &
       (ierr, options_table, noutputs, &
       operation_codes, "operation_codes")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  call CCTK_InterpGridArrays &
       (ierr, cctkGH, 3, &
       interp_handle, options_table, coord_handle, &
       npoints, coord_type, coords, &
       ninputs, inputs, &
       noutputs, output_types, outputs)
   call CCTK_INFO ("CCTK_InterpGridArrays")
  
  if (ierr /= 0) then
     if (hn > 0) then
        qlm_em_calc_error(hn) = 1
     end if
     call CCTK_WARN (1, "Interpolator failed")
     return
  end if
  
  
  
  ! Unpack the variables
  if (hn > 0) then
     
     call unpack (qlm_em_gxx    , ni, nj)
     call unpack (qlm_em_gxy    , ni, nj)
     call unpack (qlm_em_gxz    , ni, nj)
     call unpack (qlm_em_gyy    , ni, nj)
     call unpack (qlm_em_gyz    , ni, nj)
     call unpack (qlm_em_gzz    , ni, nj)
     call unpack (qlm_em_dgxxx  , ni, nj)
     call unpack (qlm_em_dgxyx  , ni, nj)
     call unpack (qlm_em_dgxzx  , ni, nj)
     call unpack (qlm_em_dgyyx  , ni, nj)
     call unpack (qlm_em_dgyzx  , ni, nj)
     call unpack (qlm_em_dgzzx  , ni, nj)
     call unpack (qlm_em_dgxxy  , ni, nj)
     call unpack (qlm_em_dgxyy  , ni, nj)
     call unpack (qlm_em_dgxzy  , ni, nj)
     call unpack (qlm_em_dgyyy  , ni, nj)
     call unpack (qlm_em_dgyzy  , ni, nj)
     call unpack (qlm_em_dgzzy  , ni, nj)
     call unpack (qlm_em_dgxxz  , ni, nj)
     call unpack (qlm_em_dgxyz  , ni, nj)
     call unpack (qlm_em_dgxzz  , ni, nj)
     call unpack (qlm_em_dgyyz  , ni, nj)
     call unpack (qlm_em_dgyzz  , ni, nj)
     call unpack (qlm_em_dgzzz  , ni, nj)
     call unpack (qlm_em_ddgxxxx, ni, nj)
     call unpack (qlm_em_ddgxyxx, ni, nj)
     call unpack (qlm_em_ddgxzxx, ni, nj)
     call unpack (qlm_em_ddgyyxx, ni, nj)
     call unpack (qlm_em_ddgyzxx, ni, nj)
     call unpack (qlm_em_ddgzzxx, ni, nj)
     call unpack (qlm_em_ddgxxxy, ni, nj)
     call unpack (qlm_em_ddgxyxy, ni, nj)
     call unpack (qlm_em_ddgxzxy, ni, nj)
     call unpack (qlm_em_ddgyyxy, ni, nj)
     call unpack (qlm_em_ddgyzxy, ni, nj)
     call unpack (qlm_em_ddgzzxy, ni, nj)
     call unpack (qlm_em_ddgxxxz, ni, nj)
     call unpack (qlm_em_ddgxyxz, ni, nj)
     call unpack (qlm_em_ddgxzxz, ni, nj)
     call unpack (qlm_em_ddgyyxz, ni, nj)
     call unpack (qlm_em_ddgyzxz, ni, nj)
     call unpack (qlm_em_ddgzzxz, ni, nj)
     call unpack (qlm_em_ddgxxyy, ni, nj)
     call unpack (qlm_em_ddgxyyy, ni, nj)
     call unpack (qlm_em_ddgxzyy, ni, nj)
     call unpack (qlm_em_ddgyyyy, ni, nj)
     call unpack (qlm_em_ddgyzyy, ni, nj)
     call unpack (qlm_em_ddgzzyy, ni, nj)
     call unpack (qlm_em_ddgxxyz, ni, nj)
     call unpack (qlm_em_ddgxyyz, ni, nj)
     call unpack (qlm_em_ddgxzyz, ni, nj)
     call unpack (qlm_em_ddgyyyz, ni, nj)
     call unpack (qlm_em_ddgyzyz, ni, nj)
     call unpack (qlm_em_ddgzzyz, ni, nj)
     call unpack (qlm_em_ddgxxzz, ni, nj)
     call unpack (qlm_em_ddgxyzz, ni, nj)
     call unpack (qlm_em_ddgxzzz, ni, nj)
     call unpack (qlm_em_ddgyyzz, ni, nj)
     call unpack (qlm_em_ddgyzzz, ni, nj)
     call unpack (qlm_em_ddgzzzz, ni, nj)
     call unpack (qlm_em_kxx    , ni, nj)
     call unpack (qlm_em_kxy    , ni, nj)
     call unpack (qlm_em_kxz    , ni, nj)
     call unpack (qlm_em_kyy    , ni, nj)
     call unpack (qlm_em_kyz    , ni, nj)
     call unpack (qlm_em_kzz    , ni, nj)
     call unpack (qlm_em_dkxxx  , ni, nj)
     call unpack (qlm_em_dkxyx  , ni, nj)
     call unpack (qlm_em_dkxzx  , ni, nj)
     call unpack (qlm_em_dkyyx  , ni, nj)
     call unpack (qlm_em_dkyzx  , ni, nj)
     call unpack (qlm_em_dkzzx  , ni, nj)
     call unpack (qlm_em_dkxxy  , ni, nj)
     call unpack (qlm_em_dkxyy  , ni, nj)
     call unpack (qlm_em_dkxzy  , ni, nj)
     call unpack (qlm_em_dkyyy  , ni, nj)
     call unpack (qlm_em_dkyzy  , ni, nj)
     call unpack (qlm_em_dkzzy  , ni, nj)
     call unpack (qlm_em_dkxxz  , ni, nj)
     call unpack (qlm_em_dkxyz  , ni, nj)
     call unpack (qlm_em_dkxzz  , ni, nj)
     call unpack (qlm_em_dkyyz  , ni, nj)
     call unpack (qlm_em_dkyzz  , ni, nj)
     call unpack (qlm_em_dkzzz  , ni, nj)
     call unpack (qlm_em_alpha  , ni, nj)
     call unpack (qlm_em_betax  , ni, nj)
     call unpack (qlm_em_betay  , ni, nj)
     call unpack (qlm_em_betaz  , ni, nj)
     if (stress_energy_state /= 0) then
        call unpack (qlm_em_ttt    , ni, nj)
        call unpack (qlm_em_ttx    , ni, nj)
        call unpack (qlm_em_tty    , ni, nj)
        call unpack (qlm_em_ttz    , ni, nj)
        call unpack (qlm_em_txx    , ni, nj)
        call unpack (qlm_em_txy    , ni, nj)
        call unpack (qlm_em_txz    , ni, nj)
        call unpack (qlm_em_tyy    , ni, nj)
        call unpack (qlm_em_tyz    , ni, nj)
        call unpack (qlm_em_tzz    , ni, nj)
     else
        qlm_em_ttt = 0
        qlm_em_ttx = 0
        qlm_em_tty = 0
        qlm_em_ttz = 0
        qlm_em_txx = 0
        qlm_em_txy = 0
        qlm_em_txz = 0
        qlm_em_tyy = 0
        qlm_em_tyz = 0
        qlm_em_tzz = 0
     end if
     if (calc_charge /= 0) then
        call unpack (qlm_em_ex    , ni, nj)
        call unpack (qlm_em_ey    , ni, nj)
        call unpack (qlm_em_ez    , ni, nj)
        call unpack (qlm_em_ax    , ni, nj)
        call unpack (qlm_em_ay    , ni, nj)
        call unpack (qlm_em_az    , ni, nj)
        call unpack (qlm_em_daxx  , ni, nj)
        call unpack (qlm_em_daxy  , ni, nj)
        call unpack (qlm_em_daxz  , ni, nj)
        call unpack (qlm_em_dayx  , ni, nj)
        call unpack (qlm_em_dayy  , ni, nj)
        call unpack (qlm_em_dayz  , ni, nj)
        call unpack (qlm_em_dazx  , ni, nj)
        call unpack (qlm_em_dazy  , ni, nj)
        call unpack (qlm_em_dazz  , ni, nj)
     else
        qlm_em_ex = 0
        qlm_em_ey = 0
        qlm_em_ez = 0
        qlm_em_ax = 0
        qlm_em_ay = 0
        qlm_em_az = 0
        qlm_em_daxx = 0
        qlm_em_daxy = 0
        qlm_em_daxz = 0
        qlm_em_dayx = 0
        qlm_em_dayy = 0
        qlm_em_dayz = 0
        qlm_em_dazx = 0
        qlm_em_dazy = 0
        qlm_em_dazz = 0

     end if
     call unpack (qlm_em_phi1     , ni, nj)
     call unpack (qlm_em_phi2     , ni, nj)

     
     
#if 0
     ! Check for poison
     call poison_check (qlm_em_gxx    , "qlm_em_gxx    ")
     call poison_check (qlm_em_gxy    , "qlm_em_gxy    ")
     call poison_check (qlm_em_gxz    , "qlm_em_gxz    ")
     call poison_check (qlm_em_gyy    , "qlm_em_gyy    ")
     call poison_check (qlm_em_gyz    , "qlm_em_gyz    ")
     call poison_check (qlm_em_gzz    , "qlm_em_gzz    ")
     call poison_check (qlm_em_dgxxx  , "qlm_em_dgxxx  ")
     call poison_check (qlm_em_dgxyx  , "qlm_em_dgxyx  ")
     call poison_check (qlm_em_dgxzx  , "qlm_em_dgxzx  ")
     call poison_check (qlm_em_dgyyx  , "qlm_em_dgyyx  ")
     call poison_check (qlm_em_dgyzx  , "qlm_em_dgyzx  ")
     call poison_check (qlm_em_dgzzx  , "qlm_em_dgzzx  ")
     call poison_check (qlm_em_dgxxy  , "qlm_em_dgxxy  ")
     call poison_check (qlm_em_dgxyy  , "qlm_em_dgxyy  ")
     call poison_check (qlm_em_dgxzy  , "qlm_em_dgxzy  ")
     call poison_check (qlm_em_dgyyy  , "qlm_em_dgyyy  ")
     call poison_check (qlm_em_dgyzy  , "qlm_em_dgyzy  ")
     call poison_check (qlm_em_dgzzy  , "qlm_em_dgzzy  ")
     call poison_check (qlm_em_dgxxz  , "qlm_em_dgxxz  ")
     call poison_check (qlm_em_dgxyz  , "qlm_em_dgxyz  ")
     call poison_check (qlm_em_dgxzz  , "qlm_em_dgxzz  ")
     call poison_check (qlm_em_dgyyz  , "qlm_em_dgyyz  ")
     call poison_check (qlm_em_dgyzz  , "qlm_em_dgyzz  ")
     call poison_check (qlm_em_dgzzz  , "qlm_em_dgzzz  ")
     call poison_check (qlm_em_ddgxxxx, "qlm_em_ddgxxxx")
     call poison_check (qlm_em_ddgxyxx, "qlm_em_ddgxyxx")
     call poison_check (qlm_em_ddgxzxx, "qlm_em_ddgxzxx")
     call poison_check (qlm_em_ddgyyxx, "qlm_em_ddgyyxx")
     call poison_check (qlm_em_ddgyzxx, "qlm_em_ddgyzxx")
     call poison_check (qlm_em_ddgzzxx, "qlm_em_ddgzzxx")
     call poison_check (qlm_em_ddgxxxy, "qlm_em_ddgxxxy")
     call poison_check (qlm_em_ddgxyxy, "qlm_em_ddgxyxy")
     call poison_check (qlm_em_ddgxzxy, "qlm_em_ddgxzxy")
     call poison_check (qlm_em_ddgyyxy, "qlm_em_ddgyyxy")
     call poison_check (qlm_em_ddgyzxy, "qlm_em_ddgyzxy")
     call poison_check (qlm_em_ddgzzxy, "qlm_em_ddgzzxy")
     call poison_check (qlm_em_ddgxxxz, "qlm_em_ddgxxxz")
     call poison_check (qlm_em_ddgxyxz, "qlm_em_ddgxyxz")
     call poison_check (qlm_em_ddgxzxz, "qlm_em_ddgxzxz")
     call poison_check (qlm_em_ddgyyxz, "qlm_em_ddgyyxz")
     call poison_check (qlm_em_ddgyzxz, "qlm_em_ddgyzxz")
     call poison_check (qlm_em_ddgzzxz, "qlm_em_ddgzzxz")
     call poison_check (qlm_em_ddgxxyy, "qlm_em_ddgxxyy")
     call poison_check (qlm_em_ddgxyyy, "qlm_em_ddgxyyy")
     call poison_check (qlm_em_ddgxzyy, "qlm_em_ddgxzyy")
     call poison_check (qlm_em_ddgyyyy, "qlm_em_ddgyyyy")
     call poison_check (qlm_em_ddgyzyy, "qlm_em_ddgyzyy")
     call poison_check (qlm_em_ddgzzyy, "qlm_em_ddgzzyy")
     call poison_check (qlm_em_ddgxxyz, "qlm_em_ddgxxyz")
     call poison_check (qlm_em_ddgxyyz, "qlm_em_ddgxyyz")
     call poison_check (qlm_em_ddgxzyz, "qlm_em_ddgxzyz")
     call poison_check (qlm_em_ddgyyyz, "qlm_em_ddgyyyz")
     call poison_check (qlm_em_ddgyzyz, "qlm_em_ddgyzyz")
     call poison_check (qlm_em_ddgzzyz, "qlm_em_ddgzzyz")
     call poison_check (qlm_em_ddgxxzz, "qlm_em_ddgxxzz")
     call poison_check (qlm_em_ddgxyzz, "qlm_em_ddgxyzz")
     call poison_check (qlm_em_ddgxzzz, "qlm_em_ddgxzzz")
     call poison_check (qlm_em_ddgyyzz, "qlm_em_ddgyyzz")
     call poison_check (qlm_em_ddgyzzz, "qlm_em_ddgyzzz")
     call poison_check (qlm_em_ddgzzzz, "qlm_em_ddgzzzz")
     call poison_check (qlm_em_kxx    , "qlm_em_kxx    ")
     call poison_check (qlm_em_kxy    , "qlm_em_kxy    ")
     call poison_check (qlm_em_kxz    , "qlm_em_kxz    ")
     call poison_check (qlm_em_kyy    , "qlm_em_kyy    ")
     call poison_check (qlm_em_kyz    , "qlm_em_kyz    ")
     call poison_check (qlm_em_kzz    , "qlm_em_kzz    ")
     call poison_check (qlm_em_dkxxx  , "qlm_em_dkxxx  ")
     call poison_check (qlm_em_dkxyx  , "qlm_em_dkxyx  ")
     call poison_check (qlm_em_dkxzx  , "qlm_em_dkxzx  ")
     call poison_check (qlm_em_dkyyx  , "qlm_em_dkyyx  ")
     call poison_check (qlm_em_dkyzx  , "qlm_em_dkyzx  ")
     call poison_check (qlm_em_dkzzx  , "qlm_em_dkzzx  ")
     call poison_check (qlm_em_dkxxy  , "qlm_em_dkxxy  ")
     call poison_check (qlm_em_dkxyy  , "qlm_em_dkxyy  ")
     call poison_check (qlm_em_dkxzy  , "qlm_em_dkxzy  ")
     call poison_check (qlm_em_dkyyy  , "qlm_em_dkyyy  ")
     call poison_check (qlm_em_dkyzy  , "qlm_em_dkyzy  ")
     call poison_check (qlm_em_dkzzy  , "qlm_em_dkzzy  ")
     call poison_check (qlm_em_dkxxz  , "qlm_em_dkxxz  ")
     call poison_check (qlm_em_dkxyz  , "qlm_em_dkxyz  ")
     call poison_check (qlm_em_dkxzz  , "qlm_em_dkxzz  ")
     call poison_check (qlm_em_dkyyz  , "qlm_em_dkyyz  ")
     call poison_check (qlm_em_dkyzz  , "qlm_em_dkyzz  ")
     call poison_check (qlm_em_dkzzz  , "qlm_em_dkzzz  ")
     call poison_check (qlm_em_alpha  , "qlm_em_alpha  ")
     call poison_check (qlm_em_betax  , "qlm_em_betax  ")
     call poison_check (qlm_em_betay  , "qlm_em_betay  ")
     call poison_check (qlm_em_betaz  , "qlm_em_betaz  ")
     call poison_check (qlm_em_ttt    , "qlm_em_ttt    ")
     call poison_check (qlm_em_ttx    , "qlm_em_ttx    ")
     call poison_check (qlm_em_tty    , "qlm_em_tty    ")
     call poison_check (qlm_em_ttz    , "qlm_em_ttz    ")
     call poison_check (qlm_em_txx    , "qlm_em_txx    ")
     call poison_check (qlm_em_txy    , "qlm_em_txy    ")
     call poison_check (qlm_em_txz    , "qlm_em_txz    ")
     call poison_check (qlm_em_tyy    , "qlm_em_tyy    ")
     call poison_check (qlm_em_tyz    , "qlm_em_tyz    ")
     call poison_check (qlm_em_tzz    , "qlm_em_tzz    ")
     call poison_check (qlm_em_ex     , "qlm_em_ex     ")
     call poison_check (qlm_em_ey     , "qlm_em_ey     ")
     call poison_check (qlm_em_ez     , "qlm_em_ez     ")
     call poison_check (qlm_em_ax     , "qlm_em_ax     ")
     call poison_check (qlm_em_ay     , "qlm_em_ay     ")
     call poison_check (qlm_em_az     , "qlm_em_az     ")
     call poison_check (qlm_em_daxx   , "qlm_em_daxx   ")
     call poison_check (qlm_em_daxy   , "qlm_em_daxy   ")
     call poison_check (qlm_em_daxz   , "qlm_em_daxz   ")
     call poison_check (qlm_em_dayx   , "qlm_em_dayx   ")
     call poison_check (qlm_em_dayy   , "qlm_em_dayy   ")
     call poison_check (qlm_em_dayz   , "qlm_em_dayz   ")
     call poison_check (qlm_em_dazx   , "qlm_em_dazx   ")
     call poison_check (qlm_em_dazy   , "qlm_em_dazy   ")
     call poison_check (qlm_em_dazz   , "qlm_em_dazz   ")
     call poison_check (qlm_em_phi1   , "qlm_em_phi1   ")
     call poison_check (qlm_em_phi2   , "qlm_em_phi2   ")
#endif
     
  end if
  

  
  ! Free interpolator options
  call Util_TableDestroy (ierr, options_table)
  
  
  
  if (hn > 0) then
     
     qlm_em_have_valid_data(hn) = 1
     
     deallocate (xcoord)
     deallocate (ycoord)
     deallocate (zcoord)
     
  end if

  
  
contains
  
  subroutine pack (arr, ni, nj)
    integer,   intent(in)    :: ni, nj
    CCTK_REAL, intent(inout) :: arr(:,:)
    CCTK_REAL :: tmp(ni,nj)
    tmp(:,:) = arr(:ni, :nj)
    call copy (arr, tmp, size(tmp))
  end subroutine pack
  
  subroutine unpack (arr, ni, nj)
    integer,   intent(in)    :: ni, nj
    CCTK_REAL, intent(inout) :: arr(:,:)
    CCTK_REAL :: tmp(ni,nj)
    call copy (tmp, arr, size(tmp))
    arr(:ni, :nj) = tmp(:,:)
    arr(ni+1:, :nj) = 0
    arr(:, nj+1:) = 0
  end subroutine unpack
  
  subroutine copy (a, b, n)
    integer,   intent(in)  :: n
    CCTK_REAL, intent(out) :: a(n)
    CCTK_REAL, intent(in)  :: b(n)
    a = b
  end subroutine copy
  
  subroutine poison (arr)
    CCTK_REAL, intent(out) :: arr(:,:)
    arr = poison_value
  end subroutine poison
  
  subroutine poison_check (arr, name)
    CCTK_REAL,    intent(in) :: arr(:,:)
    character(*), intent(in) :: name
    character*1000 :: msg
!!$    integer        :: i, j
    if (any(arr==poison_value)) then
       write (msg, '("Poison found in ",a)') trim(name)
       call CCTK_WARN (CCTK_WARN_ALERT, msg)
!!$       do j=1,size(arr,2)
!!$          do i=1,size(arr,1)
!!$             print '(2i6)', i,j
!!$          end do
!!$       end do
    end if
  end subroutine poison_check
  
end subroutine qlm_em_interpolate
