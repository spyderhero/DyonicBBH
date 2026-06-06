#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_calc_weyl_scalars (CCTK_ARGUMENTS, hn)
  use qlm_em_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_em_calc_weyl_scalars1 (CCTK_PASS_FTOF, hn)

  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_psi0(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_psi1(:,:,hn), -1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_psi2(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_psi3(:,:,hn), -1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_psi4(:,:,hn), +1)
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_i(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_j(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_s(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_sdiff(:,:,hn), +1)
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi00(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi11(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi01(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi12(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi10(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi21(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi02(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi22(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_phi20(:,:,hn), +1)
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_lambda(:,:,hn), +1)
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_lie_n_theta_l(:,:,hn), +1)
 
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_rsc(:,:,hn), +1)
  
end subroutine qlm_em_calc_weyl_scalars
