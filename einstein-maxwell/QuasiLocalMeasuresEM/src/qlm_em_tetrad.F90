#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_calc_tetrad (CCTK_ARGUMENTS, hn)
  use qlm_em_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_em_calc_tetrad1 (CCTK_PASS_FTOF, hn)
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_l0(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_l1(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_l2(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_l3(:,:,hn), +1)
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_n0(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_n1(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_n2(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_n3(:,:,hn), +1)
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_m0(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_m1(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_m2(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_m3(:,:,hn), +1)
  
end subroutine qlm_em_calc_tetrad
