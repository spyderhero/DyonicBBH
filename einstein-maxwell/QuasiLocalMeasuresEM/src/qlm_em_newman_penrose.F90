#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_calc_newman_penrose (CCTK_ARGUMENTS, hn)
  use qlm_em_boundary
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  call qlm_em_calc_newman_penrose1 (CCTK_PASS_FTOF, hn)

  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_npkappa  (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_nptau    (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_npsigma  (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_nprho    (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_npepsilon(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_npgamma  (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_npbeta   (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_npalpha  (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_nppi     (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_npnu     (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_npmu     (:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_nplambda (:,:,hn), +1)
  
end subroutine qlm_em_calc_newman_penrose
