#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "EinsteinMaxwell_FD.h"

! Current ADM geometry is refreshed by the schedule AFTER input synchronization.
! Diagnose synchronized final-state analysis times, never RK intermediate stages.
subroutine EinsteinMaxwell_calc_constraints(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  CCTK_INT i,j,k,a,b,c,m,istat,nstencil
  CCTK_REAL hh(3,3),hu(3,3),dethh,ch,lE(3),lB(3)
  CCTK_REAL d1_ch(3),d1_hh(3,3,3),d1_lE(3,3),d1_lB(3,3)
  CCTK_REAL cf1(3,3,3),cf2(3,3,3),divE,divB,jac(3,3),radius
  CCTK_REAL dx12,dy12,dz12,odx60,ody60,odz60
  logical use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11,lJ12,lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21,lJ22,lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31,lJ32,lJ33
  CCTK_POINTER lJ11_ptr,lJ12_ptr,lJ13_ptr,lJ21_ptr,lJ22_ptr,lJ23_ptr
  CCTK_POINTER lJ31_ptr,lJ32_ptr,lJ33_ptr
  pointer (lJ11_ptr,lJ11),(lJ12_ptr,lJ12),(lJ13_ptr,lJ13)
  pointer (lJ21_ptr,lJ21),(lJ22_ptr,lJ22),(lJ23_ptr,lJ23)
  pointer (lJ31_ptr,lJ31),(lJ32_ptr,lJ32),(lJ33_ptr,lJ33)

  if (compute_em_constraints == 0) return
  if (mod(cctk_iteration,em_constraints_every) /= 0) return
  nstencil=derivs_order/2
  if (any(cctk_nghostzones < nstencil)) &
    call CCTK_WARN(0,"Insufficient ghost zones for Maxwell diagnostics")
  ! ADM metric ghost data was synchronized by ConstraintsSync at this time.
  call EinsteinMaxwell_adm2bssn(CCTK_PASS_FTOF)
  call CCTK_IsFunctionAliased(istat,"MultiPatch_GetDomainSpecification")
  use_jacobian=istat /= 0
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ11_ptr,cctkGH,0,"Coordinates::J11")
    if (lJ11_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ12_ptr,cctkGH,0,"Coordinates::J12")
    if (lJ12_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ13_ptr,cctkGH,0,"Coordinates::J13")
    if (lJ13_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ21_ptr,cctkGH,0,"Coordinates::J21")
    if (lJ21_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ22_ptr,cctkGH,0,"Coordinates::J22")
    if (lJ22_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ23_ptr,cctkGH,0,"Coordinates::J23")
    if (lJ23_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ31_ptr,cctkGH,0,"Coordinates::J31")
    if (lJ31_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ32_ptr,cctkGH,0,"Coordinates::J32")
    if (lJ32_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if
  if (use_jacobian) then
    call CCTK_VarDataPtr(lJ33_ptr,cctkGH,0,"Coordinates::J33")
    if (lJ33_ptr == 0) call CCTK_WARN(0,"Missing Jacobian for Maxwell diagnostic")
  end if

  dx12=12.0d0*CCTK_DELTA_SPACE(1)
  dy12=12.0d0*CCTK_DELTA_SPACE(2)
  dz12=12.0d0*CCTK_DELTA_SPACE(3)
  odx60=1.0d0/(60.0d0*CCTK_DELTA_SPACE(1))
  ody60=1.0d0/(60.0d0*CCTK_DELTA_SPACE(2))
  odz60=1.0d0/(60.0d0*CCTK_DELTA_SPACE(3))

  ! Uncomputed boundary/ghost cells are zero placeholders, NOT valid residuals.
  EM_CE=0.0d0
  EM_CB=0.0d0
  if (compute_em_constraint_monitors /= 0) then
    EM_CE_near=0.0d0
    EM_CB_near=0.0d0
    EM_mask_near=0.0d0
    EM_CE_far=0.0d0
    EM_CB_far=0.0d0
    EM_mask_far=0.0d0
  end if

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(i,j,k,a,b,c,m,hh,hu,dethh,ch,lE,lB,d1_ch,d1_hh, &
  !$OMP d1_lE,d1_lB,cf1,cf2,divE,divB,jac,radius)
  do k=1+cctk_nghostzones(3),cctk_lsh(3)-cctk_nghostzones(3)
  do j=1+cctk_nghostzones(2),cctk_lsh(2)-cctk_nghostzones(2)
  do i=1+cctk_nghostzones(1),cctk_lsh(1)-cctk_nghostzones(1)
    ch=chi(i,j,k)
    hh(1,1)=hxx(i,j,k)
    hh(1,2)=hxy(i,j,k)
    hh(1,3)=hxz(i,j,k)
    hh(2,2)=hyy(i,j,k)
    hh(2,3)=hyz(i,j,k)
    hh(3,3)=hzz(i,j,k)
    hh(2,1)=hh(1,2)
    hh(3,1)=hh(1,3)
    hh(3,2)=hh(2,3)
    lE=(/Ex(i,j,k),Ey(i,j,k),Ez(i,j,k)/)
    lB=(/Bx(i,j,k),By(i,j,k),Bz(i,j,k)/)
    dethh =       hh(1,1) * hh(2,2) * hh(3,3)                              &
            + 2.0d0 * hh(1,2) * hh(1,3) * hh(2,3)                          &
            -        hh(1,1) * hh(2,3) ** 2                                 &
            -        hh(2,2) * hh(1,3) ** 2                                 &
            -        hh(3,3) * hh(1,2) ** 2
    hu(1,1) = (hh(2,2) * hh(3,3) - hh(2,3) ** 2     ) / dethh
    hu(2,2) = (hh(1,1) * hh(3,3) - hh(1,3) ** 2     ) / dethh
    hu(3,3) = (hh(1,1) * hh(2,2) - hh(1,2) ** 2     ) / dethh
    hu(1,2) = (hh(1,3) * hh(2,3) - hh(1,2) * hh(3,3)) / dethh
    hu(1,3) = (hh(1,2) * hh(2,3) - hh(1,3) * hh(2,2)) / dethh
    hu(2,3) = (hh(1,3) * hh(1,2) - hh(2,3) * hh(1,1)) / dethh
    hu(2,1) = hu(1,2)
    hu(3,1) = hu(1,3)
    hu(3,2) = hu(2,3)


    if (derivs_order == 4) then
      d1_ch(1) = D1X4(chi)
      d1_ch(2) = D1Y4(chi)
      d1_ch(3) = D1Z4(chi)
      d1_hh(1,1,1) = D1X4(hxx)
      d1_hh(1,1,2) = D1Y4(hxx)
      d1_hh(1,1,3) = D1Z4(hxx)
      d1_hh(1,2,1) = D1X4(hxy)
      d1_hh(1,2,2) = D1Y4(hxy)
      d1_hh(1,2,3) = D1Z4(hxy)
      d1_hh(1,3,1) = D1X4(hxz)
      d1_hh(1,3,2) = D1Y4(hxz)
      d1_hh(1,3,3) = D1Z4(hxz)
      d1_hh(2,2,1) = D1X4(hyy)
      d1_hh(2,2,2) = D1Y4(hyy)
      d1_hh(2,2,3) = D1Z4(hyy)
      d1_hh(2,3,1) = D1X4(hyz)
      d1_hh(2,3,2) = D1Y4(hyz)
      d1_hh(2,3,3) = D1Z4(hyz)
      d1_hh(3,3,1) = D1X4(hzz)
      d1_hh(3,3,2) = D1Y4(hzz)
      d1_hh(3,3,3) = D1Z4(hzz)
      d1_hh(2,1,:) = d1_hh(1,2,:)
      d1_hh(3,1,:) = d1_hh(1,3,:)
      d1_hh(3,2,:) = d1_hh(2,3,:)
      d1_lE(1,1) = D1X4(Ex)
      d1_lE(1,2) = D1Y4(Ex)
      d1_lE(1,3) = D1Z4(Ex)
      d1_lE(2,1) = D1X4(Ey)
      d1_lE(2,2) = D1Y4(Ey)
      d1_lE(2,3) = D1Z4(Ey)
      d1_lE(3,1) = D1X4(Ez)
      d1_lE(3,2) = D1Y4(Ez)
      d1_lE(3,3) = D1Z4(Ez)
      d1_lB(1,1) = D1X4(Bx)
      d1_lB(1,2) = D1Y4(Bx)
      d1_lB(1,3) = D1Z4(Bx)
      d1_lB(2,1) = D1X4(By)
      d1_lB(2,2) = D1Y4(By)
      d1_lB(2,3) = D1Z4(By)
      d1_lB(3,1) = D1X4(Bz)
      d1_lB(3,2) = D1Y4(Bz)
      d1_lB(3,3) = D1Z4(Bz)

    else if (derivs_order == 6) then
      d1_ch(1) = D1X6(chi)
      d1_ch(2) = D1Y6(chi)
      d1_ch(3) = D1Z6(chi)
      d1_hh(1,1,1) = D1X6(hxx)
      d1_hh(1,1,2) = D1Y6(hxx)
      d1_hh(1,1,3) = D1Z6(hxx)
      d1_hh(1,2,1) = D1X6(hxy)
      d1_hh(1,2,2) = D1Y6(hxy)
      d1_hh(1,2,3) = D1Z6(hxy)
      d1_hh(1,3,1) = D1X6(hxz)
      d1_hh(1,3,2) = D1Y6(hxz)
      d1_hh(1,3,3) = D1Z6(hxz)
      d1_hh(2,2,1) = D1X6(hyy)
      d1_hh(2,2,2) = D1Y6(hyy)
      d1_hh(2,2,3) = D1Z6(hyy)
      d1_hh(2,3,1) = D1X6(hyz)
      d1_hh(2,3,2) = D1Y6(hyz)
      d1_hh(2,3,3) = D1Z6(hyz)
      d1_hh(3,3,1) = D1X6(hzz)
      d1_hh(3,3,2) = D1Y6(hzz)
      d1_hh(3,3,3) = D1Z6(hzz)
      d1_hh(2,1,:) = d1_hh(1,2,:)
      d1_hh(3,1,:) = d1_hh(1,3,:)
      d1_hh(3,2,:) = d1_hh(2,3,:)
      d1_lE(1,1) = D1X6(Ex)
      d1_lE(1,2) = D1Y6(Ex)
      d1_lE(1,3) = D1Z6(Ex)
      d1_lE(2,1) = D1X6(Ey)
      d1_lE(2,2) = D1Y6(Ey)
      d1_lE(2,3) = D1Z6(Ey)
      d1_lE(3,1) = D1X6(Ez)
      d1_lE(3,2) = D1Y6(Ez)
      d1_lE(3,3) = D1Z6(Ez)
      d1_lB(1,1) = D1X6(Bx)
      d1_lB(1,2) = D1Y6(Bx)
      d1_lB(1,3) = D1Z6(Bx)
      d1_lB(2,1) = D1X6(By)
      d1_lB(2,2) = D1Y6(By)
      d1_lB(2,3) = D1Z6(By)
      d1_lB(3,1) = D1X6(Bz)
      d1_lB(3,2) = D1Y6(Bz)
      d1_lB(3,3) = D1Z6(Bz)

    else
      call CCTK_WARN(0,"Unsupported Maxwell diagnostic derivative order")
    end if
    if (use_jacobian) then
      jac(1,1)=lJ11(i,j,k)
      jac(1,2)=lJ12(i,j,k)
      jac(1,3)=lJ13(i,j,k)
      jac(2,1)=lJ21(i,j,k)
      jac(2,2)=lJ22(i,j,k)
      jac(2,3)=lJ23(i,j,k)
      jac(3,1)=lJ31(i,j,k)
      jac(3,2)=lJ32(i,j,k)
      jac(3,3)=lJ33(i,j,k)
      call EinsteinMaxwell_d1_Scalar_apply_jacobian(d1_ch,jac)
      call EinsteinMaxwell_d1_Vector_apply_jacobian(d1_lE,jac)
      call EinsteinMaxwell_d1_Vector_apply_jacobian(d1_lB,jac)
      call EinsteinMaxwell_d1_2nd_Tensor_apply_jacobian(d1_hh,jac,1)
    end if
#include "EinsteinMaxwell_Gauss.inc"
    EM_CE(i,j,k)=divE
    EM_CB(i,j,k)=divB

    ! Fixed origin-centred shells, NOT a moving apparent-horizon mask.
    if (compute_em_constraint_monitors /= 0) then
      radius=sqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)
      if (radius >= em_monitor_near_rmin .and. radius < em_monitor_near_rmax) then
        EM_mask_near(i,j,k)=1.0d0
        EM_CE_near(i,j,k)=divE
        EM_CB_near(i,j,k)=divB
      end if
      if (radius >= em_monitor_far_rmin .and. radius < em_monitor_far_rmax) then
        EM_mask_far(i,j,k)=1.0d0
        EM_CE_far(i,j,k)=divE
        EM_CB_far(i,j,k)=divB
      end if
    end if
  end do
  end do
  end do
  !$OMP END PARALLEL DO
end subroutine EinsteinMaxwell_calc_constraints

! Runs in level mode to request synchronization of the current input groups.
subroutine EinsteinMaxwell_constraints_sync(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  CCTK_INT ierr
  if (compute_em_constraints == 0) return
  if (mod(cctk_iteration,em_constraints_every) /= 0) return
  call CCTK_SyncGroup(ierr,cctkGH,"ADMBase::metric")
  if (ierr < 0) call CCTK_WARN(0,"Failed to synchronize diagnostic metric")
  call CCTK_SyncGroup(ierr,cctkGH,"ProcaBase::Ei")
  if (ierr < 0) call CCTK_WARN(0,"Failed to synchronize diagnostic electric field")
  call CCTK_SyncGroup(ierr,cctkGH,"EinsteinMaxwellEvolve::Bi")
  if (ierr < 0) call CCTK_WARN(0,"Failed to synchronize diagnostic magnetic field")
end subroutine EinsteinMaxwell_constraints_sync

subroutine EinsteinMaxwell_constraints_check(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  if (compute_em_constraint_monitors /= 0 .and. compute_em_constraints == 0) &
    call CCTK_WARN(0,"Masked Maxwell monitors require compute_em_constraints=yes")
  if (compute_em_constraints /= 0) then
    if (.not. CCTK_EQUALS(evolution_method,"EinsteinMaxwellEvolve")) &
      call CCTK_WARN(0,"Maxwell diagnostics require EinsteinMaxwell evolution")
  end if
  if (compute_em_constraint_monitors /= 0) then
    if (em_monitor_near_rmax <= em_monitor_near_rmin) &
      call CCTK_WARN(0,"Invalid near Maxwell monitoring shell")
    if (em_monitor_far_rmax <= em_monitor_far_rmin) &
      call CCTK_WARN(0,"Invalid far Maxwell monitoring shell")
  end if
end subroutine EinsteinMaxwell_constraints_check
