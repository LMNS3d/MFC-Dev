!>
!! @file m_rhs.f90
!! @brief Contains module m_rhs

#:include 'macros.fpp'

!> @brief The module contains the subroutines used to calculate the right-
!!              hand-side (RHS) in the quasi-conservative, shock- and interface-
!!              capturing finite-volume framework for the multicomponent Navier-
!!              Stokes equations supplemented by appropriate advection equations
!!              used to capture the material interfaces. The system of equations
!!              is closed by the stiffened gas equation of state, as well as any
!!              required mixture relationships. Capillarity effects are included
!!              and are modeled by the means of a volume force acting across the
!!              diffuse material interface region. The implementation details of
!!              surface tension may be found in Perigaud and Saurel (2005). Note
!!              that both viscous and surface tension effects are only available
!!              in the volume fraction model.
module m_rhs

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables
    use m_riemann_solvers      !< Exact and approximate Riemann problem solvers

    use m_cbc                  !< Characteristic boundary conditions (CBC)

    use m_bubbles              !< Bubble dynamic routines

    use m_qbmm                 !< Moment inversion

    use m_hypoelastic

    use m_monopole

    use m_viscous

    use m_ibm

    use m_nvtx

    use m_boundary_conditions

    use m_helper
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_rhs_module, &
 s_compute_rhs, &
 s_pressure_relaxation_procedure, &
 s_pressure_relaxation_Wakimura2024, &
 s_pressure_relaxation_Pelanti2014, &
 s_pressure_relaxation_Pelanti2022, &
 s_finalize_rhs_module

    !! This variable contains the WENO-reconstructed values of the cell-average
    !! conservative variables, which are located in q_cons_vf, at cell-interior
    !! Gaussian quadrature points (QP).
    type(vector_field) :: q_cons_qp !<
    !$acc declare create(q_cons_qp)

    !! The primitive variables at cell-interior Gaussian quadrature points. These
    !! are calculated from the conservative variables and gradient magnitude (GM)
    !! of the volume fractions, q_cons_qp and gm_alpha_qp, respectively.
    type(vector_field) :: q_prim_qp !<
    !$acc declare create(q_prim_qp)

    !> @name The first-order spatial derivatives of the primitive variables at cell-
    !! interior Gaussian quadrature points. These are WENO-reconstructed from
    !! their respective cell-average values, obtained through the application
    !! of the divergence theorem on the integral-average cell-boundary values
    !! of the primitive variables, located in qK_prim_n, where K = L or R.
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp)
    !$acc declare link(dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp)
#else
    type(vector_field), allocatable, dimension(:) :: dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp
    !$acc declare create(dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp)
#endif

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell-
    !! average first-order spatial derivatives of the primitive variables. The
    !! cell-average of the first-order spatial derivatives may be found in the
    !! variables dq_prim_ds_qp, where s = x, y or z.
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n)
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n)
    !$acc declare link(dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n)
    !$acc declare link(dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n)
#else
    type(vector_field), allocatable, dimension(:) :: dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n
    type(vector_field), allocatable, dimension(:) :: dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n
    !$acc declare create(dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n)
    !$acc declare create(dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n)
#endif
    !> @}

    type(vector_field) :: gm_alpha_qp  !<
    !! The gradient magnitude of the volume fractions at cell-interior Gaussian
    !! quadrature points. gm_alpha_qp is calculated from individual first-order
    !! spatial derivatives located in dq_prim_ds_qp.

    !$acc declare create(gm_alpha_qp)

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell-
    !! average gradient magnitude of volume fractions, located in gm_alpha_qp.
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), gm_alphaL_n)
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), gm_alphaR_n)
    !$acc declare link(gm_alphaL_n, gm_alphaR_n)
#else
    type(vector_field), allocatable, dimension(:) :: gm_alphaL_n
    type(vector_field), allocatable, dimension(:) :: gm_alphaR_n
    !$acc declare create(gm_alphaL_n, gm_alphaR_n)
#endif
    !> @}

    !> @name The cell-boundary values of the fluxes (src - source, gsrc - geometrical
    !! source). These are computed by applying the chosen Riemann problem solver
    !! .on the left and right cell-boundary values of the primitive variables
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), flux_n, lo_flux_n)
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), flux_src_n, lo_flux_src_n)
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), flux_gsrc_n, lo_flux_gsrc_n)
    !$acc declare link(flux_n, flux_src_n, flux_gsrc_n, lo_flux_n, lo_flux_src_n, lo_flux_gsrc_n)
#else
    type(vector_field), allocatable, dimension(:) :: flux_n, lo_flux_n
    type(vector_field), allocatable, dimension(:) :: flux_src_n, lo_flux_src_n
    type(vector_field), allocatable, dimension(:) :: flux_gsrc_n, lo_flux_gsrc_n
    !$acc declare create(flux_n, flux_src_n, flux_gsrc_n, lo_flux_n, lo_flux_src_n, lo_flux_gsrc_n)
#endif
    !> @}

    !> @name The interfacial fluctuation terms in the wave-propagation method
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), inter_fluc_m_n, lo_inter_fluc_m_n)
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), inter_fluc_p_n, lo_inter_fluc_p_n)
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), cell_fluc_n, lo_cell_fluc_n)
    !$acc declare link(inter_fluc_m_n, inter_fluc_p_n, cell_fluc_n, lo_inter_fluc_m_n, lo_inter_fluc_p_n, lo_cell_fluc_n)
#else
    type(vector_field), allocatable, dimension(:) :: inter_fluc_m_n, lo_inter_fluc_m_n
    type(vector_field), allocatable, dimension(:) :: inter_fluc_p_n, lo_inter_fluc_p_n
    type(vector_field), allocatable, dimension(:) :: cell_fluc_n, lo_cell_fluc_n
    !$acc declare create(inter_fluc_m_n, inter_fluc_p_n, cell_fluc_n, lo_inter_fluc_m_n, lo_inter_fluc_p_n, lo_cell_fluc_n)
#endif
    !> @}

    !> @name The non-conservative terms in the phasic total energy equation
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), phasic_pres_src_n, lo_phasic_pres_src_n)
    !$acc declare link(phasic_pres_src_n, lo_phasic_pres_src_n)
#else
    type(vector_field), allocatable, dimension(:) :: phasic_pres_src_n, lo_phasic_pres_src_n
    !$acc declare create(phasic_pres_src_n, lo_phasic_pres_src_n)
#endif
    !> @}

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), qL_prim, qR_prim, lo_qL_prim, lo_qR_prim)
    !$acc declare link(qL_prim, qR_prim, lo_qL_prim, lo_qR_prim)
#else
    type(vector_field), allocatable, dimension(:) :: qL_prim, qR_prim, lo_qL_prim, lo_qR_prim
    !$acc declare create(qL_prim, qR_prim, lo_qL_prim, lo_qR_prim)
#endif

    type(int_bounds_info) :: iv !< Vector field indical bounds
    !$acc declare create(iv)

    !> @name Indical bounds in the x-, y- and z-directions
    !> @{
    type(int_bounds_info) :: ix, iy, iz
    !$acc declare create(ix, iy, iz)

    type(int_bounds_info) :: is1, is2, is3
    !$acc declare create(is1, is2, is3)

    type(int_bounds_info) :: ixt, iyt, izt
    !$acc declare create(ixt, iyt, izt)

    !> @name Saved fluxes for testing
    !> @{
    type(scalar_field) :: alf_sum
    !> @}
!$acc declare create(alf_sum)

#ifdef CRAY_ACC_WAR

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), blkmod1, blkmod2, alpha1, alpha2, Kterm)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, qR_rsx_vf, qR_rsy_vf, qR_rsz_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf)
    !$acc declare link(blkmod1, blkmod2, alpha1, alpha2, Kterm)
    !$acc declare link(qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, qR_rsx_vf, qR_rsy_vf, qR_rsz_vf)
    !$acc declare link(lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf)
    !$acc declare link(dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf)

#else
    real(kind(0d0)), allocatable, dimension(:, :, :) :: blkmod1, blkmod2, alpha1, alpha2, Kterm
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, qR_rsx_vf, qR_rsy_vf, qR_rsz_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf
    !$acc declare create(blkmod1, blkmod2, alpha1, alpha2, Kterm)
    !$acc declare create(qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, qR_rsx_vf, qR_rsy_vf, qR_rsz_vf)
    !$acc declare create(lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf)
    !$acc declare create(dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf)
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), gamma_min, pres_inf)
    !$acc declare link(gamma_min, pres_inf)
#else
    real(kind(0d0)), allocatable, dimension(:) :: gamma_min, pres_inf
    !$acc declare create(gamma_min, pres_inf)
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), Res)
    !$acc declare link(Res)
#else
    real(kind(0d0)), allocatable, dimension(:, :) :: Res
    !$acc declare create(Res)
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), nbub)
    !$acc declare link(nbub)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :) :: nbub !< Bubble number density
    !$acc declare create(nbub)
#endif

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_rhs_module() ! ---------------------------------

        integer :: i, j, k, l, id !< Generic loop iterators

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================

        !$acc enter data copyin(ix, iy, iz)
        !$acc update device(ix, iy, iz)

        ixt = ix; iyt = iy; izt = iz

        @:ALLOCATE(q_cons_qp%vf(1:sys_size))
        @:ALLOCATE(q_prim_qp%vf(1:sys_size))
        ! Here, note that the momentum, the total energy, and the phasic total energy/internal energy term
        ! are assigned space independently as primitive variables. 
        ! while phasic mass, volume fraction are assigned to point to the same memory location as the conservative variables.

        do l = 1, sys_size
            @:ALLOCATE(q_cons_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        do l = mom_idx%beg, E_idx
            @:ALLOCATE(q_prim_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        do l = adv_idx%end + 1, sys_size
            @:ALLOCATE(q_prim_qp%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end do

        @:ACC_SETUP_VFs(q_cons_qp, q_prim_qp)

        do l = 1, cont_idx%end
            q_prim_qp%vf(l)%sf => q_cons_qp%vf(l)%sf
            !$acc enter data copyin(q_prim_qp%vf(l)%sf)
            !$acc enter data attach(q_prim_qp%vf(l)%sf)
        end do

        do l = adv_idx%beg, adv_idx%end
            q_prim_qp%vf(l)%sf => q_cons_qp%vf(l)%sf
            !$acc enter data copyin(q_prim_qp%vf(l)%sf)
            !$acc enter data attach(q_prim_qp%vf(l)%sf)
        end do

        ! ==================================================================

        if (qbmm) then
            @:ALLOCATE_GLOBAL(mom_sp(1:nmomsp), mom_3d(0:2, 0:2, nb))

            do i = 0, 2
                do j = 0, 2
                    do k = 1, nb
                        @:ALLOCATE(mom_3d(i, j, k)%sf( &
                                      & ix%beg:ix%end, &
                                      & iy%beg:iy%end, &
                                      & iz%beg:iz%end))
                        @:ACC_SETUP_SFs(mom_3d(i, j, k))
                    end do
                end do
            end do

            do i = 1, nmomsp
                @:ALLOCATE(mom_sp(i)%sf( &
                        & ix%beg:ix%end, &
                        & iy%beg:iy%end, &
                        & iz%beg:iz%end))
                @:ACC_SETUP_SFs(mom_sp(i))
            end do
        end if

        ! Allocation/Association of qK_cons_n and qK_prim_n ==========
        @:ALLOCATE_GLOBAL(qL_prim(1:num_dims))
        @:ALLOCATE_GLOBAL(qR_prim(1:num_dims))

        do i = 1, num_dims
            @:ALLOCATE(qL_prim(i)%vf(1:sys_size))
            @:ALLOCATE(qR_prim(i)%vf(1:sys_size))
            do l = mom_idx%beg, mom_idx%end
                @:ALLOCATE(qL_prim(i)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
                @:ALLOCATE(qR_prim(i)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
            end do
            @:ACC_SETUP_VFs(qL_prim(i), qR_prim(i))
        end do

        @:ALLOCATE_GLOBAL(lo_qL_prim(1:num_dims))
        @:ALLOCATE_GLOBAL(lo_qR_prim(1:num_dims))

        do i = 1, num_dims
            @:ALLOCATE(lo_qL_prim(i)%vf(1:sys_size))
            @:ALLOCATE(lo_qR_prim(i)%vf(1:sys_size))
            do l = mom_idx%beg, mom_idx%end
                @:ALLOCATE(lo_qL_prim(i)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
                @:ALLOCATE(lo_qR_prim(i)%vf(l)%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
            end do
            @:ACC_SETUP_VFs(lo_qL_prim(i), lo_qR_prim(i))
        end do

        if (mpp_lim .and. bubbles) then
            @:ALLOCATE(alf_sum%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end if
        ! END: Allocation/Association of qK_cons_n and qK_prim_n ======

        @:ALLOCATE_GLOBAL(qL_rsx_vf(ix%beg:ix%end, &
            iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(qR_rsx_vf(ix%beg:ix%end, &
            iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))

        @:ALLOCATE_GLOBAL(lo_qL_rsx_vf(ix%beg:ix%end, &
            iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(lo_qR_rsx_vf(ix%beg:ix%end, &
            iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))

        if (n > 0) then

            @:ALLOCATE_GLOBAL(qL_rsy_vf(iy%beg:iy%end, &
                ix%beg:ix%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE_GLOBAL(qR_rsy_vf(iy%beg:iy%end, &
                ix%beg:ix%end, iz%beg:iz%end, 1:sys_size))

            @:ALLOCATE_GLOBAL(lo_qL_rsy_vf(iy%beg:iy%end, &
                ix%beg:ix%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE_GLOBAL(lo_qR_rsy_vf(iy%beg:iy%end, &
                ix%beg:ix%end, iz%beg:iz%end, 1:sys_size))
        else
            @:ALLOCATE_GLOBAL(qL_rsy_vf(ix%beg:ix%end, &
                iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE_GLOBAL(qR_rsy_vf(ix%beg:ix%end, &
                iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))

            @:ALLOCATE_GLOBAL(lo_qL_rsy_vf(ix%beg:ix%end, &
                iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE_GLOBAL(lo_qR_rsy_vf(ix%beg:ix%end, &
                iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
        end if

        if (p > 0) then
            @:ALLOCATE_GLOBAL(qL_rsz_vf(iz%beg:iz%end, &
                iy%beg:iy%end, ix%beg:ix%end, 1:sys_size))
            @:ALLOCATE_GLOBAL(qR_rsz_vf(iz%beg:iz%end, &
                iy%beg:iy%end, ix%beg:ix%end, 1:sys_size))

            @:ALLOCATE_GLOBAL(lo_qL_rsz_vf(iz%beg:iz%end, &
                iy%beg:iy%end, ix%beg:ix%end, 1:sys_size))
            @:ALLOCATE_GLOBAL(lo_qR_rsz_vf(iz%beg:iz%end, &
                iy%beg:iy%end, ix%beg:ix%end, 1:sys_size))
        else
            @:ALLOCATE_GLOBAL(qL_rsz_vf(ix%beg:ix%end, &
                iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE_GLOBAL(qR_rsz_vf(ix%beg:ix%end, &
                iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))

            @:ALLOCATE_GLOBAL(lo_qL_rsz_vf(ix%beg:ix%end, &
                iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))
            @:ALLOCATE_GLOBAL(lo_qR_rsz_vf(ix%beg:ix%end, &
                iy%beg:iy%end, iz%beg:iz%end, 1:sys_size))

        end if

        ! Allocation of dq_prim_ds_qp ======================================

        @:ALLOCATE_GLOBAL(dq_prim_dx_qp(1:1))
        @:ALLOCATE_GLOBAL(dq_prim_dy_qp(1:1))
        @:ALLOCATE_GLOBAL(dq_prim_dz_qp(1:1))

        if (any(Re_size > 0)) then
            @:ALLOCATE(dq_prim_dx_qp(1)%vf(1:sys_size))
            @:ALLOCATE(dq_prim_dy_qp(1)%vf(1:sys_size))
            @:ALLOCATE(dq_prim_dz_qp(1)%vf(1:sys_size))
            if (any(Re_size > 0)) then

                do l = mom_idx%beg, mom_idx%end
                    @:ALLOCATE(dq_prim_dx_qp(1)%vf(l)%sf( &
                              & ix%beg:ix%end, &
                              & iy%beg:iy%end, &
                              & iz%beg:iz%end))
                end do

                @:ACC_SETUP_VFs(dq_prim_dx_qp(1))

                if (n > 0) then

                    do l = mom_idx%beg, mom_idx%end
                        @:ALLOCATE(dq_prim_dy_qp(1)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                    end do

                    @:ACC_SETUP_VFs(dq_prim_dy_qp(1))

                    if (p > 0) then

                        do l = mom_idx%beg, mom_idx%end
                            @:ALLOCATE(dq_prim_dz_qp(1)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                        end do
                        @:ACC_SETUP_VFs(dq_prim_dz_qp(1))
                    end if

                end if

            end if

        end if
        ! END: Allocation of dq_prim_ds_qp =================================

        ! Allocation/Association of dqK_prim_ds_n =======================
        @:ALLOCATE_GLOBAL(dqL_prim_dx_n(1:num_dims))
        @:ALLOCATE_GLOBAL(dqL_prim_dy_n(1:num_dims))
        @:ALLOCATE_GLOBAL(dqL_prim_dz_n(1:num_dims))
        @:ALLOCATE_GLOBAL(dqR_prim_dx_n(1:num_dims))
        @:ALLOCATE_GLOBAL(dqR_prim_dy_n(1:num_dims))
        @:ALLOCATE_GLOBAL(dqR_prim_dz_n(1:num_dims))

        if (any(Re_size > 0)) then
            do i = 1, num_dims
                @:ALLOCATE(dqL_prim_dx_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqL_prim_dy_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqL_prim_dz_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqR_prim_dx_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqR_prim_dy_n(i)%vf(1:sys_size))
                @:ALLOCATE(dqR_prim_dz_n(i)%vf(1:sys_size))

                if (any(Re_size > 0)) then

                    do l = mom_idx%beg, mom_idx%end
                        @:ALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                        @:ALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                    end do

                    if (n > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:ALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                            @:ALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                        end do
                    end if

                    if (p > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:ALLOCATE(dqL_prim_dz_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                            @:ALLOCATE(dqR_prim_dz_n(i)%vf(l)%sf( &
                                     & ix%beg:ix%end, &
                                     & iy%beg:iy%end, &
                                     & iz%beg:iz%end))
                        end do
                    end if

                end if

                @:ACC_SETUP_VFs(dqL_prim_dx_n(i), dqL_prim_dy_n(i), dqL_prim_dz_n(i))
                @:ACC_SETUP_VFs(dqR_prim_dx_n(i), dqR_prim_dy_n(i), dqR_prim_dz_n(i))
            end do
        end if
        ! END: Allocation/Association of d K_prim_ds_n ==================

        if (any(Re_size > 0)) then
            if (weno_Re_flux) then
                @:ALLOCATE_GLOBAL(dqL_rsx_vf(ix%beg:ix%end, &
                    iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                @:ALLOCATE_GLOBAL(dqR_rsx_vf(ix%beg:ix%end, &
                    iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))

                if (n > 0) then

                    @:ALLOCATE_GLOBAL(dqL_rsy_vf(iy%beg:iy%end, &
                        ix%beg:ix%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE_GLOBAL(dqR_rsy_vf(iy%beg:iy%end, &
                        ix%beg:ix%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                else
                    @:ALLOCATE_GLOBAL(dqL_rsy_vf(ix%beg:ix%end, &
                        iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE_GLOBAL(dqR_rsy_vf(ix%beg:ix%end, &
                        iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))

                end if

                if (p > 0) then
                    @:ALLOCATE_GLOBAL(dqL_rsz_vf(iz%beg:iz%end, &
                        iy%beg:iy%end, ix%beg:ix%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE_GLOBAL(dqR_rsz_vf(iz%beg:iz%end, &
                        iy%beg:iy%end, ix%beg:ix%end, mom_idx%beg:mom_idx%end))
                else
                    @:ALLOCATE_GLOBAL(dqL_rsz_vf(ix%beg:ix%end, &
                        iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))
                    @:ALLOCATE_GLOBAL(dqR_rsz_vf(ix%beg:ix%end, &
                        iy%beg:iy%end, iz%beg:iz%end, mom_idx%beg:mom_idx%end))

                end if
            end if
        end if

        ! ==================================================================

        ! Allocation of gm_alphaK_n =====================================
        @:ALLOCATE_GLOBAL(gm_alphaL_n(1:num_dims))
        @:ALLOCATE_GLOBAL(gm_alphaR_n(1:num_dims))
        ! ==================================================================

        ! Allocation/Association of flux_n, flux_src_n, and flux_gsrc_n ===
        ! and also the interfacial fluctuation terms inter_fluc_m_n and inter_fluc_p_n,
        ! and the cell-fluctuation terms cell_fluc_n, the phasic pressure source terms phasic_pres_src_n
        @:ALLOCATE_GLOBAL(cell_fluc_n(1:num_dims))
        @:ALLOCATE_GLOBAL(inter_fluc_m_n(1:num_dims))
        @:ALLOCATE_GLOBAL(inter_fluc_p_n(1:num_dims))
        @:ALLOCATE_GLOBAL(lo_cell_fluc_n(1:num_dims))
        @:ALLOCATE_GLOBAL(lo_inter_fluc_m_n(1:num_dims))
        @:ALLOCATE_GLOBAL(lo_inter_fluc_p_n(1:num_dims))
        @:ALLOCATE_GLOBAL(phasic_pres_src_n(1:num_fluids))
        @:ALLOCATE_GLOBAL(lo_phasic_pres_src_n(1:num_fluids))
        @:ALLOCATE_GLOBAL(flux_n(1:num_dims))
        @:ALLOCATE_GLOBAL(flux_src_n(1:num_dims))
        @:ALLOCATE_GLOBAL(flux_gsrc_n(1:num_dims))
        @:ALLOCATE_GLOBAL(lo_flux_n(1:num_dims))
        @:ALLOCATE_GLOBAL(lo_flux_src_n(1:num_dims))
        @:ALLOCATE_GLOBAL(lo_flux_gsrc_n(1:num_dims))

        do i = 1, num_dims

            @:ALLOCATE(cell_fluc_n(i)%vf(1:sys_size))
            @:ALLOCATE(inter_fluc_m_n(i)%vf(1:sys_size))
            @:ALLOCATE(inter_fluc_p_n(i)%vf(1:sys_size))
            @:ALLOCATE(lo_cell_fluc_n(i)%vf(1:sys_size))
            @:ALLOCATE(lo_inter_fluc_m_n(i)%vf(1:sys_size))
            @:ALLOCATE(lo_inter_fluc_p_n(i)%vf(1:sys_size))
            @:ALLOCATE(phasic_pres_src_n(i)%vf(1:sys_size))
            @:ALLOCATE(lo_phasic_pres_src_n(i)%vf(1:sys_size))
            @:ALLOCATE(flux_n(i)%vf(1:sys_size))
            @:ALLOCATE(flux_src_n(i)%vf(1:sys_size))
            @:ALLOCATE(flux_gsrc_n(i)%vf(1:sys_size))
            @:ALLOCATE(lo_flux_n(i)%vf(1:sys_size))
            @:ALLOCATE(lo_flux_src_n(i)%vf(1:sys_size))
            @:ALLOCATE(lo_flux_gsrc_n(i)%vf(1:sys_size))

            if (i == 1) then
                do l = 1, sys_size
                    @:ALLOCATE(cell_fluc_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                    @:ALLOCATE(inter_fluc_m_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                    @:ALLOCATE(inter_fluc_p_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                    @:ALLOCATE(lo_cell_fluc_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                    @:ALLOCATE(lo_inter_fluc_m_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                    @:ALLOCATE(lo_inter_fluc_p_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                    @:ALLOCATE(flux_n(i)%vf(l)%sf( &
                             & ix%beg:ix%end, &
                             & iy%beg:iy%end, &
                             & iz%beg:iz%end))
                    @:ALLOCATE(flux_gsrc_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                    @:ALLOCATE(lo_flux_n(i)%vf(l)%sf( &
                             & ix%beg:ix%end, &
                             & iy%beg:iy%end, &
                             & iz%beg:iz%end))
                    @:ALLOCATE(lo_flux_gsrc_n(i)%vf(l)%sf( &
                            & ix%beg:ix%end, &
                            & iy%beg:iy%end, &
                            & iz%beg:iz%end))
                end do
                ! The phasic pressure source terms are only needed in phasic total energy term
                do l = intxb, sys_size
                    @:ALLOCATE(phasic_pres_src_n(i)%vf(l)%sf( &
                             & ix%beg:ix%end, &
                             & iy%beg:iy%end, &
                             & iz%beg:iz%end))
                    @:ALLOCATE(lo_phasic_pres_src_n(i)%vf(l)%sf( &
                             & ix%beg:ix%end, &
                             & iy%beg:iy%end, &
                             & iz%beg:iz%end))
                end do

                if (any(Re_size > 0)) then
                    do l = mom_idx%beg, E_idx
                        @:ALLOCATE(flux_src_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                        @:ALLOCATE(lo_flux_src_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                    end do
                end if

                @:ALLOCATE(flux_src_n(i)%vf(adv_idx%beg)%sf( &
                         & ix%beg:ix%end, &
                         & iy%beg:iy%end, &
                         & iz%beg:iz%end))

                @:ALLOCATE(lo_flux_src_n(i)%vf(adv_idx%beg)%sf( &
                         & ix%beg:ix%end, &
                         & iy%beg:iy%end, &
                         & iz%beg:iz%end))

                if (riemann_solver == 1) then
                    do l = adv_idx%beg + 1, adv_idx%end
                        @:ALLOCATE(flux_src_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                        @:ALLOCATE(lo_flux_src_n(i)%vf(l)%sf( &
                                 & ix%beg:ix%end, &
                                 & iy%beg:iy%end, &
                                 & iz%beg:iz%end))
                    end do
                end if
            else
                do l = 1, sys_size
                    @:ALLOCATE(flux_gsrc_n(i)%vf(l)%sf( &
                        ix%beg:ix%end, &
                        iy%beg:iy%end, &
                        iz%beg:iz%end))
                    @:ALLOCATE(lo_flux_gsrc_n(i)%vf(l)%sf( &
                        ix%beg:ix%end, &
                        iy%beg:iy%end, &
                        iz%beg:iz%end))
                end do
            end if

            @:ACC_SETUP_VFs(cell_fluc_n(i), inter_fluc_m_n(i), inter_fluc_p_n(i), &
                        lo_cell_fluc_n(i), lo_inter_fluc_m_n(i), lo_inter_fluc_p_n(i), &
                        flux_n(i), flux_src_n(i), flux_gsrc_n(i), phasic_pres_src_n(i), &
                        lo_flux_n(i), lo_flux_src_n(i), lo_flux_gsrc_n(i), lo_phasic_pres_src_n(i))

            ! We can see that the i-direction's flux is the same as the first flux 
            ! since they share the same memory location. Confusing! Is it right?
            if (i == 1) then
                if (riemann_solver /= 1) then
                    do l = adv_idx%beg + 1, adv_idx%end
                        flux_src_n(i)%vf(l)%sf => flux_src_n(i)%vf(adv_idx%beg)%sf
                        lo_flux_src_n(i)%vf(l)%sf => lo_flux_src_n(i)%vf(adv_idx%beg)%sf

                        !$acc enter data attach(flux_src_n(i)%vf(l)%sf, lo_flux_src_n(i)%vf(l)%sf)
                    end do
                end if
            else
                do l = 1, sys_size
                    cell_fluc_n(i)%vf(l)%sf => cell_fluc_n(1)%vf(l)%sf
                    inter_fluc_m_n(i)%vf(l)%sf => inter_fluc_m_n(1)%vf(l)%sf
                    inter_fluc_p_n(i)%vf(l)%sf => inter_fluc_p_n(1)%vf(l)%sf
                    lo_cell_fluc_n(i)%vf(l)%sf => lo_cell_fluc_n(1)%vf(l)%sf
                    lo_inter_fluc_m_n(i)%vf(l)%sf => lo_inter_fluc_m_n(1)%vf(l)%sf
                    lo_inter_fluc_p_n(i)%vf(l)%sf => lo_inter_fluc_p_n(1)%vf(l)%sf
                    flux_n(i)%vf(l)%sf => flux_n(1)%vf(l)%sf
                    flux_src_n(i)%vf(l)%sf => flux_src_n(1)%vf(l)%sf
                    phasic_pres_src_n(i)%vf(l)%sf => phasic_pres_src_n(1)%vf(l)%sf
                    lo_phasic_pres_src_n(i)%vf(l)%sf => lo_phasic_pres_src_n(1)%vf(l)%sf
                    lo_flux_n(i)%vf(l)%sf => lo_flux_n(1)%vf(l)%sf
                    lo_flux_src_n(i)%vf(l)%sf => lo_flux_src_n(1)%vf(l)%sf

                    !$acc enter data attach(cell_fluc_n(i)%vf(l)%sf,inter_fluc_m_n(i)%vf(l)%sf,inter_fluc_p_n(i)%vf(l)%sf,&
                    !$acc  lo_cell_fluc_n(i)%vf(l)%sf, lo_inter_fluc_m_n(i)%vf(l)%sf, lo_inter_fluc_p_n(i)%vf(l)%sf, &
                    !&acc  flux_n(i)%vf(l)%sf,flux_src_n(i)%vf(l)%sf, phasic_pres_src_n(i)%vf(l)%sf,&
                    !&acc  lo_flux_n(i)%vf(l)%sf,lo_flux_src_n(i)%vf(l)%sf, lo_phasic_pres_src_n(i)%vf(l)%sf)
                end do
            end if
        end do

        ! END: Allocation/Association of flux_n, flux_src_n, and flux_gsrc_n ===

        if (alt_soundspeed) then
            @:ALLOCATE_GLOBAL(blkmod1(0:m, 0:n, 0:p), blkmod2(0:m, 0:n, 0:p), alpha1(0:m, 0:n, 0:p), alpha2(0:m, 0:n, 0:p), Kterm(0:m, 0:n, 0:p))
        end if

        @:ALLOCATE_GLOBAL(gamma_min(1:num_fluids), pres_inf(1:num_fluids))

        do i = 1, num_fluids
            gamma_min(i) = 1d0/fluid_pp(i)%gamma + 1d0
            pres_inf(i) = fluid_pp(i)%pi_inf/(1d0 + fluid_pp(i)%gamma)
        end do
        !$acc update device(gamma_min, pres_inf)

        if (any(Re_size > 0)) then
            @:ALLOCATE_GLOBAL(Res(1:2, 1:maxval(Re_size)))
        end if

        if (any(Re_size > 0)) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            !$acc update device(Res, Re_idx, Re_size)
        end if

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        if (riemann_solver == 1) then
            s_riemann_solver => s_hll_riemann_solver
        elseif (riemann_solver == 2) then
            s_riemann_solver => s_hllc_riemann_solver
        end if

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables
        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        else if (bubbles) then          ! Volume fraction for bubbles
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

        !$acc parallel loop collapse(4) gang vector default(present)
        do id = 1, num_dims
            do i = 1, sys_size
                do l = startz, p - startz
                    do k = starty, n - starty
                        do j = startx, m - startx
                            flux_gsrc_n(id)%vf(i)%sf(j, k, l) = 0d0
                            lo_flux_gsrc_n(id)%vf(i)%sf(j, k, l) = 0d0
                        end do
                    end do
                end do
            end do
        end do

        if (bubbles) then
            @:ALLOCATE_GLOBAL(nbub(0:m, 0:n, 0:p))
        end if

    end subroutine s_initialize_rhs_module ! -------------------------------

    subroutine s_compute_rhs(q_cons_vf, q_prim_vf, rhs_vf, pb, rhs_pb, mv, rhs_mv, t_step, MOOD_button, lo_rhs_vf) ! -------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf, lo_rhs_vf
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(INOUT) :: pb, mv
        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:, 1:), intent(INOUT) :: rhs_pb, rhs_mv
        integer, intent(IN) :: t_step
        logical, intent(IN) :: MOOD_button

        real(kind(0d0)) :: top, bottom  !< Numerator and denominator when evaluating flux limiter function
        real(kind(0d0)), dimension(num_fluids) :: myalpha_rho, myalpha

        real(kind(0d0)) :: tmp1, tmp2, tmp3, tmp4, &
                           c_gas, c_liquid, &
                           Cpbw, Cpinf, Cpinf_dot, &
                           myH, myHdot, rddot, alf_gas

        real(kind(0d0)) :: n_tait, B_tait, angle, angle_z

        real(kind(0d0)), dimension(nb) :: Rtmp, Vtmp
        real(kind(0d0)) :: myR, myV, alf, myP, myRho, R2Vav
        real(kind(0d0)), dimension(0:m, 0:n, 0:p) :: nbub
        integer :: ndirs

        real(kind(0d0)) :: mytime, sound
        real(kind(0d0)) :: start, finish
        real(kind(0d0)) :: s2, const_sos, s1

        integer :: i, j, k, l, q, ii, id !< Generic loop iterators
        integer :: term_index
        logical :: wave_propagation

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================

        !$acc update device(ix, iy, iz)

        ! Association/Population of Working Variables ======================
        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        q_cons_qp%vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        ! ==================================================================

        ! Converting Conservative to Primitive Variables ==================

        if (mpp_lim .and. bubbles) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end
                        alf_sum%sf(j, k, l) = 0d0
                        !$acc loop seq
                        do i = advxb, advxe - 1
                            alf_sum%sf(j, k, l) = alf_sum%sf(j, k, l) + q_cons_qp%vf(i)%sf(j, k, l)
                        end do
                        !$acc loop seq
                        do i = advxb, advxe - 1
                            q_cons_qp%vf(i)%sf(j, k, l) = q_cons_qp%vf(i)%sf(j, k, l)*(1.d0 - q_cons_qp%vf(alf_idx)%sf(j, k, l)) &
                                                          /alf_sum%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        end if

        call nvtxStartRange("RHS-CONVERT")
        call s_convert_conservative_to_primitive_variables( &
            q_cons_qp%vf, &
            q_prim_qp%vf, &
            gm_alpha_qp%vf, &
            ix, iy, iz)
        call nvtxEndRange

        call nvtxStartRange("RHS-MPI")
        call s_populate_primitive_variables_buffers(q_prim_qp%vf, pb, mv)
        call nvtxEndRange

        if (t_step == t_step_stop) return
        ! ==================================================================

        if (qbmm) call s_mom_inv(q_cons_qp%vf, q_prim_qp%vf, mom_sp, mom_3d, pb, rhs_pb, mv, rhs_mv, ix, iy, iz, nbub)

        call nvtxStartRange("Viscous")
        if (any(Re_size > 0)) call s_get_viscous(qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                                                 dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, &
                                                 qL_prim, &
                                                 qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                                                 dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n, &
                                                 qR_prim, &
                                                 q_prim_qp, &
                                                 dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, &
                                                 ix, iy, iz)
        call nvtxEndRange()

        ! Dimensional Splitting Loop =======================================

        do id = 1, num_dims

            ! Configuring Coordinate Direction Indexes ======================
            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

            if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            ! ===============================================================
            ! Reconstructing Primitive/Conservative Variables ===============

            if (all(Re_size == 0)) then
                iv%beg = 1; iv%end = sys_size
                call nvtxStartRange("RHS-WENO")
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(1:sys_size), &
                    qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                    qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                    id, .false.)
                if (MOOD_button) then
                    call s_reconstruct_cell_boundary_values( &
                        q_prim_qp%vf(1:sys_size), &
                        lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, &
                        lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf, &
                        id, .true.)
                end if
                call nvtxEndRange
            else
                call nvtxStartRange("RHS-WENO")
                iv%beg = 1; iv%end = contxe
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                    qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                    id, .false.)
                if (MOOD_button) then
                    call s_reconstruct_cell_boundary_values( &
                        q_prim_qp%vf(iv%beg:iv%end), &
                        lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, &
                        lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf, &
                        id, .true.)
                end if

                iv%beg = E_idx; iv%end = E_idx
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                    qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                    id, .false.)
                if (MOOD_button) then
                    call s_reconstruct_cell_boundary_values( &
                        q_prim_qp%vf(iv%beg:iv%end), &
                        lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, &
                        lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf, &
                        id, .true.)
                end if

                iv%beg = advxb; iv%end = advxe
                call s_reconstruct_cell_boundary_values( &
                    q_prim_qp%vf(iv%beg:iv%end), &
                    qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                    qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                    id, .false.)
                if (MOOD_button) then
                    call s_reconstruct_cell_boundary_values( &
                        q_prim_qp%vf(iv%beg:iv%end), &
                        lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, &
                        lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf, &
                        id, .true.)
                end if

                if (bubbles) then
                    iv%beg = bubxb; iv%end = bubxe
                    call s_reconstruct_cell_boundary_values( &
                        q_prim_qp%vf(iv%beg:iv%end), &
                        qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                        qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                        id, .false.)
                    if (MOOD_button) then
                        call s_reconstruct_cell_boundary_values( &
                            q_prim_qp%vf(iv%beg:iv%end), &
                            lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, &
                            lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf, &
                            id, .true.)
                    end if
                end if

                if (adv_n) then
                    iv%beg = n_idx; iv%end = n_idx
                    call s_reconstruct_cell_boundary_values( &
                        q_prim_qp%vf(iv%beg:iv%end), &
                        qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                        qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                        id, .false.)
                    if (MOOD_button) then
                        call s_reconstruct_cell_boundary_values( &
                            q_prim_qp%vf(iv%beg:iv%end), &
                            lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, &
                            lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf, &
                            id, .true.)
                    end if
                end if

                iv%beg = mom_idx%beg; iv%end = mom_idx%end
                if (weno_Re_flux) then
                    call s_reconstruct_cell_boundary_values_visc_deriv( &
                        dq_prim_dx_qp(1)%vf(iv%beg:iv%end), &
                        dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, &
                        dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf, &
                        id, dqL_prim_dx_n(id)%vf(iv%beg:iv%end), dqR_prim_dx_n(id)%vf(iv%beg:iv%end), &
                        ix, iy, iz)
                    if (n > 0) then
                        call s_reconstruct_cell_boundary_values_visc_deriv( &
                            dq_prim_dy_qp(1)%vf(iv%beg:iv%end), &
                            dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, &
                            dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf, &
                            id, dqL_prim_dy_n(id)%vf(iv%beg:iv%end), dqR_prim_dy_n(id)%vf(iv%beg:iv%end), &
                            ix, iy, iz)
                        if (p > 0) then
                            call s_reconstruct_cell_boundary_values_visc_deriv( &
                                dq_prim_dz_qp(1)%vf(iv%beg:iv%end), &
                                dqL_rsx_vf, dqL_rsy_vf, dqL_rsz_vf, &
                                dqR_rsx_vf, dqR_rsy_vf, dqR_rsz_vf, &
                                id, dqL_prim_dz_n(id)%vf(iv%beg:iv%end), dqR_prim_dz_n(id)%vf(iv%beg:iv%end), &
                                ix, iy, iz)
                        end if
                    end if
                end if
                call nvtxEndRange
            end if

            ! Configuring Coordinate Direction Indexes ======================
            if (id == 1) then
                ix%beg = -1; iy%beg = 0; iz%beg = 0
            elseif (id == 2) then
                ix%beg = 0; iy%beg = -1; iz%beg = 0
            else
                ix%beg = 0; iy%beg = 0; iz%beg = -1
            end if
            ix%end = m; iy%end = n; iz%end = p
            ! ===============================================================

            ! Computing Riemann Solver Flux and Source Flux =================
            wave_propagation = .false.
            if (MOOD_button) then
                call s_riemann_solver(lo_qR_rsx_vf, lo_qR_rsy_vf, lo_qR_rsz_vf, &
                                  dqR_prim_dx_n(id)%vf, &
                                  dqR_prim_dy_n(id)%vf, &
                                  dqR_prim_dz_n(id)%vf, &
                                  lo_qR_prim(id)%vf, &
                                  lo_qL_rsx_vf, lo_qL_rsy_vf, lo_qL_rsz_vf, &
                                  dqL_prim_dx_n(id)%vf, &
                                  dqL_prim_dy_n(id)%vf, &
                                  dqL_prim_dz_n(id)%vf, &
                                  lo_qL_prim(id)%vf, &
                                  q_prim_qp%vf, &
                                  lo_flux_n(id)%vf, &
                                  lo_flux_src_n(id)%vf, &
                                  lo_flux_gsrc_n(id)%vf, &
                                  id, ix, iy, iz, &
                                  lo_inter_fluc_m_n(id)%vf, &
                                  lo_inter_fluc_p_n(id)%vf, &
                                  lo_cell_fluc_n(id)%vf, &
                                  lo_phasic_pres_src_n(id)%vf, &
                                  wave_propagation)
            end if
            call s_riemann_solver(qR_rsx_vf, qR_rsy_vf, qR_rsz_vf, &
                                  dqR_prim_dx_n(id)%vf, &
                                  dqR_prim_dy_n(id)%vf, &
                                  dqR_prim_dz_n(id)%vf, &
                                  qR_prim(id)%vf, &
                                  qL_rsx_vf, qL_rsy_vf, qL_rsz_vf, &
                                  dqL_prim_dx_n(id)%vf, &
                                  dqL_prim_dy_n(id)%vf, &
                                  dqL_prim_dz_n(id)%vf, &
                                  qL_prim(id)%vf, &
                                  q_prim_qp%vf, &
                                  flux_n(id)%vf, &
                                  flux_src_n(id)%vf, &
                                  flux_gsrc_n(id)%vf, &
                                  id, ix, iy, iz, &
                                  inter_fluc_m_n(id)%vf, &
                                  inter_fluc_p_n(id)%vf, &
                                  cell_fluc_n(id)%vf, &
                                  phasic_pres_src_n(id)%vf, &
                                  wave_propagation)
            call nvtxEndRange

            ! ===============================================================

            ! call nvtxStartRange("RHS_Flux_Add")
            ! call nvtxEndRange

            ! Additional physics and source terms ==============================

            ! RHS addition for advection source
            call nvtxStartRange("RHS_advection_source")
            call s_compute_advection_source_term(id, &
                                                 rhs_vf, &
                                                 q_cons_qp, &
                                                 q_prim_qp, &
                                                 flux_src_n(id), &
                                                 lo_rhs_vf, &
                                                 MOOD_button)
            call nvtxEndRange()

            ! RHS additions for hypoelasticity
            call nvtxStartRange("RHS_Hypoelasticity")
            if (hypoelasticity) call s_compute_hypoelastic_rhs(id, &
                                                               q_prim_qp%vf, &
                                                               rhs_vf)
            call nvtxEndRange

            ! RHS additions for viscosity
            call nvtxStartRange("RHS_viscous")
            if (any(Re_size > 0d0)) call s_compute_viscous_rhs(id, &
                                                               q_prim_qp%vf, &
                                                               rhs_vf, &
                                                               flux_src_n(id)%vf, &
                                                               dq_prim_dx_qp(1)%vf, &
                                                               dq_prim_dy_qp(1)%vf, &
                                                               dq_prim_dz_qp(1)%vf, &
                                                               ixt, iyt, izt)
            call nvtxEndRange

            ! RHS additions for sub-grid bubbles
            call nvtxStartRange("RHS_bubbles")
            if (bubbles) call s_compute_bubbles_rhs(id, &
                                                    q_prim_qp%vf)
            call nvtxEndRange

            ! RHS additions for qbmm bubbles
            call nvtxStartRange("RHS_qbmm")
            if (qbmm) call s_compute_qbmm_rhs(id, &
                                              q_cons_qp%vf, &
                                              q_prim_qp%vf, &
                                              rhs_vf, &
                                              flux_n(id)%vf, &
                                              pb, &
                                              rhs_pb, &
                                              mv, &
                                              rhs_mv)
            call nvtxEndRange
            ! END: Additional physics and source terms =========================

        end do

        if (ib) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (ib_markers%sf(j, k, l) /= 0) then
                            do i = 1, sys_size
                                rhs_vf(i)%sf(j, k, l) = 0d0
                            end do
                        end if
                    end do
                end do
            end do
        end if

        ! END: Dimensional Splitting Loop =================================

        ! Additional Physics and Source Temrs ==================================
        ! Additions for monopole
        call nvtxStartRange("RHS_monopole")
        if (monopole) call s_monopole_calculations(q_cons_qp%vf(1:sys_size), &
                                                   q_prim_qp%vf(1:sys_size), &
                                                   t_step, &
                                                   num_dims, &
                                                   rhs_vf)
        call nvtxEndRange

        ! Add bubles source term
        call nvtxStartRange("RHS_bubbles")
        if (bubbles .and. (.not. adap_dt) .and. (.not. qbmm)) call s_compute_bubble_source( &
            q_cons_qp%vf(1:sys_size), &
            q_prim_qp%vf(1:sys_size), &
            t_step, &
            rhs_vf)
        call nvtxEndRange
        ! END: Additional pphysics and source terms ============================

        if (run_time_info .or. probe_wrt .or. ib) then

            ix%beg = -buff_size; iy%beg = 0; iz%beg = 0
            if (n > 0) iy%beg = -buff_size; 
            if (p > 0) iz%beg = -buff_size; 
            ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
            !$acc update device(ix, iy, iz)

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = iz%beg, iz%end
                    do k = iy%beg, iy%end
                        do j = ix%beg, ix%end
                            q_prim_vf(i)%sf(j, k, l) = q_prim_qp%vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do

        end if

        ! ==================================================================

    end subroutine s_compute_rhs ! -----------------------------------------

    subroutine s_compute_advection_source_term(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, &
                                                    lo_rhs_vf, MOOD_button)

        type(vector_field), intent(INOUT) :: q_cons_vf
        type(vector_field), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf, lo_rhs_vf
        type(vector_field), intent(INOUT) :: flux_src_n_vf
        logical, intent(IN) :: MOOD_button

        integer, intent(in) :: idir
        integer :: i, j, k, l, q
        logical :: wave_propagation, nonconservative_term_energy, nonconservative_term_advection
        real(kind(0d0)) :: sum_nonconservative_term, rho, lo_sum_nonconservative_term
        real(kind(0d0)), dimension(1:num_fluids) :: Y_mf

        wave_propagation = .false.
        nonconservative_term_energy = .true.
        nonconservative_term_advection = .true.
        if (wave_propagation) then
            nonconservative_term_energy = .false.
            nonconservative_term_advection = .false.
        end if

        if (alt_soundspeed) then
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        blkmod1(j, k, l) = ((gammas(1) + 1d0)*q_prim_vf%vf(E_idx)%sf(j, k, l) + &
                                            pi_infs(1))/gammas(1)
                        blkmod2(j, k, l) = ((gammas(2) + 1d0)*q_prim_vf%vf(E_idx)%sf(j, k, l) + &
                                            pi_infs(2))/gammas(2)
                        alpha1(j, k, l) = q_cons_vf%vf(advxb)%sf(j, k, l)

                        if (bubbles) then
                            alpha2(j, k, l) = q_cons_vf%vf(alf_idx - 1)%sf(j, k, l)
                        else
                            alpha2(j, k, l) = q_cons_vf%vf(advxe)%sf(j, k, l)
                        end if

                        Kterm(j, k, l) = alpha1(j, k, l)*alpha2(j, k, l)*(blkmod2(j, k, l) - blkmod1(j, k, l))/ &
                                         (alpha1(j, k, l)*blkmod2(j, k, l) + alpha2(j, k, l)*blkmod1(j, k, l))
                    end do
                end do
            end do
        end if

        if (idir == 1) then

            if (bc_x%beg <= -5 .and. bc_x%beg >= -13) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, &
                           flux_src_n(idir)%vf, idir, -1, ix, iy, iz)
                if (MOOD_button) then
                    call s_cbc(q_prim_vf%vf, lo_flux_n(idir)%vf, &
                               lo_flux_src_n(idir)%vf, idir, -1, ix, iy, iz)
                end if
            end if

            if (bc_x%end <= -5 .and. bc_x%end >= -13) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, &
                           flux_src_n(idir)%vf, idir, 1, ix, iy, iz)
                if (MOOD_button) then
                    call s_cbc(q_prim_vf%vf, lo_flux_n(idir)%vf, &
                               lo_flux_src_n(idir)%vf, idir, 1, ix, iy, iz)
                end if
            end if

            !$acc parallel loop collapse(4) gang vector default(present)
            do j = 1, sys_size
                do q = 0, p
                    do l = 0, n
                        do k = 0, m

                            if (.not. wave_propagation) then
                                ! if we use the original HLLC solver, we use the following code
                                rhs_vf(j)%sf(k, l, q) = 1d0/dx(k)* &
                                                        (flux_n(1)%vf(j)%sf(k - 1, l, q) &
                                                        - flux_n(1)%vf(j)%sf(k, l, q))
                                ! print *, "original HLLC rhs:", rhs_vf(j)%sf(k, l, q)
                            else
                                ! if we use the path-conservative HLLC solver, which is based on 
                                ! the wave-propagation method, we use the following code
                                rhs_vf(j)%sf(k, l, q) = -1d0/dx(k)* &
                                                        (inter_fluc_p_n(1)%vf(j)%sf(k-1,l,q) &
                                                        + inter_fluc_m_n(1)%vf(j)%sf(k,l,q)  &
                                                        + cell_fluc_n(1)%vf(j)%sf(k,l,q))
                                ! if (j < 4)  print *, "wave-propagation HLLC rhs:",j, inter_fluc_p_n(1)%vf(j)%sf(k-1,l,q), &
                                !             inter_fluc_m_n(1)%vf(j)%sf(k,l,q), cell_fluc_n(1)%vf(j)%sf(k,l,q)
                            end if

                            if (MOOD_button) then
                                if (.not. wave_propagation) then
                                    lo_rhs_vf(j)%sf(k, l, q) = 1d0/dx(k)* &
                                                            (lo_flux_n(1)%vf(j)%sf(k - 1, l, q) &
                                                            - lo_flux_n(1)%vf(j)%sf(k, l, q))
                                else
                                    lo_rhs_vf(j)%sf(k, l, q) = -1d0/dx(k)* &
                                                            (lo_inter_fluc_p_n(1)%vf(j)%sf(k-1,l,q) &
                                                            + lo_inter_fluc_m_n(1)%vf(j)%sf(k,l,q)  &
                                                            + lo_cell_fluc_n(1)%vf(j)%sf(k,l,q))
                                end if
                            end if

                        end do
                    end do
                end do
            end do

            if (nonconservative_term_energy) then
                ! This part is about the nonconservative term of the internal energy equation in the Saural's six-equation model(2009) alpha*p*dalpha/dx
                ! while it is not used in Pelanti's six-equation model(2014)
                if (model_eqns == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                rho = 0d0
                                do i = 1, num_fluids
                                    rho = rho + q_prim_vf%vf(i + contxb - 1)%sf(j, k, l)
                                end do
                                do i = 1, num_fluids
                                    Y_mf(i) = q_prim_vf%vf(i + contxb - 1)%sf(j, k, l)/max(rho, sgm_eps)
                                end do
                                sum_nonconservative_term = -q_prim_vf%vf(contxe + idir)%sf(j, k, l)* &
                                    (Y_mf(2)*(phasic_pres_src_n(1)%vf(intxb)%sf(j,k,l) - &
                                            phasic_pres_src_n(1)%vf(intxb)%sf(j-1,k,l)) -&
                                     Y_mf(1)*(phasic_pres_src_n(1)%vf(intxe)%sf(j,k,l) - &
                                            phasic_pres_src_n(1)%vf(intxe)%sf(j-1,k,l)))
                                do i = 1, num_fluids
                                    ! ! model using phasic internal energy
                                    ! rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                    !     rhs_vf(i + intxb - 1)%sf(j, k, l) - 1d0/dx(j)* &
                                    !     q_cons_vf%vf(i + advxb - 1)%sf(j, k, l)* &
                                    !     q_prim_vf%vf(E_idx)%sf(j, k, l)* &
                                    !     (flux_src_n(1)%vf(advxb)%sf(j, k, l) - &
                                    !     flux_src_n(1)%vf(advxb)%sf(j - 1, k, l))

                                    ! print *, "sum_nonconservative_term", sum_nonconservative_term
                                    ! model using phasic total energy
                                    rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                        rhs_vf(i + intxb - 1)%sf(j, k, l) + (-1)**i * 1d0/dx(j)* &
                                        sum_nonconservative_term
                                end do

                                if (MOOD_button) then
                                    lo_sum_nonconservative_term = -q_prim_vf%vf(contxe + idir)%sf(j, k, l)* &
                                        (Y_mf(2)*(lo_phasic_pres_src_n(1)%vf(intxb)%sf(j,k,l) - &
                                                lo_phasic_pres_src_n(1)%vf(intxb)%sf(j-1,k,l)) -&
                                         Y_mf(1)*(lo_phasic_pres_src_n(1)%vf(intxe)%sf(j,k,l) - &
                                                lo_phasic_pres_src_n(1)%vf(intxe)%sf(j-1,k,l)))
                                    do i = 1, num_fluids
                                        lo_rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                            lo_rhs_vf(i + intxb - 1)%sf(j, k, l) + (-1)**i * 1d0/dx(j)* &
                                            lo_sum_nonconservative_term
                                    end do
                                end if

                            end do
                        end do
                    end do
                end if
            end if

            ! need to check whether we need to calculate the volume fraction source term or not in PC HLLC solver
            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do j = advxb, advxe
                    do q = 0, p
                        do l = 0, n
                            do k = 0, m
                                rhs_vf(j)%sf(k, l, q) = &
                                    rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                    q_prim_vf%vf(contxe + idir)%sf(k, l, q)* &
                                    (flux_src_n(1)%vf(j)%sf(k - 1, l, q) &
                                    - flux_src_n(1)%vf(j)%sf(k, l, q))
                            end do
                        end do
                    end do
                end do
            else
                if (alt_soundspeed) then
                    do j = advxb, advxe
                        if ((j == advxe) .and. (bubbles .neqv. .true.)) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do q = 0, p
                                do l = 0, n
                                    do k = 0, m
                                        rhs_vf(j)%sf(k, l, q) = &
                                            rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                            (q_cons_vf%vf(j)%sf(k, l, q) - Kterm(k, l, q))* &
                                            (flux_src_n(1)%vf(j)%sf(k, l, q) &
                                            - flux_src_n(1)%vf(j)%sf(k - 1, l, q))
                                    end do
                                end do
                            end do
                        else if ((j == advxb) .and. (bubbles .neqv. .true.)) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do q = 0, p
                                do l = 0, n
                                    do k = 0, m
                                        rhs_vf(j)%sf(k, l, q) = &
                                            rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                            (q_cons_vf%vf(j)%sf(k, l, q) + Kterm(k, l, q))* &
                                            (flux_src_n(1)%vf(j)%sf(k, l, q) &
                                            - flux_src_n(1)%vf(j)%sf(k - 1, l, q))
                                    end do
                                end do
                            end do
                        end if
                    end do
                else
                    if (nonconservative_term_advection) then
                        !$acc parallel loop collapse(4) gang vector default(present)
                        do j = advxb, advxe 
                            do q = 0, p
                                do l = 0, n
                                    do k = 0, m
                                                
                                        ! print *, "rhs before adding source term", rhs_vf(j)%sf(k, l, q)
                                        rhs_vf(j)%sf(k, l, q) = &
                                            rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                            q_cons_vf%vf(j)%sf(k, l, q)* &
                                            (flux_src_n(1)%vf(advxb)%sf(k, l, q) &
                                            - flux_src_n(1)%vf(advxb)%sf(k - 1, l, q))
                                        ! print *, "source term of volume fraction in rhs:", 1d0/dx(k)* &
                                        !     q_cons_vf%vf(j)%sf(k, l, q)* &
                                        !     (flux_src_n(1)%vf(advxb)%sf(k, l, q) &
                                        !     - flux_src_n(1)%vf(advxb)%sf(k - 1, l, q))
                                        ! print *, "rhs after adding source term", rhs_vf(j)%sf(k, l, q)
                                        if (MOOD_button) then
                                            lo_rhs_vf(j)%sf(k, l, q) = &
                                                lo_rhs_vf(j)%sf(k, l, q) + 1d0/dx(k)* &
                                                q_cons_vf%vf(j)%sf(k, l, q)* &
                                                (lo_flux_src_n(1)%vf(advxb)%sf(k, l, q) &
                                                - lo_flux_src_n(1)%vf(advxb)%sf(k - 1, l, q))
                                        end if
                                        
                                    end do
                                end do
                            end do
                        end do
                    end if
                end if
            end if

            ! if (.not. wave_propagation) then
            !     do q = 0, p
            !         do l = 0, n
            !             do k = 0, m
            !                 if (abs(rhs_vf(advxb)%sf(k, l, q) + rhs_vf(advxe)%sf(k, l, q)) > 1d-10) then
            !                     print *, "sum of volume fraction not equal 1 with error", (rhs_vf(advxb)%sf(k, l, q) + rhs_vf(advxe)%sf(k, l, q))
            !                     print *, "sum of volume fraction not equal 1 with error with low-order scheme", &
            !                                 (lo_rhs_vf(advxb)%sf(k, l, q) + lo_rhs_vf(advxe)%sf(k, l, q))
            !                     print *, "sum of rhs in volume fraction equation", 1d0/dx(k)* &
            !                     (q_cons_vf%vf(advxb)%sf(k, l, q) + q_cons_vf%vf(advxe)%sf(k, l, q))* &
            !                     (flux_src_n(1)%vf(advxb)%sf(k, l, q) &
            !                     - flux_src_n(1)%vf(advxb)%sf(k - 1, l, q)) + 1d0/dx(k)* &
            !                     (flux_n(1)%vf(advxb)%sf(k - 1, l, q) - flux_n(1)%vf(advxb)%sf(k, l, q) + &
            !                     flux_n(1)%vf(advxe)%sf(k - 1, l, q) - flux_n(1)%vf(advxe)%sf(k, l, q))
            !                 end if
            !             end do
            !         end do
            !     end do
            ! end if

        elseif (idir == 2) then
            ! RHS Contribution in y-direction ===============================
            ! Applying the Riemann fluxes

            if (bc_y%beg <= -5 .and. bc_y%beg >= -13) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, &
                           flux_src_n(idir)%vf, idir, -1, ix, iy, iz)
                if (MOOD_button) then
                    call s_cbc(q_prim_vf%vf, lo_flux_n(idir)%vf, &
                               lo_flux_src_n(idir)%vf, idir, -1, ix, iy, iz)
                end if
            end if

            if (bc_y%end <= -5 .and. bc_y%end >= -13) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, &
                           flux_src_n(idir)%vf, idir, 1, ix, iy, iz)
                if (MOOD_button) then
                    call s_cbc(q_prim_vf%vf, lo_flux_n(idir)%vf, &
                               lo_flux_src_n(idir)%vf, idir, 1, ix, iy, iz)
                end if
            end if

            !$acc parallel loop collapse(4) gang vector default(present)
            do j = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do q = 0, m
                            if (.not. wave_propagation) then 
                                ! if original HLLC solver is used, we use the following code
                                rhs_vf(j)%sf(q, k, l) = &
                                    rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                    (flux_n(2)%vf(j)%sf(q, k - 1, l) &
                                    - flux_n(2)%vf(j)%sf(q, k, l))
                            else
                                ! if path-conservative HLLC solver is used, we use the following code
                                rhs_vf(j)%sf(q, k, l) = &
                                    rhs_vf(j)%sf(q, k, l) - 1d0/dy(k)* &
                                    (inter_fluc_p_n(2)%vf(j)%sf(q,k-1,l) &
                                    + inter_fluc_m_n(2)%vf(j)%sf(q,k,l)  &
                                    + cell_fluc_n(2)%vf(j)%sf(q,k,l))
                            end if

                            if (MOOD_button) then
                                if (.not. wave_propagation) then
                                    lo_rhs_vf(j)%sf(q, k, l) = &
                                        lo_rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                        (lo_flux_n(2)%vf(j)%sf(q, k - 1, l) &
                                        - lo_flux_n(2)%vf(j)%sf(q, k, l))
                                else
                                    lo_rhs_vf(j)%sf(q, k, l) = &
                                        lo_rhs_vf(j)%sf(q, k, l) - 1d0/dy(k)* &
                                        (lo_inter_fluc_p_n(2)%vf(j)%sf(q,k-1,l) &
                                        + lo_inter_fluc_m_n(2)%vf(j)%sf(q,k,l)  &
                                        + lo_cell_fluc_n(2)%vf(j)%sf(q,k,l))
                                end if
                            end if
                        end do
                    end do
                end do
            end do

            if (nonconservative_term_energy) then
                ! ! Not used in Pelanti's six-equation model(2014)
                if (model_eqns == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                rho = 0d0
                                do i = 1, num_fluids
                                    rho = rho + q_prim_vf%vf(i + contxb - 1)%sf(j, k, l)
                                end do
                                do i = 1, num_fluids
                                    Y_mf(i) = q_prim_vf%vf(i + contxb - 1)%sf(j, k, l)/max(rho, sgm_eps)
                                end do
                                sum_nonconservative_term = -q_prim_vf%vf(contxe + idir)%sf(j, k, l)* &
                                    (Y_mf(2)*(phasic_pres_src_n(2)%vf(intxb)%sf(j,k,l) - &
                                            phasic_pres_src_n(2)%vf(intxb)%sf(j,k-1,l)) -&
                                     Y_mf(1)*(phasic_pres_src_n(2)%vf(intxe)%sf(j,k,l) - &
                                            phasic_pres_src_n(2)%vf(intxe)%sf(j,k-1,l)))
                                do i = 1, num_fluids
                                    ! ! model using phasic internal energy
                                    ! rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                    !     rhs_vf(i + intxb - 1)%sf(j, k, l) - 1d0/dy(k)* &
                                    !     q_cons_vf%vf(i + advxb - 1)%sf(j, k, l)* &
                                    !     q_prim_vf%vf(E_idx)%sf(j, k, l)* &
                                    !     (flux_src_n(2)%vf(advxb)%sf(j, k, l) - &
                                    !     flux_src_n(2)%vf(advxb)%sf(j, k - 1, l))
                                    ! model using phasic total energy
                                    rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                        rhs_vf(i + intxb - 1)%sf(j, k, l) + (-1)**i * 1d0/dy(k)* &
                                        sum_nonconservative_term
                                end do

                                if (MOOD_button) then
                                    lo_sum_nonconservative_term = -q_prim_vf%vf(contxe + idir)%sf(j, k, l)* &
                                        (Y_mf(2)*(lo_phasic_pres_src_n(2)%vf(intxb)%sf(j,k,l) - &
                                                lo_phasic_pres_src_n(2)%vf(intxb)%sf(j,k-1,l)) -&
                                         Y_mf(1)*(lo_phasic_pres_src_n(2)%vf(intxe)%sf(j,k,l) - &
                                                lo_phasic_pres_src_n(2)%vf(intxe)%sf(j,k-1,l)))
                                    do i = 1, num_fluids
                                        lo_rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                            lo_rhs_vf(i + intxb - 1)%sf(j, k, l) + (-1)**i * 1d0/dy(k)* &
                                            lo_sum_nonconservative_term
                                    end do
                                end if
                            end do
                        end do
                    end do

                    if (cyl_coord) then
                        !$acc parallel loop collapse(4) gang vector default(present)
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    do i = 1, num_fluids
                                        rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                            rhs_vf(i + intxb - 1)%sf(j, k, l) - 5d-1/y_cc(k)* &
                                            q_cons_vf%vf(i + advxb - 1)%sf(j, k, l)* &
                                            q_prim_vf%vf(E_idx)%sf(j, k, l)* &
                                            (flux_src_n(2)%vf(advxb)%sf(j, k, l) + &
                                            flux_src_n(2)%vf(advxb)%sf(j, k - 1, l))
                                    end do
                                end do
                            end do
                        end do
                    end if
                end if

            end if

            if (cyl_coord) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do q = 0, m
                                rhs_vf(j)%sf(q, k, l) = &
                                    rhs_vf(j)%sf(q, k, l) - 5d-1/y_cc(k)* &
                                    (flux_gsrc_n(2)%vf(j)%sf(q, k - 1, l) &
                                    + flux_gsrc_n(2)%vf(j)%sf(q, k, l))
                            end do
                        end do
                    end do
                end do
            end if

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do j = advxb, advxe
                    do l = 0, p
                        do k = 0, n
                            do q = 0, m
                                rhs_vf(j)%sf(q, k, l) = &
                                    rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                    q_prim_vf%vf(contxe + idir)%sf(q, k, l)* &
                                    (flux_src_n(2)%vf(j)%sf(q, k - 1, l) &
                                    - flux_src_n(2)%vf(j)%sf(q, k, l))
                            end do
                        end do
                    end do
                end do
            else

                if (alt_soundspeed) then
                    do j = advxb, advxe
                        if ((j == advxe) .and. (bubbles .neqv. .true.)) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do l = 0, p
                                do k = 0, n
                                    do q = 0, m
                                        rhs_vf(j)%sf(q, k, l) = &
                                            rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                            (q_cons_vf%vf(j)%sf(q, k, l) - Kterm(q, k, l))* &
                                            (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                            - flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                    end do
                                end do
                            end do
                            if (cyl_coord) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do l = 0, p
                                    do k = 0, n
                                        do q = 0, m
                                            rhs_vf(j)%sf(q, k, l) = &
                                                rhs_vf(j)%sf(q, k, l) - &
                                                (Kterm(q, k, l)/2d0/y_cc(k))* &
                                                (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                + flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                        end do
                                    end do
                                end do
                            end if
                        else if ((j == advxb) .and. (bubbles .neqv. .true.)) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do l = 0, p
                                do k = 0, n
                                    do q = 0, m
                                        rhs_vf(j)%sf(q, k, l) = &
                                            rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                            (q_cons_vf%vf(j)%sf(q, k, l) + Kterm(q, k, l))* &
                                            (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                            - flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                    end do
                                end do
                            end do
                            if (cyl_coord) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do l = 0, p
                                    do k = 0, n
                                        do q = 0, m
                                            rhs_vf(j)%sf(q, k, l) = &
                                                rhs_vf(j)%sf(q, k, l) + &
                                                (Kterm(q, k, l)/2d0/y_cc(k))* &
                                                (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                + flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                        end do
                                    end do
                                end do
                            end if
                        end if
                    end do
                else
                    if (nonconservative_term_advection) then
                        !$acc parallel loop collapse(4) gang vector default(present)
                        do j = advxb, advxe
                            do l = 0, p
                                do k = 0, n
                                    do q = 0, m
                                        rhs_vf(j)%sf(q, k, l) = &
                                            rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                            q_cons_vf%vf(j)%sf(q, k, l)* &
                                            (flux_src_n(2)%vf(advxb)%sf(q, k, l) &
                                            - flux_src_n(2)%vf(advxb)%sf(q, k - 1, l))

                                        if (MOOD_button) then
                                            lo_rhs_vf(j)%sf(q, k, l) = &
                                                lo_rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                                q_cons_vf%vf(j)%sf(q, k, l)* &
                                                (lo_flux_src_n(2)%vf(advxb)%sf(q, k, l) &
                                                - lo_flux_src_n(2)%vf(advxb)%sf(q, k - 1, l))
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end if
                end if
            end if

        elseif (idir == 3) then
            ! RHS Contribution in z-direction ===============================

            ! Applying the Riemann fluxes

            if (bc_z%beg <= -5 .and. bc_z%beg >= -13) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, &
                           flux_src_n(idir)%vf, idir, -1, ix, iy, iz)
                if (MOOD_button) then
                    call s_cbc(q_prim_vf%vf, lo_flux_n(idir)%vf, &
                               lo_flux_src_n(idir)%vf, idir, -1, ix, iy, iz)
                end if
            end if

            if (bc_z%end <= -5 .and. bc_z%end >= -13) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, &
                           flux_src_n(idir)%vf, idir, 1, ix, iy, iz)
                if (MOOD_button) then
                    call s_cbc(q_prim_vf%vf, lo_flux_n(idir)%vf, &
                               lo_flux_src_n(idir)%vf, idir, 1, ix, iy, iz)
                end if
            end if

            if (grid_geometry == 3) then ! Cylindrical Coordinates
                !$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, sys_size
                    do k = 0, p
                        do q = 0, n
                            do l = 0, m
                                rhs_vf(j)%sf(l, q, k) = &
                                    rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)/y_cc(q)* &
                                    q_prim_vf%vf(contxe + idir)%sf(l, q, k)* &
                                    (flux_n(3)%vf(j)%sf(l, q, k - 1) &
                                     - flux_n(3)%vf(j)%sf(l, q, k))
                            end do
                        end do
                    end do
                end do

                !$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, sys_size
                    do k = 0, p
                        do q = 0, n
                            do l = 0, m
                                rhs_vf(j)%sf(l, q, k) = &
                                    rhs_vf(j)%sf(l, q, k) - 5d-1/y_cc(q)* &
                                    (flux_gsrc_n(3)%vf(j)%sf(l, q, k - 1) &
                                     - flux_gsrc_n(3)%vf(j)%sf(l, q, k))
                            end do
                        end do
                    end do
                end do

            else ! Cartesian Coordinates
                !$acc parallel loop collapse(4) gang vector default(present)
                do j = 1, sys_size
                    do k = 0, p
                        do q = 0, n
                            do l = 0, m
                                if (.not. wave_propagation) then
                                    ! if original HLLC solver is used, we use the following code
                                    rhs_vf(j)%sf(l, q, k) = &
                                        rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                        (flux_n(3)%vf(j)%sf(l, q, k - 1) &
                                        - flux_n(3)%vf(j)%sf(l, q, k))
                                else
                                    ! if path-conservative HLLC solver is used, we use the following code
                                    rhs_vf(j)%sf(l, q, k) = &
                                        rhs_vf(j)%sf(l, q, k) - 1d0/dz(k)* &
                                        (inter_fluc_p_n(3)%vf(j)%sf(l,q,k-1) &
                                        + inter_fluc_m_n(3)%vf(j)%sf(l,q,k)  &
                                        + cell_fluc_n(3)%vf(j)%sf(l,q,k))
                                end if

                                if (MOOD_button) then
                                    if (.not. wave_propagation) then
                                        lo_rhs_vf(j)%sf(l, q, k) = &
                                            lo_rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                            (lo_flux_n(3)%vf(j)%sf(l, q, k - 1) &
                                            - lo_flux_n(3)%vf(j)%sf(l, q, k))
                                    else
                                        lo_rhs_vf(j)%sf(l, q, k) = &
                                            lo_rhs_vf(j)%sf(l, q, k) - 1d0/dz(k)* &
                                            (lo_inter_fluc_p_n(3)%vf(j)%sf(l,q,k-1) &
                                            + lo_inter_fluc_m_n(3)%vf(j)%sf(l,q,k)  &
                                            + lo_cell_fluc_n(3)%vf(j)%sf(l,q,k))
                                    end if
                                end if
                            end do
                        end do
                    end do
                end do
            end if

            if (nonconservative_term_energy) then
                ! Not used in Pelanti's six-equation model(2014)
                if (model_eqns == 3) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                rho = 0d0
                                do i = 1, num_fluids
                                    rho = rho + q_prim_vf%vf(i + contxb - 1)%sf(j, k, l)
                                end do
                                do i = 1, num_fluids
                                    Y_mf(i) = q_prim_vf%vf(i + contxb - 1)%sf(j, k, l)/max(rho, sgm_eps)
                                end do
                                sum_nonconservative_term = -q_prim_vf%vf(contxe + idir)%sf(j, k, l)* &
                                    (Y_mf(2)*(phasic_pres_src_n(3)%vf(intxb)%sf(j,k,l) - &
                                            phasic_pres_src_n(3)%vf(intxb)%sf(j,k,l-1)) -&
                                     Y_mf(1)*(phasic_pres_src_n(3)%vf(intxe)%sf(j,k,l) - &
                                            phasic_pres_src_n(3)%vf(intxe)%sf(j,k,l-1)))
                                do i = 1, num_fluids
                                    ! ! model using phasic internal energy
                                    ! rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                    !     rhs_vf(i + intxb - 1)%sf(j, k, l) - 1d0/dz(l)* &
                                    !     q_cons_vf%vf(i + advxb - 1)%sf(j, k, l)* &
                                    !     q_prim_vf%vf(E_idx)%sf(j, k, l)* &
                                    !     (flux_src_n(3)%vf(advxb)%sf(j, k, l) - &
                                    !     flux_src_n(3)%vf(advxb)%sf(j, k, l - 1))
                                    ! model using phasic total energy
                                    rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                        rhs_vf(i + intxb - 1)%sf(j, k, l) + (-1)**i * 1d0/dz(l)* &
                                        sum_nonconservative_term
                                end do
                                if (MOOD_button) then
                                    lo_sum_nonconservative_term = -q_prim_vf%vf(contxe + idir)%sf(j, k, l)* &
                                        (Y_mf(2)*(lo_phasic_pres_src_n(3)%vf(intxb)%sf(j,k,l) - &
                                                lo_phasic_pres_src_n(3)%vf(intxb)%sf(j,k,l-1)) -&
                                         Y_mf(1)*(lo_phasic_pres_src_n(3)%vf(intxe)%sf(j,k,l) - &
                                                lo_phasic_pres_src_n(3)%vf(intxe)%sf(j,k,l-1)))
                                    do i = 1, num_fluids
                                        lo_rhs_vf(i + intxb - 1)%sf(j, k, l) = &
                                            lo_rhs_vf(i + intxb - 1)%sf(j, k, l) + (-1)**i * 1d0/dz(l)* &
                                            lo_sum_nonconservative_term
                                    end do
                                end if
                            end do
                        end do
                    end do
                end if

            end if

            if (grid_geometry == 3) then
                if (riemann_solver == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do j = advxb, advxe
                        do l = 0, p
                            do k = 0, n
                                do q = 0, m
                                    rhs_vf(j)%sf(q, k, l) = &
                                        rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                        q_prim_vf%vf(contxe + idir)%sf(q, k, l)* &
                                        (flux_src_n(2)%vf(j)%sf(q, k - 1, l) &
                                        - flux_src_n(2)%vf(j)%sf(q, k, l))
                                end do
                            end do
                        end do
                    end do
                else

                    if (alt_soundspeed) then
                        do j = advxb, advxe
                            if ((j == advxe) .and. (bubbles .neqv. .true.)) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do l = 0, p
                                    do k = 0, n
                                        do q = 0, m
                                            rhs_vf(j)%sf(q, k, l) = &
                                                rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                                (q_cons_vf%vf(j)%sf(q, k, l) - Kterm(q, k, l))* &
                                                (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                - flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                        end do
                                    end do
                                end do
                                if (cyl_coord) then
                                    !$acc parallel loop collapse(3) gang vector default(present)
                                    do l = 0, p
                                        do k = 0, n
                                            do q = 0, m
                                                rhs_vf(j)%sf(q, k, l) = &
                                                    rhs_vf(j)%sf(q, k, l) - &
                                                    (Kterm(q, k, l)/2d0/y_cc(k))* &
                                                    (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                    + flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                            end do
                                        end do
                                    end do
                                end if
                            else if ((j == advxb) .and. (bubbles .neqv. .true.)) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do l = 0, p
                                    do k = 0, n
                                        do q = 0, m
                                            rhs_vf(j)%sf(q, k, l) = &
                                                rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                                (q_cons_vf%vf(j)%sf(q, k, l) + Kterm(q, k, l))* &
                                                (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                - flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                        end do
                                    end do
                                end do
                                if (cyl_coord) then
                                    !$acc parallel loop collapse(3) gang vector default(present)
                                    do l = 0, p
                                        do k = 0, n
                                            do q = 0, m
                                                rhs_vf(j)%sf(q, k, l) = &
                                                    rhs_vf(j)%sf(q, k, l) + &
                                                    (Kterm(q, k, l)/2d0/y_cc(k))* &
                                                    (flux_src_n(2)%vf(j)%sf(q, k, l) &
                                                    + flux_src_n(2)%vf(j)%sf(q, k - 1, l))
                                            end do
                                        end do
                                    end do
                                end if
                            end if
                        end do
                    else
                        if (nonconservative_term_advection) then
                            !$acc parallel loop collapse(4) gang vector default(present)
                            do j = advxb, advxe
                                do l = 0, p
                                    do k = 0, n
                                        do q = 0, m
                                            rhs_vf(j)%sf(q, k, l) = &
                                                rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                                q_cons_vf%vf(j)%sf(q, k, l)* &
                                                (flux_src_n(2)%vf(advxb)%sf(q, k, l) &
                                                - flux_src_n(2)%vf(advxb)%sf(q, k - 1, l))
                                            if (MOOD_button) then
                                                lo_rhs_vf(j)%sf(q, k, l) = &
                                                    lo_rhs_vf(j)%sf(q, k, l) + 1d0/dy(k)* &
                                                    q_cons_vf%vf(j)%sf(q, k, l)* &
                                                    (lo_flux_src_n(2)%vf(advxb)%sf(q, k, l) &
                                                    - lo_flux_src_n(2)%vf(advxb)%sf(q, k - 1, l))
                                            end if
                                        end do
                                    end do
                                end do
                            end do
                        end if
                    end if
                end if
            else
                if (riemann_solver == 1) then
                    !$acc parallel loop collapse(4) gang vector default(present)
                    do j = advxb, advxe
                        do k = 0, p
                            do q = 0, n
                                do l = 0, m
                                    rhs_vf(j)%sf(l, q, k) = &
                                        rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                        q_prim_vf%vf(contxe + idir)%sf(l, q, k)* &
                                        (flux_src_n(3)%vf(j)%sf(l, q, k - 1) &
                                        - flux_src_n(3)%vf(j)%sf(l, q, k))
                                end do
                            end do
                        end do
                    end do
                else
                    if (alt_soundspeed) then
                        do j = advxb, advxe
                            if ((j == advxe) .and. (bubbles .neqv. .true.)) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do k = 0, p
                                    do q = 0, n
                                        do l = 0, m
                                            rhs_vf(j)%sf(l, q, k) = &
                                                rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                                (q_cons_vf%vf(j)%sf(l, q, k) - Kterm(l, q, k))* &
                                                (flux_src_n(3)%vf(j)%sf(l, q, k) &
                                                - flux_src_n(3)%vf(j)%sf(l, q, k - 1))
                                        end do
                                    end do
                                end do
                            else if ((j == advxb) .and. (bubbles .neqv. .true.)) then
                                !$acc parallel loop collapse(3) gang vector default(present)
                                do k = 0, p
                                    do q = 0, n
                                        do l = 0, m
                                            rhs_vf(j)%sf(l, q, k) = &
                                                rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                                (q_cons_vf%vf(j)%sf(l, q, k) + Kterm(l, q, k))* &
                                                (flux_src_n(3)%vf(j)%sf(l, q, k) &
                                                - flux_src_n(3)%vf(j)%sf(l, q, k - 1))
                                        end do
                                    end do
                                end do
                            end if
                        end do
                    else
                        if (nonconservative_term_advection) then
                            !$acc parallel loop collapse(4) gang vector default(present)
                            do j = advxb, advxe
                                do k = 0, p
                                    do q = 0, n
                                        do l = 0, m
                                            rhs_vf(j)%sf(l, q, k) = &
                                                rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                                q_cons_vf%vf(j)%sf(l, q, k)* &
                                                (flux_src_n(3)%vf(advxb)%sf(l, q, k) &
                                                - flux_src_n(3)%vf(advxb)%sf(l, q, k - 1))
                                            if (MOOD_button) then
                                                lo_rhs_vf(j)%sf(l, q, k) = &
                                                    lo_rhs_vf(j)%sf(l, q, k) + 1d0/dz(k)* &
                                                    q_cons_vf%vf(j)%sf(l, q, k)* &
                                                    (lo_flux_src_n(3)%vf(advxb)%sf(l, q, k) &
                                                    - lo_flux_src_n(3)%vf(advxb)%sf(l, q, k - 1))
                                            end if
                                        end do
                                    end do
                                end do
                            end do
                        end if
                    end if
                end if
            end if
        end if ! id loop

    end subroutine s_compute_advection_source_term

    !>  The purpose of this procedure is to relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure with anyrate, from which the corresponding volume
        !!      fraction of each phase are recomputed. 
        !!      Method corresponding to paper `Marica Pelanti 2022 IJMF 6eq phase change arbitrary-rate`
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_pressure_relaxation_Pelanti2022(q_cons_vf) ! ----------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf

        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        real(kind(0d0)) :: pres_relax
        real(kind(0d0)), dimension(num_fluids) :: pres_K_init
        real(kind(0d0)) :: f_pres
        real(kind(0d0)) :: df_pres
        real(kind(0d0)), dimension(num_fluids) :: rho_K_s
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho
        real(kind(0d0)), dimension(num_fluids) :: alpha, rho_e, rho_K
        real(kind(0d0)), dimension(num_fluids) :: xi_p, xi_p0, Gamma_pelanti, chi_pelanti
        real(kind(0d0)) :: sum_alpha
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: mQ
        real(kind(0d0)), dimension(2) :: Re
        real(kind(0d0)) :: a0, b0, d0, delta0
        real(kind(0d0)) :: p_I, deleta_p0

        integer :: i, j, k, l, q, iter !< Generic loop iterators
        integer :: relax !< Relaxation procedure determination variable
        logical :: real_in_hyperbolic_step
        real(kind(0d0)), dimension(1:sys_size) :: initial_q_cons_vf

        !$acc parallel loop collapse(3) gang vector private(pres_K_init, rho_K_s, alpha_rho, alpha, Re, pres_relax)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Check whether the solution from the hyperbolic step is real
                    real_in_hyperbolic_step = .true.
                    do i = 1, sys_size
                        initial_q_cons_vf(i) = q_cons_vf(i)%sf(j, k, l)
                        if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                            real_in_hyperbolic_step = .false.
                            call s_mpi_abort('NAD criteria violated in hyperbolic step')
                        end if
                    end do

                    ! Numerical correction of the volume fractions
                    if (mpp_lim) then
                        sum_alpha = 0d0

                        !$acc loop seq
                        do i = 1, num_fluids
                            if ((q_cons_vf(i + contxb - 1)%sf(j, k, l) <= 0d0) .or. &
                                (q_cons_vf(i + advxb - 1)%sf(j, k, l) <= 0d0)) then
                                q_cons_vf(i + contxb - 1)%sf(j, k, l) = 1d-8
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 1d-8
                                q_cons_vf(i + intxb - 1)%sf(j, k, l) = 1d-8
                            end if

                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) >= 1d0) &
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 1d0 - 1d-8
                            sum_alpha = sum_alpha + q_cons_vf(i + advxb - 1)%sf(j, k, l)
                        end do

                        !$acc loop seq
                        do i = 1, num_fluids
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)/sum_alpha
                        end do
                    end if

                    ! Pressures relaxation procedure ===================================
                    rho = 0d0
                    !$acc loop seq
                    do i = 1, num_fluids
                        rho = rho + q_cons_vf(i)%sf(j, k, l)
                    end do

                    dyn_pres = 0d0

                    !$acc loop seq
                    do i = momxb, momxe
                        dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j, k, l)* &
                                q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                    end do
                    ! Initial partial pressures
                    !$acc loop seq
                    do i = 1, num_fluids
                        pres_K_init(i) = &
                            ((q_cons_vf(i + intxb - 1)%sf(j, k, l) - &
                                q_cons_vf(i)%sf(j, k, l)*qvs(i) - &
                                dyn_pres/max(rho, sgm_eps)*q_cons_vf(i)%sf(j, k, l))/ &
                                max(q_cons_vf(i + advxb - 1)%sf(j, k, l), sgm_eps) &
                                - pi_infs(i))/gammas(i)
                    end do

                    !$acc loop seq
                    do i = 1, num_fluids
                        alpha(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)
                        alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                    end do

                    p_I = alpha(1)*pres_K_init(1) + alpha(2)*pres_K_init(2)

                    ! print *, p_I

                    !$acc loop seq
                    do i = 1, num_fluids
                        rho_K(i) = alpha_rho(i)/alpha(i)
                        Gamma_pelanti(i) = 1d0/gammas(i)
                        chi_pelanti(i) = - qvs(i)/gammas(i)
                        rho_e(i) = (q_cons_vf(i + intxb - 1)%sf(j, k, l) - &
                                    dyn_pres/max(rho, sgm_eps)*alpha_rho(i))/alpha(i)
                        xi_p0(i) = (Gamma_pelanti(i)*(rho_e(i) + p_I) + chi_pelanti(i)*rho_K(i))/alpha(i)
                        xi_p(i) = xi_p0(i)
                    end do
                    alpha(1) = alpha(1) - (pres_K_init(2) - pres_K_init(1))/(xi_p0(1) + xi_p0(2))
                    ! print *, p_I, alpha(1)
                    alpha(2) = 1d0 - alpha(1)

                    ! Pressure calculation ==================================
                    gamma = 0d0
                    pi_inf = 0d0
                    mQ = 0d0

                    sum_alpha = 0d0

                    if (mpp_lim) then
                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho(i) = max(1d-8, alpha_rho(i))
                            alpha(i) = min(max(1d-8, alpha(i)), 1d0 - 1d-8)
                            sum_alpha = sum_alpha + alpha(i)
                        end do

                        alpha = alpha/max(sum_alpha, sgm_eps)

                    end if

                    !$acc loop seq
                    do i = 1, num_fluids
                        gamma = gamma + alpha(i)*gammas(i)
                        pi_inf = pi_inf + alpha(i)*pi_infs(i)
                        mQ = mQ + alpha_rho(i)*qvs(i)
                    end do

                    pres_relax = (q_cons_vf(E_idx)%sf(j, k, l) - dyn_pres - mQ &
                                 - pi_inf)/gamma

                    !$acc loop seq
                    do i = 1, num_fluids
                        q_cons_vf(i + advxb - 1)%sf(j, k, l) = alpha(i)
                        q_cons_vf(i + intxb - 1)%sf(j, k, l) = &
                            alpha(i)*(gammas(i)*pres_relax + pi_infs(i)) + &
                            q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i) + &
                            dyn_pres/max(rho, sgm_eps)*q_cons_vf(i)%sf(j, k, l)
                    end do
                    ! ==================================================================
                    do i = 1, sys_size
                        if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l)) .and. real_in_hyperbolic_step) then
                            print *, 'troubled cell caused by NAD in relaxation step', i
                            ! print *, 'relax =', relax
                            ! print *, 'a0', a0
                            ! print *, 'b0', b0
                            ! print *, 'd0', d0
                            print *, 'pres_relax', pres_relax
                            print *, 'pres_K_init', pres_K_init
                            ! print *, 'alpha0', q_cons_vf(advxb)%sf(j, k, l), q_cons_vf(advxe)%sf(j, k, l)
                            do q = 1, sys_size
                                print *, 'initial_q_cons_vf', initial_q_cons_vf(q)
                                print *, 'q_cons_vf', q_cons_vf(q)%sf(j, k, l)
                            end do
                            call s_mpi_abort('NAD criteria violated in pressure relaxation step')
                        end if
                    end do
                end do
            end do
        end do

    end subroutine s_pressure_relaxation_Pelanti2022 ! -----------------------


    !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal total energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed.
        !!      Method corresponding to paper `Hiro Wakimura 2024 JCP 6eq phase change BVD`
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_pressure_relaxation_Wakimura2024(q_cons_vf) ! ----------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf

        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        real(kind(0d0)) :: pres_relax
        real(kind(0d0)), dimension(num_fluids) :: pres_K_init
        real(kind(0d0)) :: f_pres
        real(kind(0d0)) :: df_pres
        real(kind(0d0)), dimension(num_fluids) :: rho_K_s
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho
        real(kind(0d0)), dimension(num_fluids) :: alpha
        real(kind(0d0)) :: sum_alpha
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: mQ
        real(kind(0d0)), dimension(2) :: Re
        real(kind(0d0)) :: a0, b0, d0, delta0

        integer :: i, j, k, l, q, iter !< Generic loop iterators
        integer :: relax !< Relaxation procedure determination variable
        logical :: real_in_hyperbolic_step
        real(kind(0d0)), dimension(1:sys_size) :: initial_q_cons_vf

        !$acc parallel loop collapse(3) gang vector private(pres_K_init, rho_K_s, alpha_rho, alpha, Re, pres_relax)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Check whether the solution from the hyperbolic step is real
                    real_in_hyperbolic_step = .true.
                    do i = 1, sys_size
                        initial_q_cons_vf(i) = q_cons_vf(i)%sf(j, k, l)
                        if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                            real_in_hyperbolic_step = .false.
                            print *, 'troubled cell caused by NAD in hyperbolic step', i
                            do q = 1, sys_size
                                print *, 'initial_q_cons_vf', initial_q_cons_vf(q)
                            end do
                            call s_mpi_abort('NAD criteria violated')
                        end if
                    end do

                    ! Numerical correction of the volume fractions
                    if (mpp_lim) then
                        sum_alpha = 0d0

                        !$acc loop seq
                        do i = 1, num_fluids
                            if ((q_cons_vf(i + contxb - 1)%sf(j, k, l) < 0d0) .or. &
                                (q_cons_vf(i + advxb - 1)%sf(j, k, l) < 0d0)) then
                                q_cons_vf(i + contxb - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + intxb - 1)%sf(j, k, l) = 0d0
                            end if

                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > 1d0) &
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 1d0
                            sum_alpha = sum_alpha + q_cons_vf(i + advxb - 1)%sf(j, k, l)
                        end do

                        !$acc loop seq
                        do i = 1, num_fluids
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)/sum_alpha
                        end do
                    end if

                    ! Pressures relaxation procedure ===================================

                    ! Is the pressure relaxation procedure necessary?
                    relax = 1

                    !$acc loop seq
                    do i = 1, num_fluids
                        if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > (1d0 - sgm_eps)) relax = 0
                    end do

                    if (relax == 1) then
                        rho = 0d0
                        !$acc loop seq
                        do i = 1, num_fluids
                            rho = rho + q_cons_vf(i)%sf(j, k, l)
                        end do

                        dyn_pres = 0d0

                        !$acc loop seq
                        do i = momxb, momxe
                            dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j, k, l)* &
                                    q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                        end do
                        ! Initial partial pressures
                        !$acc loop seq
                        do i = 1, num_fluids
                            pres_K_init(i) = &
                                ((q_cons_vf(i + intxb - 1)%sf(j, k, l) - &
                                    q_cons_vf(i)%sf(j, k, l)*qvs(i) - &
                                    dyn_pres/max(rho, sgm_eps)*q_cons_vf(i)%sf(j, k, l))/ &
                                    max(q_cons_vf(i + advxb - 1)%sf(j, k, l), sgm_eps) &
                                    - pi_infs(i))/gammas(i)
                        end do
                        ! coefficients of the quadratic equation
                        a0 = gs_min(1)*q_cons_vf(advxe)%sf(j, k, l) + &
                                gs_min(2)*q_cons_vf(advxb)%sf(j, k, l)
                        b0 = -gs_min(1)*q_cons_vf(advxe)%sf(j, k, l)*(pres_K_init(2) - ps_inf(1)) &
                                -gs_min(2)*q_cons_vf(advxb)%sf(j, k, l)*(pres_K_init(1) - ps_inf(2))
                        d0 = -gs_min(1)*q_cons_vf(advxe)%sf(j, k, l)*pres_K_init(2)*ps_inf(1) &
                                -gs_min(2)*q_cons_vf(advxb)%sf(j, k, l)*pres_K_init(1)*ps_inf(2)
                        delta0 = b0**2 - 4d0*a0*d0
                        ! solution of the quadratic equation
                        if (b0 > 0) then
                            pres_relax = 2d0*d0/(-b0 - dsqrt(delta0))
                        else
                            pres_relax = (-b0 + dsqrt(delta0))/(2d0*a0)
                        end if

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)* &
                                    ((gs_min(i) - 1d0)*pres_relax + pres_K_init(i) + gs_min(i)*ps_inf(i))/ &
                                    ((gs_min(i) - 1d0)*pres_relax + pres_relax + gs_min(i)*ps_inf(i))
                            alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                        end do

                        if (mpp_lim) then
                            sum_alpha = 0d0
                            !$acc loop seq
                            do i = 1, num_fluids
                                alpha_rho(i) = max(0d0, alpha_rho(i))
                                alpha(i) = min(max(0d0, alpha(i)), 1d0)
                                sum_alpha = sum_alpha + alpha(i)
                            end do

                            alpha = alpha/max(sum_alpha, sgm_eps)

                        end if

                        !$acc loop seq
                        do i = 1, num_fluids
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = alpha(i)
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = &
                                    alpha(i)*(gammas(i)*pres_relax + pi_infs(i)) + &
                                    alpha_rho(i)*qvs(i) + &
                                    dyn_pres*alpha_rho(i)/max(rho, sgm_eps)
                        end do
                    end if

                    do i = 1, sys_size
                        if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l)) .and. real_in_hyperbolic_step) then
                            print *, 'troubled cell caused by NAD in relaxation step', i
                            ! print *, 'relax =', relax
                            ! print *, 'a0', a0
                            ! print *, 'b0', b0
                            ! print *, 'd0', d0
                            print *, 'delta0', delta0
                            print *, 'pres_relax', pres_relax
                            print *, 'pres_K_init', pres_K_init
                            ! do q = 1, sys_size
                            !     print *, 'initial_q_cons_vf', initial_q_cons_vf(q)
                            !     print *, 'q_cons_vf', q_cons_vf(q)%sf(j, k, l)
                            ! end do
                            call s_mpi_abort('NAD criteria violated')
                        end if
                    end do
                    ! ==================================================================
                end do
            end do
        end do

    end subroutine s_pressure_relaxation_Wakimura2024 ! -----------------------

    !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal total energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed.
        !!      Method corresponding to paper `Marica Pelanti 2014 JCP mixture-energy-consistent 6eq`
        !!  @note This method is different from the one in Wakimura 2024 in terms of the definition of interfacial pressure
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_pressure_relaxation_Pelanti2014(q_cons_vf) ! ----------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf

        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        real(kind(0d0)) :: pres_relax
        real(kind(0d0)), dimension(num_fluids) :: pres_K_init
        real(kind(0d0)) :: f_pres
        real(kind(0d0)) :: df_pres
        real(kind(0d0)), dimension(num_fluids) :: rho_K_s
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho
        real(kind(0d0)), dimension(num_fluids) :: alpha, alpha0, rhok0, Zk0, ck0
        real(kind(0d0)) :: sum_alpha
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: mQ
        real(kind(0d0)), dimension(2) :: Re
        real(kind(0d0)) :: a0, b0, d0, delta0
        real(kind(0d0)) :: C1, C2
        real(kind(0d0)) :: p_I0, p_I_avg, p_I

        integer :: i, j, k, l, q, iter !< Generic loop iterators
        integer :: relax !< Relaxation procedure determination variable
        logical :: real_in_hyperbolic_step
        real(kind(0d0)), dimension(1:sys_size) :: initial_q_cons_vf

        !$acc parallel loop collapse(3) gang vector private(pres_K_init, rho_K_s, alpha_rho, alpha, Re, pres_relax)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Check whether the solution from the hyperbolic step is real
                    real_in_hyperbolic_step = .true.
                    do i = 1, sys_size
                        initial_q_cons_vf(i) = q_cons_vf(i)%sf(j, k, l)
                        if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l))) then
                            real_in_hyperbolic_step = .false.
                        end if
                    end do

                    ! Numerical correction of the volume fractions
                    if (mpp_lim) then
                        sum_alpha = 0d0

                        !$acc loop seq
                        do i = 1, num_fluids
                            if ((q_cons_vf(i + contxb - 1)%sf(j, k, l) < 0d0) .or. &
                                (q_cons_vf(i + advxb - 1)%sf(j, k, l) < 0d0)) then
                                q_cons_vf(i + contxb - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + intxb - 1)%sf(j, k, l) = 0d0
                            end if

                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > 1d0) &
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 1d0
                            sum_alpha = sum_alpha + q_cons_vf(i + advxb - 1)%sf(j, k, l)
                        end do

                        !$acc loop seq
                        do i = 1, num_fluids
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)/sum_alpha
                        end do
                    end if

                    ! Pressures relaxation procedure ===================================

                    ! Is the pressure relaxation procedure necessary?
                    relax = 1

                    !$acc loop seq
                    do i = 1, num_fluids
                        if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > (1d0 - sgm_eps)) relax = 0
                    end do

                    if (relax == 1) then
                        rho = 0d0
                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                            rho = rho + alpha_rho(i)
                            alpha0(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)
                            rhok0(i) = alpha_rho(i)/max(alpha0(i), sgm_eps)
                        end do

                        dyn_pres = 0d0

                        !$acc loop seq
                        do i = momxb, momxe
                            dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j, k, l)* &
                                    q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                        end do
                        ! Initial partial pressures
                        !$acc loop seq
                        do i = 1, num_fluids
                            pres_K_init(i) = &
                                ((q_cons_vf(i + intxb - 1)%sf(j, k, l) - &
                                    alpha_rho(i)*qvs(i) - &
                                    dyn_pres/max(rho, sgm_eps)*alpha_rho(i))/ &
                                    max(alpha0(i), sgm_eps) &
                                    - pi_infs(i))/gammas(i)

                            ck0(i) = dsqrt(gs_min(i)*(pres_K_init(i) + ps_inf(i))/max(rhok0(i),sgm_eps))
                            
                            Zk0(i) = rhok0(i)*ck0(i)
                        end do
                        if (ieee_is_nan(ck0(1)) .or. ieee_is_nan(ck0(2)) ) then
                            print *, 'ck0', ck0
                            print *, 'alpha0', alpha0
                            print *, 'alpha_rho', alpha_rho
                            print *, 'rho', rho
                            print *, 'rhok0', rhok0
                            print *, 'pres_K_init', pres_K_init
                            print *, 'p + pi', pres_K_init + ps_inf
                            call s_mpi_abort('ck0 is NaN')
                        end if
                        ! coefficients of the quadratic equation
                        ! p_I0 = max(pres_K_init(1), pres_K_init(2))
                        p_I0 = (Zk0(2)*pres_K_init(1) + Zk0(1)*pres_K_init(2))/ &
                                (Zk0(1) + Zk0(2))
                        
                        C1 = 2d0*gs_min(1)*ps_inf(1) + (gs_min(1) - 1d0)*p_I0
                        C2 = 2d0*gs_min(2)*ps_inf(2) + (gs_min(2) - 1d0)*p_I0

                        a0 = 1d0 + &
                                gs_min(1)*alpha0(2) + &
                                gs_min(2)*alpha0(1)
                        b0 = C1*alpha0(2) + C2*alpha0(1) &
                                - (1d0 + gs_min(1))*alpha0(2)*pres_K_init(2) &
                                - (1d0 + gs_min(2))*alpha0(1)*pres_K_init(1)
                        d0 = -(C2*alpha0(1)*pres_K_init(1) + &
                               C1*alpha0(2)*pres_K_init(2))

                        delta0 = b0**2 - 4d0*a0*d0
                        ! solution of the quadratic equation
                        if (b0 > 0) then
                            pres_relax = 2d0*d0/(-b0 - dsqrt(delta0))
                        else
                            pres_relax = (-b0 + dsqrt(delta0))/(2d0*a0)
                        end if

                        p_I = pres_relax
                        p_I_avg = (p_I0 + p_I)/2d0

                        ! calculate the volume fractions based on the relaxed pressure
                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha(i) = alpha0(i)* &
                                    ((gs_min(i) - 1d0)*p_I_avg + pres_K_init(i) + gs_min(i)*ps_inf(i))/ &
                                    ((gs_min(i) - 1d0)*p_I_avg + pres_relax + gs_min(i)*ps_inf(i))
                        end do

                        if (mpp_lim) then
                            sum_alpha = 0d0
                            !$acc loop seq
                            do i = 1, num_fluids
                                alpha_rho(i) = max(0d0, alpha_rho(i))
                                alpha(i) = min(max(0d0, alpha(i)), 1d0)
                                sum_alpha = sum_alpha + alpha(i)
                            end do

                            alpha = alpha/max(sum_alpha, sgm_eps)

                        end if

                        !$acc loop seq
                        do i = 1, num_fluids
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = alpha(i)
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = &
                                    alpha(i)*(gammas(i)*pres_relax + pi_infs(i)) + &
                                    alpha_rho(i)*qvs(i) + &
                                    dyn_pres*alpha_rho(i)/max(rho, sgm_eps)
                        end do
                    end if

                    do i = 1, sys_size
                        if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l)) .and. real_in_hyperbolic_step) then
                            print *, 'troubled cell caused by NAD in relaxation step', i
                            ! print *, 'relax =', relax
                            ! print *, 'a0', a0
                            ! print *, 'b0', b0
                            ! print *, 'd0', d0
                            print *, 'delta0', delta0
                            print *, 'pres_relax', pres_relax
                            print *, 'pres_K_init', pres_K_init
                            ! do q = 1, sys_size
                            !     print *, 'initial_q_cons_vf', initial_q_cons_vf(q)
                            !     print *, 'q_cons_vf', q_cons_vf(q)%sf(j, k, l)
                            ! end do
                            call s_mpi_abort('NAD criteria violated')
                        end if
                    end do
                    ! ==================================================================
                end do
            end do
        end do

    end subroutine s_pressure_relaxation_Pelanti2014 ! -----------------------

    !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_pressure_relaxation_procedure(q_cons_vf) ! ----------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf

        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        real(kind(0d0)) :: pres_relax
        real(kind(0d0)), dimension(num_fluids) :: pres_K_init
        real(kind(0d0)) :: f_pres
        real(kind(0d0)) :: df_pres
        real(kind(0d0)), dimension(num_fluids) :: rho_K_s
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho
        real(kind(0d0)), dimension(num_fluids) :: alpha
        real(kind(0d0)) :: sum_alpha
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: mQ
        real(kind(0d0)), dimension(2) :: Re
        logical :: enery_correction

        integer :: i, j, k, l, q, iter !< Generic loop iterators
        integer :: relax !< Relaxation procedure determination variable

        !$acc parallel loop collapse(3) gang vector private(pres_K_init, rho_K_s, alpha_rho, alpha, Re, pres_relax)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Numerical correction of the volume fractions
                    if (mpp_lim) then
                        sum_alpha = 0d0

                        !$acc loop seq
                        do i = 1, num_fluids
                            if ((q_cons_vf(i + contxb - 1)%sf(j, k, l) < 0d0) .or. &
                                (q_cons_vf(i + advxb - 1)%sf(j, k, l) < 0d0)) then
                                q_cons_vf(i + contxb - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 0d0
                                q_cons_vf(i + intxb - 1)%sf(j, k, l) = 0d0
                            end if

                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > 1d0) &
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = 1d0
                            sum_alpha = sum_alpha + q_cons_vf(i + advxb - 1)%sf(j, k, l)
                        end do

                        !$acc loop seq
                        do i = 1, num_fluids
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)/sum_alpha
                        end do
                    end if

                    ! Pressures relaxation procedure ===================================

                    ! Is the pressure relaxation procedure necessary?
                    relax = 1

                    ! !$acc loop seq
                    ! do i = 1, num_fluids
                    !     if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > (1d0 - sgm_eps)) relax = 0
                    ! end do

                    if (relax == 1) then
                        ! Initial state
                        pres_relax = 0d0

                        rho = 0d0
                        !$acc loop seq
                        do i = 1, num_fluids
                            rho = rho + q_cons_vf(i)%sf(j, k, l)
                        end do

                        dyn_pres = 0d0

                        !$acc loop seq
                        do i = momxb, momxe
                            dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j, k, l)* &
                                    q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                        end do


                        !$acc loop seq
                        do i = 1, num_fluids
                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) then
                                pres_K_init(i) = &
                                    ((q_cons_vf(i + intxb - 1)%sf(j, k, l) - &
                                     q_cons_vf(i)%sf(j, k, l)*qvs(i) - &
                                     dyn_pres/max(rho, sgm_eps)*q_cons_vf(i)%sf(j, k, l))/ &
                                     max(q_cons_vf(i + advxb - 1)%sf(j, k, l), sgm_eps) &
                                     - pi_infs(i))/gammas(i)

                                if (pres_K_init(i) <= -(1d0 - 1d-8)*pres_inf(i) + 1d-8) &
                                    pres_K_init(i) = -(1d0 - 1d-8)*pres_inf(i) + 1d-8
                            else
                                pres_K_init(i) = 0d0
                            end if
                            pres_relax = pres_relax + q_cons_vf(i + advxb - 1)%sf(j, k, l)*pres_K_init(i)
                        end do

                        ! Iterative process for relaxed pressure determination
                        f_pres = 1d-9
                        df_pres = 1d9

                        !$acc loop seq
                        do i = 1, num_fluids
                            rho_K_s(i) = 0d0
                        end do

                        !$acc loop seq
                        do iter = 0, 99

                            if (DABS(f_pres) > 1d-10) then
                                pres_relax = pres_relax - f_pres/df_pres

                                ! Physical pressure
                                do i = 1, num_fluids
                                    if (pres_relax <= -(1d0 - 1d-8)*pres_inf(i) + 1d-8) &
                                        pres_relax = -(1d0 - 1d-8)*pres_inf(i) + 1d0
                                end do

                                ! Newton-Raphson method
                                f_pres = -1d0
                                df_pres = 0d0

                                !$acc loop seq
                                do i = 1, num_fluids
                                    if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) then
                                        rho_K_s(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/ &
                                                     max(q_cons_vf(i + advxb - 1)%sf(j, k, l), sgm_eps) &
                                                     *((pres_relax + pres_inf(i))/(pres_K_init(i) + &
                                                                                   pres_inf(i)))**(1d0/gamma_min(i))

                                        f_pres = f_pres + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                                                 /rho_K_s(i)

                                        df_pres = df_pres - q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                                                  /(gamma_min(i)*rho_K_s(i)*(pres_relax + pres_inf(i)))
                                    end if
                                end do
                            end if

                        end do

                        ! Cell update of the volume fraction
                        !$acc loop seq
                        do i = 1, num_fluids
                            if (q_cons_vf(i + advxb - 1)%sf(j, k, l) > sgm_eps) &
                                q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                                                                       /rho_K_s(i)
                        end do
                    end if

                    ! ==================================================================
                    enery_correction = .true.
                    if (enery_correction) then
                        ! Mixture-total-energy correction ==================================

                        ! The mixture-total-energy correction of the mixture pressure P is not necessary here
                        ! because the primitive variables are directly recovered later on by the conservative
                        ! variables (see s_convert_conservative_to_primitive_variables called in s_compute_rhs).
                        ! However, the internal-energy equations should be reset with the corresponding mixture
                        ! pressure from the correction. This step is carried out below.

                        !$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                            alpha(i) = q_cons_vf(E_idx + i)%sf(j, k, l)
                        end do

                        if (bubbles) then
                            rho = 0d0
                            gamma = 0d0
                            pi_inf = 0d0

                            if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho = rho + alpha_rho(i)
                                    gamma = gamma + alpha(i)*gammas(i)
                                    pi_inf = pi_inf + alpha(i)*pi_infs(i)
                                end do
                            else if ((model_eqns == 2) .and. (num_fluids > 2)) then
                                !$acc loop seq
                                do i = 1, num_fluids - 1
                                    rho = rho + alpha_rho(i)
                                    gamma = gamma + alpha(i)*gammas(i)
                                    pi_inf = pi_inf + alpha(i)*pi_infs(i)
                                end do
                            else
                                rho = alpha_rho(1)
                                gamma = gammas(1)
                                pi_inf = pi_infs(1)
                            end if
                        else
                            rho = 0d0
                            gamma = 0d0
                            pi_inf = 0d0
                            mQ = 0d0

                            sum_alpha = 0d0

                            if (mpp_lim) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho(i) = max(0d0, alpha_rho(i))
                                    alpha(i) = min(max(0d0, alpha(i)), 1d0)
                                    sum_alpha = sum_alpha + alpha(i)
                                end do

                                alpha = alpha/max(sum_alpha, sgm_eps)

                            end if

                            !$acc loop seq
                            do i = 1, num_fluids
                                rho = rho + alpha_rho(i)
                                gamma = gamma + alpha(i)*gammas(i)
                                pi_inf = pi_inf + alpha(i)*pi_infs(i)
                                mQ = mQ + alpha_rho(i)*qvs(i)
                            end do

                            if (any(Re_size > 0)) then
                                !$acc loop seq
                                do i = 1, 2
                                    Re(i) = dflt_real

                                    if (Re_size(i) > 0) Re(i) = 0d0
                                    !$acc loop seq
                                    do q = 1, Re_size(i)
                                        Re(i) = alpha(Re_idx(i, q))/Res(i, q) &
                                                + Re(i)
                                    end do

                                    Re(i) = 1d0/max(Re(i), sgm_eps)

                                end do
                            end if
                        end if

                        dyn_pres = 0d0

                        !$acc loop seq
                        do i = momxb, momxe
                            dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j, k, l)* &
                                    q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                        end do

                        pres_relax = (q_cons_vf(E_idx)%sf(j, k, l) - dyn_pres - mQ &
                                    - pi_inf)/gamma

                        !$acc loop seq
                        do i = 1, num_fluids
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = &
                                q_cons_vf(i + advxb - 1)%sf(j, k, l)* &
                                (gammas(i)*pres_relax + pi_infs(i)) + &
                                q_cons_vf(i)%sf(j, k, l)*qvs(i) + &
                                dyn_pres*q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
                        end do
                        ! ==================================================================
                    end if
                end do
            end do
        end do

    end subroutine s_pressure_relaxation_procedure ! -----------------------

    !>  The purpose of this subroutine is to WENO-reconstruct the
        !!      left and the right cell-boundary values, including values
        !!      at the Gaussian quadrature points, from the cell-averaged
        !!      variables.
        !!  @param v_vf Cell-average variables
        !!  @param vL_qp Left WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param vR_qp Right WENO-reconstructed, cell-boundary values including
        !!          the values at the quadrature points, of the cell-average variables
        !!  @param norm_dir Splitting coordinate direction
    subroutine s_reconstruct_cell_boundary_values(v_vf, vL_x, vL_y, vL_z, vR_x, vR_y, vR_z, & ! -
                                                  norm_dir, low_order)

        type(scalar_field), dimension(iv%beg:iv%end), intent(IN) :: v_vf

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: vL_x, vL_y, vL_z, vR_x, vR_y, vR_z

        integer, intent(IN) :: norm_dir

        logical, intent(IN) :: low_order

        integer :: weno_dir !< Coordinate direction of the WENO reconstruction

        integer :: i, j, k, l
        ! Reconstruction in s1-direction ===================================

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            weno_dir = 1; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            weno_dir = 2; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        else
            is1 = iz; is2 = iy; is3 = ix
            weno_dir = 3; is1%beg = is1%beg + weno_polyn
            is1%end = is1%end - weno_polyn

        end if

        if (n > 0) then
            if (p > 0) then

                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, iv%beg:iv%end), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, iv%beg:iv%end), &
                            norm_dir, weno_dir, &
                            is1, is2, is3, low_order)
            else
                call s_weno(v_vf(iv%beg:iv%end), &
                            vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, iv%beg:iv%end), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, iv%beg:iv%end), vR_z(:, :, :, :), &
                            norm_dir, weno_dir, &
                            is1, is2, is3, low_order)
            end if
        else

            call s_weno(v_vf(iv%beg:iv%end), &
                        vL_x(:, :, :, iv%beg:iv%end), vL_y(:, :, :, :), vL_z(:, :, :, :), vR_x(:, :, :, iv%beg:iv%end), vR_y(:, :, :, :), vR_z(:, :, :, :), &
                        norm_dir, weno_dir, &
                        is1, is2, is3, low_order)
        end if

        ! ==================================================================
    end subroutine s_reconstruct_cell_boundary_values ! --------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_rhs_module() ! -----------------------------------

        integer :: i, j, k, l !< Generic loop iterators

        do j = cont_idx%beg, cont_idx%end
            !$acc exit data detach(q_prim_qp%vf(j)%sf)
            nullify (q_prim_qp%vf(j)%sf)
        end do

        do j = adv_idx%beg, adv_idx%end
            !$acc exit data detach(q_prim_qp%vf(j)%sf)
            nullify (q_prim_qp%vf(j)%sf)
        end do

        do j = mom_idx%beg, E_idx
            @:DEALLOCATE(q_cons_qp%vf(j)%sf)
            @:DEALLOCATE(q_prim_qp%vf(j)%sf)
        end do

        @:DEALLOCATE(q_cons_qp%vf, q_prim_qp%vf)
        @:DEALLOCATE_GLOBAL(qL_rsx_vf, qR_rsx_vf)
        @:DEALLOCATE_GLOBAL(lo_qL_rsx_vf, lo_qR_rsx_vf)

        if (n > 0) then
            @:DEALLOCATE_GLOBAL(qL_rsy_vf, qR_rsy_vf)
            @:DEALLOCATE_GLOBAL(lo_qL_rsy_vf, lo_qR_rsy_vf)
        end if

        if (p > 0) then
            @:DEALLOCATE_GLOBAL(qL_rsz_vf, qR_rsz_vf)
            @:DEALLOCATE_GLOBAL(lo_qL_rsz_vf, lo_qR_rsz_vf)
        end if

        if (any(Re_size > 0) .and. weno_Re_flux) then
            @:DEALLOCATE_GLOBAL(dqL_rsx_vf, dqR_rsx_vf)

            if (n > 0) then
                @:DEALLOCATE_GLOBAL(dqL_rsy_vf, dqR_rsy_vf)
            end if

            if (p > 0) then
                @:DEALLOCATE_GLOBAL(dqL_rsz_vf, dqR_rsz_vf)
            end if
        end if

        if (mpp_lim .and. bubbles) then
            !deallocate(alf_sum%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
            !$acc exit data delete(alf_sum%sf(ix%beg:ix%end, iy%beg:iy%end, iz%beg:iz%end))
        end if

        if (any(Re_size > 0)) then
            do l = mom_idx%beg, mom_idx%end
                @:DEALLOCATE(dq_prim_dx_qp(1)%vf(l)%sf)
            end do

            if (n > 0) then

                do l = mom_idx%beg, mom_idx%end
                    @:DEALLOCATE(dq_prim_dy_qp(1)%vf(l)%sf)
                end do

                if (p > 0) then
                    do l = mom_idx%beg, mom_idx%end
                        @:DEALLOCATE(dq_prim_dz_qp(1)%vf(l)%sf)
                    end do
                end if

            end if

            @:DEALLOCATE(dq_prim_dx_qp(1)%vf)
            @:DEALLOCATE(dq_prim_dy_qp(1)%vf)
            @:DEALLOCATE(dq_prim_dz_qp(1)%vf)
        end if

        if (any(Re_size > 0)) then
            do i = num_dims, 1, -1
                if (any(Re_size > 0)) then

                    do l = mom_idx%beg, mom_idx%end
                        @:DEALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf)
                        @:DEALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf)
                    end do

                    if (n > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:DEALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf)
                            @:DEALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf)
                        end do
                    end if

                    if (p > 0) then
                        do l = mom_idx%beg, mom_idx%end
                            @:DEALLOCATE(dqL_prim_dz_n(i)%vf(l)%sf)
                            @:DEALLOCATE(dqR_prim_dz_n(i)%vf(l)%sf)
                        end do
                    end if

                end if

                @:DEALLOCATE(dqL_prim_dx_n(i)%vf)
                @:DEALLOCATE(dqL_prim_dy_n(i)%vf)
                @:DEALLOCATE(dqL_prim_dz_n(i)%vf)
                @:DEALLOCATE(dqR_prim_dx_n(i)%vf)
                @:DEALLOCATE(dqR_prim_dy_n(i)%vf)
                @:DEALLOCATE(dqR_prim_dz_n(i)%vf)
            end do
        end if

        @:DEALLOCATE_GLOBAL(dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n)
        @:DEALLOCATE_GLOBAL(dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n)

        do i = num_dims, 1, -1
            if (i /= 1) then
                do l = 1, sys_size
                    nullify (cell_fluc_n(i)%vf(l)%sf)
                    nullify (inter_fluc_m_n(i)%vf(l)%sf)
                    nullify (inter_fluc_p_n(i)%vf(l)%sf)
                    nullify (lo_cell_fluc_n(i)%vf(l)%sf)
                    nullify (lo_inter_fluc_m_n(i)%vf(l)%sf)
                    nullify (lo_inter_fluc_p_n(i)%vf(l)%sf)
                    nullify (phasic_pres_src_n(i)%vf(l)%sf)
                    nullify (lo_phasic_pres_src_n(i)%vf(l)%sf)
                    nullify (flux_n(i)%vf(l)%sf)
                    nullify (flux_src_n(i)%vf(l)%sf)
                    nullify (lo_flux_n(i)%vf(l)%sf)
                    nullify (lo_flux_src_n(i)%vf(l)%sf)
                    @:DEALLOCATE(flux_gsrc_n(i)%vf(l)%sf)
                    @:DEALLOCATE(lo_flux_gsrc_n(i)%vf(l)%sf)
                end do
            else
                do l = 1, sys_size
                    @:DEALLOCATE(cell_fluc_n(i)%vf(l)%sf)
                    @:DEALLOCATE(inter_fluc_m_n(i)%vf(l)%sf)
                    @:DEALLOCATE(inter_fluc_p_n(i)%vf(l)%sf)
                    @:DEALLOCATE(lo_cell_fluc_n(i)%vf(l)%sf)
                    @:DEALLOCATE(lo_inter_fluc_m_n(i)%vf(l)%sf)
                    @:DEALLOCATE(lo_inter_fluc_p_n(i)%vf(l)%sf)
                    @:DEALLOCATE(flux_n(i)%vf(l)%sf)
                    @:DEALLOCATE(flux_gsrc_n(i)%vf(l)%sf)
                    @:DEALLOCATE(lo_flux_n(i)%vf(l)%sf)
                    @:DEALLOCATE(lo_flux_gsrc_n(i)%vf(l)%sf)
                end do

                if (any(Re_size > 0)) then
                    do l = mom_idx%beg, E_idx
                        @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                        @:DEALLOCATE(lo_flux_src_n(i)%vf(l)%sf)
                    end do
                end if

                do l = intxb, sys_size
                    @:DEALLOCATE(phasic_pres_src_n(i)%vf(l)%sf)
                    @:DEALLOCATE(lo_phasic_pres_src_n(i)%vf(l)%sf)
                end do

                if (riemann_solver == 1) then
                    do l = adv_idx%beg + 1, adv_idx%end
                        @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                        @:DEALLOCATE(lo_flux_src_n(i)%vf(l)%sf)
                    end do
                else
                    do l = adv_idx%beg + 1, adv_idx%end
                        nullify (flux_src_n(i)%vf(l)%sf)
                        nullify (lo_flux_src_n(i)%vf(l)%sf)
                    end do
                end if

                @:DEALLOCATE(flux_src_n(i)%vf(adv_idx%beg)%sf)
                @:DEALLOCATE(lo_flux_src_n(i)%vf(adv_idx%beg)%sf)
            end if

            @:DEALLOCATE(cell_fluc_n(i)%vf, inter_fluc_m_n(i)%vf, inter_fluc_p_n(i)%vf)
            @:DEALLOCATE(lo_cell_fluc_n(i)%vf, lo_inter_fluc_m_n(i)%vf, lo_inter_fluc_p_n(i)%vf)
            @:DEALLOCATE(phasic_pres_src_n(i)%vf, lo_phasic_pres_src_n(i)%vf)
            @:DEALLOCATE(flux_n(i)%vf, flux_src_n(i)%vf, flux_gsrc_n(i)%vf)
            @:DEALLOCATE(lo_flux_n(i)%vf, lo_flux_src_n(i)%vf, lo_flux_gsrc_n(i)%vf)
        end do

        @:DEALLOCATE_GLOBAL(cell_fluc_n, inter_fluc_m_n, inter_fluc_p_n)
        @:DEALLOCATE_GLOBAL(lo_cell_fluc_n, lo_inter_fluc_m_n, lo_inter_fluc_p_n)
        @:DEALLOCATE_GLOBAL(flux_n, flux_src_n, flux_gsrc_n, phasic_pres_src_n)
        @:DEALLOCATE_GLOBAL(lo_flux_n, lo_flux_src_n, lo_flux_gsrc_n, lo_phasic_pres_src_n)

        s_riemann_solver => null()
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_rhs_module ! ---------------------------------

end module m_rhs
