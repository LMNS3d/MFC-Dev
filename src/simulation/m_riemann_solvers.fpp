!>
!! @file m_riemann_solvers.f90
!! @brief Contains module m_riemann_solvers

!> @brief This module features a database of approximate and exact Riemann
!!              problem solvers for the Navier-Stokes system of equations, which
!!              is supplemented by appropriate advection equations that are used
!!              to capture the material interfaces. The closure of the system is
!!              achieved by the stiffened gas equation of state and any required
!!              mixture relations. Surface tension effects are accounted for and
!!              are modeled by means of a volume force acting across the diffuse
!!              material interface region. The implementation details of viscous
!!              and capillary effects, into the Riemann solvers, may be found in
!!              Perigaud and Saurel (2005). Note that both effects are available
!!              only in the volume fraction model. At this time, the approximate
!!              and exact Riemann solvers that are listed below are available:
!!                  1) Harten-Lax-van Leer (HLL)
!!                  2) Harten-Lax-van Leer-Contact (HLLC)
!!                  3) Exact

#:include 'macros.fpp'
#:include 'inline_riemann.fpp'
#:include 'inline_conversions.fpp'

module m_riemann_solvers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_bubbles              !< To get the bubble wall pressure function

    use ieee_arithmetic
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_riemann_solvers_module, &
 s_riemann_solver, &
 s_hll_riemann_solver, &
 s_hllc_riemann_solver, &
 s_finalize_riemann_solvers_module

    abstract interface ! =======================================================

        !> Abstract interface to the subroutines that are utilized to compute the
        !! Riemann problem solution. For additional information please reference:
        !!                        1) s_hll_riemann_solver
        !!                        2) s_hllc_riemann_solver
        !!                        3) s_exact_riemann_solver
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param dqL_prim_dx_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqL_prim_dy_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqL_prim_dz_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqR_prim_dy_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqR_prim_dz_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param gm_alphaL_vf  Left averaged gradient magnitude
        !!  @param gm_alphaR_vf Right averaged gradient magnitude
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
        subroutine s_abstract_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
                                             dqL_prim_dy_vf, &
                                             dqL_prim_dz_vf, &
                                             qL_prim_vf, &
                                             qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
                                             dqR_prim_dy_vf, &
                                             dqR_prim_dz_vf, &
                                             qR_prim_vf, &
                                             q_prim_vf, &
                                             flux_vf, flux_src_vf, &
                                             flux_gsrc_vf, &
                                             norm_dir, ix, iy, iz, &
                                             inter_fluc_m_vf, inter_fluc_p_vf, &
                                             cell_fluc_vf, &
                                             phasic_pres_src_vf, &
                                             wave_propagation)

            import :: scalar_field, int_bounds_info, sys_size, startx, starty, startz

            real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
            type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

            type(scalar_field), allocatable, dimension(:), intent(INOUT) :: qL_prim_vf, qR_prim_vf

            type(scalar_field), &
                allocatable, dimension(:), &
                intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                                 dqL_prim_dy_vf, dqR_prim_dy_vf, &
                                 dqL_prim_dz_vf, dqR_prim_dz_vf

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf, inter_fluc_m_vf, inter_fluc_p_vf, cell_fluc_vf, phasic_pres_src_vf

            integer, intent(IN) :: norm_dir

            type(int_bounds_info), intent(IN) :: ix, iy, iz
            logical, intent(IN) :: wave_propagation

        end subroutine s_abstract_riemann_solver

        !> The abstract interface to the subroutines that are utilized to compute
        !! the viscous source fluxes for either Cartesian or cylindrical geometries.
        !! For more information please refer to:
        !!      1) s_compute_cartesian_viscous_source_flux
        !!      2) s_compute_cylindrical_viscous_source_flux
        subroutine s_compute_abstract_viscous_source_flux(velL_vf, & ! -------------
                                                          dvelL_dx_vf, &
                                                          dvelL_dy_vf, &
                                                          dvelL_dz_vf, &
                                                          velR_vf, &
                                                          dvelR_dx_vf, &
                                                          dvelR_dy_vf, &
                                                          dvelR_dz_vf, &
                                                          flux_src_vf, &
                                                          norm_dir, &
                                                          ix, iy, iz)

            import :: scalar_field, int_bounds_info, num_dims, sys_size

            type(scalar_field), &
                dimension(num_dims), &
                intent(IN) :: velL_vf, velR_vf, &
                              dvelL_dx_vf, dvelR_dx_vf, &
                              dvelL_dy_vf, dvelR_dy_vf, &
                              dvelL_dz_vf, dvelR_dz_vf

            type(scalar_field), &
                dimension(sys_size), &
                intent(INOUT) :: flux_src_vf

            integer, intent(IN) :: norm_dir

            type(int_bounds_info), intent(IN) :: ix, iy, iz

        end subroutine s_compute_abstract_viscous_source_flux

    end interface ! ============================================================

    !> The cell-boundary values of the fluxes (src - source) that are computed
    !! through the chosen Riemann problem solver, and the direct evaluation of
    !! source terms, by using the left and right states given in qK_prim_rs_vf,
    !! dqK_prim_ds_vf where ds = dx, dy or dz.
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_rsx_vf, flux_src_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_rsy_vf, flux_src_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_rsz_vf, flux_src_rsz_vf)
    !$acc declare link( flux_rsx_vf, flux_src_rsx_vf, flux_rsy_vf,  &
    !$acc   flux_src_rsy_vf, flux_rsz_vf, flux_src_rsz_vf )
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_rsx_vf, flux_src_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_rsy_vf, flux_src_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_rsz_vf, flux_src_rsz_vf
    !$acc declare create( flux_rsx_vf, flux_src_rsx_vf, flux_rsy_vf,  &
    !$acc   flux_src_rsy_vf, flux_rsz_vf, flux_src_rsz_vf )
#endif
    !> @}

    !> The interfacial fluctuations of flow used in the wave propagation method,
    !! where the left-going and right-going fluctuations are stored in the
    !! minus (m) and plus (p) directions, respectively.
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), inter_fluc_rsx_m, inter_fluc_rsx_p, cell_fluc_rsx)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), inter_fluc_rsy_m, inter_fluc_rsy_p, cell_fluc_rsy)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), inter_fluc_rsz_m, inter_fluc_rsz_p, cell_fluc_rsz)
    !$acc declare link( inter_fluc_rsx_m, inter_fluc_rsx_p, inter_fluc_rsy_m, inter_fluc_rsy_p,  &
    !$acc   inter_fluc_rsz_m, inter_fluc_rsz_p, cell_fluc_rsx, cell_fluc_rsy, cell_fluc_rsz )
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: inter_fluc_rsx_m, inter_fluc_rsx_p, cell_fluc_rsx
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: inter_fluc_rsy_m, inter_fluc_rsy_p, cell_fluc_rsy
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: inter_fluc_rsz_m, inter_fluc_rsz_p, cell_fluc_rsz
    !$acc declare create( inter_fluc_rsx_m, inter_fluc_rsx_p, inter_fluc_rsy_m, inter_fluc_rsy_p,  &
    !$acc   inter_fluc_rsz_m, inter_fluc_rsz_p, cell_fluc_rsx, cell_fluc_rsy, cell_fluc_rsz )
#endif
    !> @}

    !> The cell-boundary values of the geometrical source flux that are computed
    !! through the chosen Riemann problem solver by using the left and right
    !! states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.
    !> @{
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_gsrc_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_gsrc_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), flux_gsrc_rsz_vf)
    !$acc declare link( flux_gsrc_rsx_vf, flux_gsrc_rsy_vf, flux_gsrc_rsz_vf )
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsx_vf !<
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsy_vf !<
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_gsrc_rsz_vf !<
    !$acc declare create( flux_gsrc_rsx_vf, flux_gsrc_rsy_vf, flux_gsrc_rsz_vf )
#endif
    !> @}

    ! The cell-boundary values of the velocity. vel_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the source flux.
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), vel_src_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), vel_src_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), vel_src_rsz_vf)
    !$acc declare link(vel_src_rsx_vf, vel_src_rsy_vf, vel_src_rsz_vf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: vel_src_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: vel_src_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: vel_src_rsz_vf
    !$acc declare create(vel_src_rsx_vf, vel_src_rsy_vf, vel_src_rsz_vf)
#endif

! The cell-boundary values of the phasic pressure (alphak*pk). pha_pres_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the non-conservative term of phasic total energy equation.
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), pha_pres_src_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), pha_pres_src_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), pha_pres_src_rsz_vf)
    !$acc declare link(pha_pres_src_rsx_vf, pha_pres_src_rsy_vf, pha_pres_src_rsz_vf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: pha_pres_src_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: pha_pres_src_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: pha_pres_src_rsz_vf
    !$acc declare create(pha_pres_src_rsx_vf, pha_pres_src_rsy_vf, pha_pres_src_rsz_vf)
#endif

! The cell-averaged values of the primitive variables utilized for pp limiter.
#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), q_prim_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), q_prim_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), q_prim_rsz_vf)
    !$acc declare link(q_prim_rsx_vf, q_prim_rsy_vf, q_prim_rsz_vf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsz_vf
    !$acc declare create(q_prim_rsx_vf, q_prim_rsy_vf, q_prim_rsz_vf)
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), mom_sp_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), mom_sp_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), mom_sp_rsz_vf)
    !$acc declare link(mom_sp_rsx_vf, mom_sp_rsy_vf, mom_sp_rsz_vf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: mom_sp_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: mom_sp_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: mom_sp_rsz_vf
    !$acc declare create(mom_sp_rsx_vf, mom_sp_rsy_vf, mom_sp_rsz_vf)
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), Re_avg_rsx_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), Re_avg_rsy_vf)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :), Re_avg_rsz_vf)
    !$acc declare link(Re_avg_rsx_vf, Re_avg_rsy_vf, Re_avg_rsz_vf)
#else
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: Re_avg_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: Re_avg_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: Re_avg_rsz_vf
    !$acc declare link(Re_avg_rsx_vf, Re_avg_rsy_vf, Re_avg_rsz_vf)
#endif

    procedure(s_abstract_riemann_solver), &
        pointer :: s_riemann_solver => null() !<
    !! Pointer to the procedure that is utilized to calculate either the HLL,
    !! HLLC or exact intercell fluxes, based on the choice of Riemann solver

    procedure(s_compute_abstract_viscous_source_flux), &
        pointer :: s_compute_viscous_source_flux => null() !<
    !! Pointer to the subroutine that is utilized to compute the viscous source
    !! flux for either Cartesian or cylindrical geometries.

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(int_bounds_info) :: is1, is2, is3
    type(int_bounds_info) :: isx, isy, isz
    !> @}

!$acc declare create(is1, is2, is3, isx, isy, isz)

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:),  Gs)
    !$acc declare link(Gs)
#else
    real(kind(0d0)), allocatable, dimension(:) :: Gs
    !$acc declare create(Gs)
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), Res)
    !$acc declare link(Res)
#else
    real(kind(0d0)), allocatable, dimension(:, :) :: Res
    !$acc declare create(Res)
#endif

contains

    @:s_compute_speed_of_sound()

    subroutine s_hll_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, & ! -------
                                    dqL_prim_dy_vf, &
                                    dqL_prim_dz_vf, &
                                    qL_prim_vf, &
                                    qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
                                    dqR_prim_dy_vf, &
                                    dqR_prim_dz_vf, &
                                    qR_prim_vf, &
                                    q_prim_vf, &
                                    flux_vf, flux_src_vf, &
                                    flux_gsrc_vf, &
                                    norm_dir, ix, iy, iz, &
                                    inter_fluc_m_vf, inter_fluc_p_vf, &
                                    cell_fluc_vf, &
                                    phasic_pres_src_vf, &
                                    wave_propagation)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), allocatable, dimension(:), intent(INOUT) :: qL_prim_vf, qR_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf

        ! Intercell fluxes
        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf, inter_fluc_m_vf, inter_fluc_p_vf, cell_fluc_vf, phasic_pres_src_vf

        integer, intent(IN) :: norm_dir
        type(int_bounds_info), intent(IN) :: ix, iy, iz
        logical, intent(IN) :: wave_propagation

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
        real(kind(0d0)) :: rho_L, rho_R
        real(kind(0d0)), dimension(num_dims) :: vel_L, vel_R
        real(kind(0d0)) :: pres_L, pres_R
        real(kind(0d0)) :: E_L, E_R
        real(kind(0d0)) :: H_L, H_R
        real(kind(0d0)), dimension(num_fluids) :: alpha_L, alpha_R
        real(kind(0d0)) :: Y_L, Y_R
        real(kind(0d0)) :: gamma_L, gamma_R
        real(kind(0d0)) :: pi_inf_L, pi_inf_R
        real(kind(0d0)) :: qv_L, qv_R
        real(kind(0d0)) :: c_L, c_R
        real(kind(0d0)), dimension(6) :: tau_e_L, tau_e_R
        real(kind(0d0)) :: G_L, G_R
        real(kind(0d0)), dimension(2) :: Re_L, Re_R

        real(kind(0d0)) :: rho_avg
        real(kind(0d0)), dimension(num_dims) :: vel_avg
        real(kind(0d0)) :: H_avg
        real(kind(0d0)) :: gamma_avg
        real(kind(0d0)) :: c_avg

        real(kind(0d0)) :: s_L, s_R, s_M, s_P, s_S
        real(kind(0d0)) :: xi_L, xi_R !< Left and right wave speeds functions
        real(kind(0d0)) :: xi_M, xi_P

        real(kind(0d0)) :: nbub_L, nbub_R
        real(kind(0d0)) :: ptilde_L, ptilde_R
        real(kind(0d0)) :: vel_L_rms, vel_R_rms, vel_avg_rms
        real(kind(0d0)) :: blkmod1, blkmod2
        real(kind(0d0)) :: rho_Star, E_Star, p_Star, p_K_Star, rho_K_Star
        real(kind(0d0)) :: Ms_L, Ms_R, pres_SL, pres_SR
        real(kind(0d0)) :: alpha_L_sum, alpha_R_sum

        integer :: i, j, k, l, q !< Generic loop iterators

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions
        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            qL_prim_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            qR_prim_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
        call s_initialize_riemann_solver( &
            q_prim_vf, &
            flux_vf, flux_src_vf, &
            flux_gsrc_vf, &
            norm_dir, ix, iy, iz)
        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]

            if (norm_dir == ${NORM_DIR}$) then
                !$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, vel_avg, tau_e_L, tau_e_R, G_L, G_R, Re_L, Re_R, &
                !$acc rho_avg, h_avg, gamma_avg, s_L, s_R, s_S)
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            !$acc loop seq
                            do i = 1, contxe
                                alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                            end do

                            !$acc loop seq
                            do i = 1, num_dims
                                vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                            end do

                            vel_L_rms = 0d0; vel_R_rms = 0d0

                            !$acc loop seq
                            do i = 1, num_dims
                                vel_L_rms = vel_L_rms + vel_L(i)**2d0
                                vel_R_rms = vel_R_rms + vel_R(i)**2d0
                            end do

                            !$acc loop seq
                            do i = 1, num_fluids
                                alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                            end do

                            pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                            pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                            rho_L = 0d0
                            gamma_L = 0d0
                            pi_inf_L = 0d0
                            qv_L = 0d0

                            rho_R = 0d0
                            gamma_R = 0d0
                            pi_inf_R = 0d0
                            qv_R = 0d0

                            alpha_L_sum = 0d0
                            alpha_R_sum = 0d0

                            if (mpp_lim) then
                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_L(i) = max(0d0, alpha_rho_L(i))
                                    alpha_L(i) = min(max(0d0, alpha_L(i)), 1d0)
                                    alpha_L_sum = alpha_L_sum + alpha_L(i)
                                end do

                                alpha_L = alpha_L/max(alpha_L_sum, sgm_eps)

                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_rho_R(i) = max(0d0, alpha_rho_R(i))
                                    alpha_R(i) = min(max(0d0, alpha_R(i)), 1d0)
                                    alpha_R_sum = alpha_R_sum + alpha_R(i)
                                end do

                                alpha_R = alpha_R/max(alpha_R_sum, sgm_eps)
                            end if

                            !$acc loop seq
                            do i = 1, num_fluids
                                rho_L = rho_L + alpha_rho_L(i)
                                gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)
                                qv_L = qv_L + alpha_rho_L(i)*qvs(i)

                                rho_R = rho_R + alpha_rho_R(i)
                                gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                            end do

                            if (any(Re_size > 0)) then
                                !$acc loop seq
                                do i = 1, 2
                                    Re_L(i) = dflt_real

                                    if (Re_size(i) > 0) Re_L(i) = 0d0

                                    !$acc loop seq
                                    do q = 1, Re_size(i)
                                        Re_L(i) = alpha_L(Re_idx(i, q))/Res(i, q) &
                                                  + Re_L(i)
                                    end do

                                    Re_L(i) = 1d0/max(Re_L(i), sgm_eps)

                                end do

                                !$acc loop seq
                                do i = 1, 2
                                    Re_R(i) = dflt_real

                                    if (Re_size(i) > 0) Re_R(i) = 0d0

                                    !$acc loop seq
                                    do q = 1, Re_size(i)
                                        Re_R(i) = alpha_R(Re_idx(i, q))/Res(i, q) &
                                                  + Re_R(i)
                                    end do

                                    Re_R(i) = 1d0/max(Re_R(i), sgm_eps)
                                end do
                            end if

                            E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*vel_L_rms + qv_L
                            E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*vel_R_rms + qv_R

                            H_L = (E_L + pres_L)/rho_L
                            H_R = (E_R + pres_R)/rho_R

                            if (hypoelasticity) then
                                !$acc loop seq
                                do i = 1, strxe - strxb + 1
                                    tau_e_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, strxb - 1 + i)
                                    tau_e_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, strxb - 1 + i)
                                end do

                                G_L = 0d0
                                G_R = 0d0

                                !$acc loop seq
                                do i = 1, num_fluids
                                    G_L = G_L + alpha_L(i)*Gs(i)
                                    G_R = G_R + alpha_R(i)*Gs(i)
                                end do

                                do i = 1, strxe - strxb + 1
                                    ! Elastic contribution to energy if G large enough
                                    !TODO take out if statement if stable without
                                    if ((G_L > 1000) .and. (G_R > 1000)) then
                                        E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4d0*G_L)
                                        E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4d0*G_R)
                                        ! Additional terms in 2D and 3D
                                        if ((i == 2) .or. (i == 4) .or. (i == 5)) then
                                            E_L = E_L + (tau_e_L(i)*tau_e_L(i))/(4d0*G_L)
                                            E_R = E_R + (tau_e_R(i)*tau_e_R(i))/(4d0*G_R)
                                        end if
                                    end if
                                end do
                            end if

                            @:compute_average_state()

                            call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                          vel_L_rms, c_L)

                            call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                          vel_R_rms, c_R)

                            !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                            ! variables are placeholders to call the subroutine.

                            call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                          vel_avg_rms, c_avg)

                            if (any(Re_size > 0)) then
                                !$acc loop seq
                                do i = 1, 2
                                    Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
                                end do
                            end if

                            if (wave_speeds == 1) then
                                if (hypoelasticity) then
                                    s_L = min(vel_L(dir_idx(1)) - sqrt(c_L*c_L + &
                                                                       (((4d0*G_L)/3d0) + &
                                                                        tau_e_L(dir_idx_tau(1)))/rho_L) &
                                              , vel_R(dir_idx(1)) - sqrt(c_R*c_R + &
                                                                         (((4d0*G_R)/3d0) + &
                                                                          tau_e_R(dir_idx_tau(1)))/rho_R))
                                    s_R = max(vel_R(dir_idx(1)) + sqrt(c_R*c_R + &
                                                                       (((4d0*G_R)/3d0) + &
                                                                        tau_e_R(dir_idx_tau(1)))/rho_R) &
                                              , vel_L(dir_idx(1)) + sqrt(c_L*c_L + &
                                                                         (((4d0*G_L)/3d0) + &
                                                                          tau_e_L(dir_idx_tau(1)))/rho_L))
                                else
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)
                                end if

                                s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                       (s_L - vel_L(dir_idx(1))) - &
                                       rho_R*vel_R(dir_idx(1))* &
                                       (s_R - vel_R(dir_idx(1)))) &
                                      /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                                        rho_R*(s_R - vel_R(dir_idx(1))))
                            elseif (wave_speeds == 2) then
                                pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
                                                (vel_L(dir_idx(1)) - &
                                                 vel_R(dir_idx(1))))

                                pres_SR = pres_SL

                                Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
                                                     (pres_SL/pres_L - 1d0)*pres_L/ &
                                                     ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
                                Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
                                                     (pres_SR/pres_R - 1d0)*pres_R/ &
                                                     ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

                                s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                s_S = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                            (pres_L - pres_R)/ &
                                            (rho_avg*c_avg))
                            end if

                            s_M = min(0d0, s_L); s_P = max(0d0, s_R)

                            xi_M = (5d-1 + sign(5d-1, s_L)) &
                                   + (5d-1 - sign(5d-1, s_L)) &
                                   *(5d-1 + sign(5d-1, s_R))
                            xi_P = (5d-1 - sign(5d-1, s_R)) &
                                   + (5d-1 - sign(5d-1, s_L)) &
                                   *(5d-1 + sign(5d-1, s_R))

                            ! Mass
                            !$acc loop seq
                            do i = 1, contxe
                                flux_rs${XYZ}$_vf(j, k, l, i) = &
                                    (s_M*alpha_rho_R(i)*vel_R(dir_idx(1)) &
                                     - s_P*alpha_rho_L(i)*vel_L(dir_idx(1)) &
                                     + s_M*s_P*(alpha_rho_L(i) &
                                                - alpha_rho_R(i))) &
                                    /(s_M - s_P)
                            end do

                            ! Momentum
                            if (bubbles) then
                                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(i)) &
                                              + dir_flg(dir_idx(i))*(pres_R - ptilde_R)) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(i)) &
                                                + dir_flg(dir_idx(i))*(pres_L - ptilde_L)) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                    - rho_R*vel_R(dir_idx(i)))) &
                                        /(s_M - s_P)
                                end do
                            else if (hypoelasticity) then
                                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(i)) &
                                              + dir_flg(dir_idx(i))*pres_R &
                                              - tau_e_R(dir_idx_tau(i))) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(i)) &
                                                + dir_flg(dir_idx(i))*pres_L &
                                                - tau_e_L(dir_idx_tau(i))) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                    - rho_R*vel_R(dir_idx(i)))) &
                                        /(s_M - s_P)
                                end do
                            else
                                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *vel_R(dir_idx(i)) &
                                              + dir_flg(dir_idx(i))*pres_R) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *vel_L(dir_idx(i)) &
                                                + dir_flg(dir_idx(i))*pres_L) &
                                         + s_M*s_P*(rho_L*vel_L(dir_idx(i)) &
                                                    - rho_R*vel_R(dir_idx(i)))) &
                                        /(s_M - s_P)
                                end do
                            end if

                            ! Energy
                            if (bubbles) then
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    (s_M*vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) &
                                     - s_P*vel_L(dir_idx(1))*(E_L + pres_L - ptilde_L) &
                                     + s_M*s_P*(E_L - E_R)) &
                                    /(s_M - s_P)
                            else if (hypoelasticity) then
                                !TODO: simplify this so it's not split into 3
                                if (num_dims == 1) then
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                        (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) &
                                              - (tau_e_R(dir_idx_tau(1))*vel_R(dir_idx(1)))) &
                                         - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) &
                                                - (tau_e_L(dir_idx_tau(1))*vel_L(dir_idx(1)))) &
                                         + s_M*s_P*(E_L - E_R)) &
                                        /(s_M - s_P)
                                else if (num_dims == 2) then
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                        (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) &
                                              - (tau_e_R(dir_idx_tau(1))*vel_R(dir_idx(1))) &
                                              - (tau_e_R(dir_idx_tau(2))*vel_R(dir_idx(2)))) &
                                         - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) &
                                                - (tau_e_L(dir_idx_tau(1))*vel_L(dir_idx(1))) &
                                                - (tau_e_L(dir_idx_tau(2))*vel_L(dir_idx(2)))) &
                                         + s_M*s_P*(E_L - E_R)) &
                                        /(s_M - s_P)
                                else if (num_dims == 3) then
                                    flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                        (s_M*(vel_R(dir_idx(1))*(E_R + pres_R) &
                                              - (tau_e_R(dir_idx_tau(1))*vel_R(dir_idx(1))) &
                                              - (tau_e_R(dir_idx_tau(2))*vel_R(dir_idx(2))) &
                                              - (tau_e_R(dir_idx_tau(3))*vel_R(dir_idx(3)))) &
                                         - s_P*(vel_L(dir_idx(1))*(E_L + pres_L) &
                                                - (tau_e_L(dir_idx_tau(1))*vel_L(dir_idx(1))) &
                                                - (tau_e_L(dir_idx_tau(2))*vel_L(dir_idx(2))) &
                                                - (tau_e_L(dir_idx_tau(3))*vel_L(dir_idx(3)))) &
                                         + s_M*s_P*(E_L - E_R)) &
                                        /(s_M - s_P)
                                end if
                            else
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    (s_M*vel_R(dir_idx(1))*(E_R + pres_R) &
                                     - s_P*vel_L(dir_idx(1))*(E_L + pres_L) &
                                     + s_M*s_P*(E_L - E_R)) &
                                    /(s_M - s_P)
                            end if

                            ! Elastic Stresses
                            if (hypoelasticity) then
                                do i = 1, strxe - strxb + 1 !TODO: this indexing may be slow
                                    flux_rs${XYZ}$_vf(j, k, l, strxb - 1 + i) = &
                                        (s_M*(rho_R*vel_R(dir_idx(1)) &
                                              *tau_e_R(i)) &
                                         - s_P*(rho_L*vel_L(dir_idx(1)) &
                                                *tau_e_L(i)) &
                                         + s_M*s_P*(rho_L*tau_e_L(i) &
                                                    - rho_R*tau_e_R(i))) &
                                        /(s_M - s_P)
                                end do
                            end if

                            ! Advection
                            !$acc loop seq
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf(j, k, l, i) = &
                                    (qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                     - qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)) &
                                    *s_M*s_P/(s_M - s_P)
                                flux_src_rs${XYZ}$_vf(j, k, l, i) = &
                                    (s_M*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                     - s_P*qL_prim_rs${XYZ}$_vf(j, k, l, i)) &
                                    /(s_M - s_P)
                            end do

                            ! Div(U)?
                            !$acc loop seq
                            do i = 1, num_dims
                                vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = &
                                    (xi_M*(rho_L*vel_L(dir_idx(i))* &
                                           (s_L - vel_L(dir_idx(1))) - &
                                           pres_L*dir_flg(dir_idx(i))) - &
                                     xi_P*(rho_R*vel_R(dir_idx(i))* &
                                           (s_R - vel_R(dir_idx(1))) - &
                                           pres_R*dir_flg(dir_idx(i)))) &
                                    /(xi_M*rho_L*(s_L - vel_L(dir_idx(1))) - &
                                      xi_P*rho_R*(s_R - vel_R(dir_idx(1))))
                            end do

                            if (bubbles) then
                                ! From HLLC: Kills mass transport @ bubble gas density
                                if (num_fluids > 1) then
                                    flux_rs${XYZ}$_vf(j, k, l, contxe) = 0d0
                                end if
                            end if
                        end do
                    end do
                end do
            end if

        #:endfor

        if (any(Re_size > 0)) then
            if (weno_Re_flux) then

                call s_compute_viscous_source_flux( &
                    qL_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    dqL_prim_dz_vf(momxb:momxe), &
                    qR_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    dqR_prim_dz_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy, iz)
            else
                call s_compute_viscous_source_flux( &
                    q_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    dqL_prim_dz_vf(momxb:momxe), &
                    q_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    dqR_prim_dz_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy, iz)
            end if
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz, &
                                       inter_fluc_m_vf, inter_fluc_p_vf, &
                                       cell_fluc_vf, &
                                       phasic_pres_src_vf, &
                                       wave_propagation)

    end subroutine s_hll_riemann_solver

    !> This procedure is the implementation of the Harten, Lax,
        !!      van Leer, and contact (HLLC) approximate Riemann solver,
        !!      see Toro (1999) and Johnsen (2007). The viscous and the
        !!      surface tension effects have been included by modifying
        !!      the exact Riemann solver of Perigaud and Saurel (2005).
        !!  @param qL_prim_vf The left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param dqL_prim_dx_vf The left WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqL_prim_dy_vf The left WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqL_prim_dz_vf The left WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqR_prim_dy_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqR_prim_dz_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param gm_alphaL_vf Left averaged gradient magnitude
        !!  @param gm_alphaR_vf Right averaged gradient magnitude
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_hllc_riemann_solver(qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, & ! ------
                                     dqL_prim_dy_vf, &
                                     dqL_prim_dz_vf, &
                                     qL_prim_vf, &
                                     qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
                                     dqR_prim_dy_vf, &
                                     dqR_prim_dz_vf, &
                                     qR_prim_vf, &
                                     q_prim_vf, &
                                     flux_vf, flux_src_vf, &
                                     flux_gsrc_vf, &
                                     norm_dir, ix, iy, iz, &
                                     inter_fluc_m_vf, inter_fluc_p_vf, &
                                     cell_fluc_vf, &
                                     phasic_pres_src_vf, &
                                     wave_propagation)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), allocatable, dimension(:), intent(INOUT) :: qL_prim_vf, qR_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf

        ! Intercell fluxes
        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf, inter_fluc_m_vf, inter_fluc_p_vf, cell_fluc_vf, phasic_pres_src_vf

        integer, intent(IN) :: norm_dir
        type(int_bounds_info), intent(IN) :: ix, iy, iz
        logical, intent(IN) :: wave_propagation

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_L, alpha_rho_R
        real(kind(0d0)) :: rho_L, rho_R
        real(kind(0d0)), dimension(num_dims) :: vel_L, vel_R
        real(kind(0d0)) :: pres_L, pres_R
        real(kind(0d0)) :: E_L, E_R
        real(kind(0d0)) :: H_L, H_R
        real(kind(0d0)), dimension(num_fluids) :: alpha_L, alpha_R
        real(kind(0d0)), dimension(num_fluids) :: rho_K_L, rho_K_R
        real(kind(0d0)) :: Y_L, Y_R
        real(kind(0d0)) :: gamma_L, gamma_R
        real(kind(0d0)) :: pi_inf_L, pi_inf_R
        real(kind(0d0)) :: qv_L, qv_R
        real(kind(0d0)) :: c_L, c_R
        real(kind(0d0)) :: rho_c2_L, rho_c2_R
        real(kind(0d0)), dimension(2) :: Re_L, Re_R

        real(kind(0d0)) :: rho_avg
        real(kind(0d0)), dimension(num_dims) :: vel_avg
        real(kind(0d0)) :: H_avg
        real(kind(0d0)) :: gamma_avg
        real(kind(0d0)) :: c_avg

        real(kind(0d0)) :: s_L, s_R, s_M, s_P, s_S
        real(kind(0d0)) :: xi_L, xi_R !< Left and right wave speeds functions
        real(kind(0d0)) :: xi_M, xi_P

        real(kind(0d0)) :: nbub_L, nbub_R
        real(kind(0d0)), dimension(nb) :: R0_L, R0_R
        real(kind(0d0)), dimension(nb) :: V0_L, V0_R
        real(kind(0d0)), dimension(nb) :: P0_L, P0_R
        real(kind(0d0)), dimension(nb) :: pbw_L, pbw_R
        real(kind(0d0)), dimension(nb, nmom) :: moms_L, moms_R
        real(kind(0d0)) :: ptilde_L, ptilde_R

        real(kind(0d0)) :: alpha_L_sum, alpha_R_sum, nbub_L_denom, nbub_R_denom

        real(kind(0d0)) :: PbwR3Lbar, Pbwr3Rbar
        real(kind(0d0)) :: R3Lbar, R3Rbar
        real(kind(0d0)) :: R3V2Lbar, R3V2Rbar

        real(kind(0d0)) :: vel_L_rms, vel_R_rms, vel_avg_rms
        real(kind(0d0)) :: blkmod1, blkmod2
        real(kind(0d0)) :: rho_Star, E_Star, p_Star, p_K_Star, rho_K_Star
        real(kind(0d0)) :: pres_SL, pres_SR, Ms_L, Ms_R
        real(kind(0d0)) :: start, finish
        integer :: i, j, k, l, q !< Generic loop iterators
        integer :: idx1, idxi

        !!! Variables for pp-limiter
        integer :: pp_button !< Button to activate the pp-limiter
        real(kind(0d0)) :: pp1_thetaL, pp1_thetaR, pp1_theta_mL
        real(kind(0d0)) :: pp2_thetaL, pp2_thetaR, pp2_theta_mL
        real(kind(0d0)) :: pp3_thetaL, pp3_thetaR, pp3_theta_mL
        real(kind(0d0)) :: pp4_thetaL, pp4_thetaR, pp4_theta_mL
        real(kind(0d0)) :: pp_thetaL, pp_thetaR, pp_theta_mL
        !!! Variables for one-order scheme in pp-limiter
        real(kind(0d0)), dimension(num_fluids) :: lo_alpha_rho_L, lo_alpha_rho_R
        real(kind(0d0)), dimension(num_fluids) :: lo_rho_K_L, lo_rho_K_R
        real(kind(0d0)) :: lo_rho_L, lo_rho_R
        real(kind(0d0)), dimension(num_dims) :: lo_vel_L, lo_vel_R
        real(kind(0d0)) :: lo_vel_L_rms, lo_vel_R_rms
        real(kind(0d0)) :: lo_pres_L, lo_pres_R
        real(kind(0d0)) :: lo_E_L, lo_E_R
        real(kind(0d0)) :: lo_H_L, lo_H_R
        real(kind(0d0)), dimension(num_fluids) :: lo_alpha_L, lo_alpha_R
        real(kind(0d0)), dimension(num_fluids) :: lo_Ek_L, lo_Ek_R
        real(kind(0d0)), dimension(num_fluids) :: Ek_SL, Ek_SR
        real(kind(0d0)), dimension(1:sys_size) :: lo_F_L, lo_F_R, F_HLLC
        real(kind(0d0)) :: lo_Y_L, lo_Y_R
        real(kind(0d0)) :: lo_gamma_L, lo_gamma_R
        real(kind(0d0)) :: lo_pi_inf_L, lo_pi_inf_R
        real(kind(0d0)) :: lo_qv_L, lo_qv_R
        real(kind(0d0)) :: lo_c_L, lo_c_R
        real(kind(0d0)) :: lo_rho_c2_L, lo_rho_c2_R

        !!! Variables for path-conservative HLLC RS
        real(kind(0d0)), dimension(num_fluids) :: Ek_PC_SL, Ek_PC_SR
        real(kind(0d0)), dimension(num_fluids) :: Ek_L, Ek_R
        real(kind(0d0)) :: E_PC_SL, E_PC_SR
        real(kind(0d0)), dimension(num_fluids) :: Y_PC_L, Y_PC_R, kappa_PC, Pi_PC
        real(kind(0d0)) :: u_hat_1_L, u_hat_2_L, u_hat_3_L
        real(kind(0d0)) :: u_hat_1_R, u_hat_2_R, u_hat_3_R
        integer :: m_PC
        real(kind(0d0)) :: Sigma_234_L, Sigma_234_R
        real(kind(0d0)) :: a_PC_L, b_PC_L, c_PC_L, d_PC_L
        real(kind(0d0)) :: a_PC_R, b_PC_R, c_PC_R, d_PC_R
        real(kind(0d0)) :: rhou_p, rhou_m, rho_p, rho_m
        logical :: consider_nonconservative

        !!! Variables for the cell total fluctuation in wave-propagation method
        real(kind(0d0)) :: s_L_cell, s_R_cell, s_S_cell
        real(kind(0d0)), dimension(1:sys_size) :: q_cell_L, q_cell_R
        real(kind(0d0)), dimension(1:sys_size) :: q_cell_SL, q_cell_SR
        ! cell L and R variables for the total fluctuation
        real(kind(0d0)), dimension(1:num_fluids) :: alpha_rho_L_cell, alpha_rho_R_cell
        real(kind(0d0)) :: alpha_L_sum_cell
        real(kind(0d0)) :: rho_L_cell, rho_R_cell
        real(kind(0d0)), dimension(1:num_dims) :: vel_L_cell, vel_R_cell
        real(kind(0d0)) :: pres_L_cell, pres_R_cell
        real(kind(0d0)) :: E_L_cell, E_R_cell, E_SL_cell, E_SR_cell
        real(kind(0d0)), dimension(1:num_fluids) :: Ek_SL_cell, Ek_SR_cell
        real(kind(0d0)), dimension(1:num_fluids) :: Ek_L_cell, Ek_R_cell
        real(kind(0d0)) :: H_L_cell, H_R_cell
        real(kind(0d0)), dimension(1:num_fluids) :: alpha_L_cell, alpha_R_cell
        real(kind(0d0)), dimension(1:num_fluids) :: rho_K_L_cell, rho_K_R_cell
        real(kind(0d0)) :: Y_L_cell, Y_R_cell
        real(kind(0d0)) :: gamma_L_cell, gamma_R_cell
        real(kind(0d0)) :: pi_inf_L_cell, pi_inf_R_cell
        real(kind(0d0)) :: qv_L_cell, qv_R_cell
        real(kind(0d0)) :: c_L_cell, c_R_cell
        real(kind(0d0)) :: rho_c2_L_cell, rho_c2_R_cell
        real(kind(0d0)) :: vel_L_rms_cell, vel_R_rms_cell
        real(kind(0d0)) :: xi_L_cell, xi_R_cell

        ! check phasic variables
        real(kind(0d0)), dimension(1:num_fluids) :: phasic_rho_c2_L, phasic_rho_c2_R
        real(kind(0d0)), dimension(1:num_fluids) :: lo_phasic_rho_c2_L, lo_phasic_rho_c2_R

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions

        call s_populate_riemann_states_variables_buffers( &
            qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
            dqL_prim_dy_vf, &
            dqL_prim_dz_vf, &
            qL_prim_vf, &
            qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
            dqR_prim_dy_vf, &
            dqR_prim_dz_vf, &
            qR_prim_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction

        call s_initialize_riemann_solver( &
            q_prim_vf, &
            flux_vf, flux_src_vf, &
            flux_gsrc_vf, &
            norm_dir, ix, iy, iz)

        #:for NORM_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]

            if (norm_dir == ${NORM_DIR}$) then
                if (model_eqns == 3) then
                    !ME3

                    !$acc parallel loop collapse(3) gang vector default(present) private(vel_L, vel_R, Re_L, Re_R, &
                    !$acc rho_avg, h_avg, gamma_avg, s_L, s_R, s_S, vel_avg_rms, alpha_L, alpha_R)

                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                                pres_L_cell = qR_prim_rs${XYZ}$_vf(j, k, l, E_idx)

                                gamma_L_cell = 0d0
                                pi_inf_L_cell = 0d0

                                rho_L = 0d0
                                gamma_L = 0d0
                                pi_inf_L = 0d0
                                qv_L = 0d0

                                rho_R = 0d0
                                gamma_R = 0d0
                                pi_inf_R = 0d0
                                qv_R = 0d0

                                !$acc loop seq
                                do i = 1, num_fluids

                                    gamma_L_cell = gamma_L_cell + qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                    pi_inf_L_cell = pi_inf_L_cell + qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)

                                    rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                    pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                    qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)

                                    rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                    pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                    qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                end do

                                ! variables for low-order scheme
                                lo_pres_L = q_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                lo_pres_R = q_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)
                                lo_rho_L = 0d0
                                lo_gamma_L = 0d0
                                lo_pi_inf_L = 0d0
                                lo_qv_L = 0d0
                                lo_rho_R = 0d0
                                lo_gamma_R = 0d0
                                lo_pi_inf_R = 0d0
                                lo_qv_R = 0d0

                                !$acc loop seq
                                do i = 1, num_fluids
                                    lo_alpha_L(i) = q_prim_rs${XYZ}$_vf(j, k, l, advxb + i - 1)
                                    lo_alpha_rho_L(i) = q_prim_rs${XYZ}$_vf(j, k, l, i)
                                    lo_rho_K_L(i) = lo_alpha_rho_L(i)/max(lo_alpha_L(i), sgm_eps)
                                    lo_rho_L = lo_rho_L + q_prim_rs${XYZ}$_vf(j, k, l, i)
                                    lo_gamma_L = lo_gamma_L + q_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                    lo_pi_inf_L = lo_pi_inf_L + q_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                    lo_qv_L = lo_qv_L + q_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)
                                    
                                    lo_alpha_R(i) = q_prim_rs${XYZ}$_vf(j + 1, k, l, advxb + i - 1)
                                    lo_alpha_rho_R(i) = q_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    lo_rho_K_R(i) = lo_alpha_rho_R(i)/max(lo_alpha_R(i), sgm_eps)
                                    lo_rho_R = lo_rho_R + q_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    lo_gamma_R = lo_gamma_R + q_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                    lo_pi_inf_R = lo_pi_inf_R + q_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                    lo_qv_R = lo_qv_R + q_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                end do

                                lo_vel_L_rms = 0d0; lo_vel_R_rms = 0d0

                                !$acc loop seq
                                do i = 1, num_dims
                                    lo_vel_L(i) = q_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    lo_vel_R(i) = q_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                    lo_vel_L_rms = lo_vel_L_rms + lo_vel_L(i)**2d0
                                    lo_vel_R_rms = lo_vel_R_rms + lo_vel_R(i)**2d0
                                end do
                                !$acc loop seq
                                do i = 1, num_fluids
                                    lo_Ek_L(i) = gammas(i)*lo_pres_L + pi_infs(i) + 5d-1*lo_rho_K_L(i)*lo_vel_L_rms + &
                                                    lo_rho_K_L(i)*qvs(i)
                                    lo_Ek_R(i) = gammas(i)*lo_pres_R + pi_infs(i) + 5d-1*lo_rho_K_R(i)*lo_vel_R_rms + &
                                                    lo_rho_K_R(i)*qvs(i)
                                end do

                                lo_E_L = lo_gamma_L*lo_pres_L + lo_pi_inf_L + 5d-1*lo_rho_L*lo_vel_L_rms + lo_qv_L

                                lo_E_R = lo_gamma_R*lo_pres_R + lo_pi_inf_R + 5d-1*lo_rho_R*lo_vel_R_rms + lo_qv_R

                                pp_button = 1
                                if (pp_button == 1)then

                                    pp1_thetaL = 1d0; pp1_thetaR = 1d0; pp2_thetaL = 1d0; pp2_thetaR = 1d0
                                    pp3_thetaL = 1d0; pp3_thetaR = 1d0; pp4_thetaL = 1d0; pp4_thetaR = 1d0
                                    pp1_theta_mL = 1d0; pp2_theta_mL = 1d0; pp3_theta_mL = 1d0; pp4_theta_mL = 1d0
                                    !!! check alpha1_rho1 and alpha2_rho2 =======================================
                                    do i = 1, cont_idx%end
        
                                        if (q_prim_rs${XYZ}$_vf( j ,k,l,i) <1d-11 &
                                            .or. ieee_is_nan(q_prim_rs${XYZ}$_vf( j ,k,l,i)) &
                                            )then
                                            pp1_thetaL = 0d0*pp1_thetaL 
                                        else
                                            ! if (qL_prim_rs${XYZ}$_vf(j,k,l,i) < 1d-11)then
                                            !     pp1_thetaL = (q_prim_rs${XYZ}$_vf( j ,k,l,i)-1d-11) / &
                                            !         (q_prim_rs${XYZ}$_vf( j ,k,l,i) - qL_prim_rs${XYZ}$_vf( j ,k,l,i))*pp1_thetaL 
        
                                            if (qL_prim_rs${XYZ}$_vf(j,k,l,i) < 1d-11 &
                                                .or. ieee_is_nan(qL_prim_rs${XYZ}$_vf(j,k,l,i)) &
                                                )then
        
                                                pp1_thetaL = 0d0*pp1_thetaL

                                            endif

                                            if (wave_propagation) then
                                                if (qR_prim_rs${XYZ}$_vf(j,k,l,i) < 1d-11 &
                                                    .or. ieee_is_nan(qR_prim_rs${XYZ}$_vf(j,k,l,i)) &
                                                    )then
            
                                                    pp1_theta_mL = 0d0*pp1_theta_mL

                                                endif 
                                            end if
                                        endif
        
                                        if (q_prim_rs${XYZ}$_vf( j+1 ,k,l,i) <1d-11 &
                                            .or. ieee_is_nan(q_prim_rs${XYZ}$_vf( j+1 ,k,l,i)) &
                                            )then
                                            pp1_thetaR = 0d0*pp1_thetaR
                                        else
                                            ! if (qR_prim_rs${XYZ}$_vf(j+1,k,l,i) < 1d-11)then
                                            !     pp1_thetaR = (q_prim_rs${XYZ}$_vf( j+1 ,k,l,i)-1d-11) / &
                                            !         (q_prim_rs${XYZ}$_vf( j+1 ,k,l,i) - qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i))*pp1_thetaR
        
                                            if (qR_prim_rs${XYZ}$_vf(j+1,k,l,i) < 1d-11 &
                                                .or. ieee_is_nan(qR_prim_rs${XYZ}$_vf(j+1,k,l,i)) &
                                                )then
                                                pp1_thetaR = 0d0*pp1_thetaR

                                            endif
                                        endif
                                        ! qL_prim_rs${XYZ}$_vf( j ,k,l,i) = (1-pp1_thetaL)*q_prim_rs${XYZ}$_vf( j ,k,l,i) + &
                                        !             pp1_thetaL*qL_prim_rs${XYZ}$_vf( j ,k,l,i)
                                        ! qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i) = (1-pp1_thetaR)*q_prim_rs${XYZ}$_vf( j+1 ,k,l,i) + &
                                        !             pp1_thetaR*qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i)
                                    enddo
                                    !!! --------------------------------------------------------------------------
        
        
                                    !!! check alpha1 and alpha2 =================================================
                                    do i = adv_idx%beg, adv_idx%end
        
                                        if (q_prim_rs${XYZ}$_vf( j ,k,l,i) <1d-11 &
                                            .or. ieee_is_nan(q_prim_rs${XYZ}$_vf( j ,k,l,i)) &
                                            )then
                                            pp2_thetaL = 0d0*pp2_thetaL
                                        else
                                            ! if (qL_prim_rs${XYZ}$_vf(j,k,l,i) < 1d-11)then
                                            !     pp2_thetaL = (q_prim_rs${XYZ}$_vf( j ,k,l,i)-1d-11) / &
                                            !         (q_prim_rs${XYZ}$_vf( j ,k,l,i) - qL_prim_rs${XYZ}$_vf( j ,k,l,i))*pp2_thetaL
        
                                            if (qL_prim_rs${XYZ}$_vf(j,k,l,i) < 1d-11 &
                                                .or. ieee_is_nan(qL_prim_rs${XYZ}$_vf(j,k,l,i)) &
                                                )then
                                                pp2_thetaL = 0d0*pp2_thetaL
        
                                            endif

                                            if (wave_propagation) then
                                                if (qR_prim_rs${XYZ}$_vf(j,k,l,i) < 1d-11 &
                                                    .or. ieee_is_nan(qR_prim_rs${XYZ}$_vf(j,k,l,i)) &
                                                    )then
                                                    pp2_theta_mL = 0d0*pp2_theta_mL
            
                                                endif
                                            end if
                                        endif
        
                                        if (q_prim_rs${XYZ}$_vf( j+1 ,k,l,i) < 1d-11 &
                                            .or. ieee_is_nan(q_prim_rs${XYZ}$_vf( j+1 ,k,l,i)) &
                                            )then
                                            pp2_thetaR = 0d0*pp2_thetaR
                                        else
                                            ! if (qR_prim_rs${XYZ}$_vf(j+1,k,l,i) < 1d-11)then
                                            !     pp2_thetaR = (q_prim_rs${XYZ}$_vf( j+1 ,k,l,i) -1d-11)/ &
                                            !         (q_prim_rs${XYZ}$_vf( j+1 ,k,l,i) - qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i))*pp2_thetaR
        
                                            if (qR_prim_rs${XYZ}$_vf(j+1,k,l,i) < 1d-11 &
                                                .or. ieee_is_nan(qR_prim_rs${XYZ}$_vf(j+1,k,l,i)) &
                                                )then
                                                pp2_thetaR = 0d0*pp2_thetaR
        
                                            endif
                                        endif
                                        
                                        ! qL_prim_rs${XYZ}$_vf( j ,k,l,i) = (1-pp2_thetaL)*q_prim_rs${XYZ}$_vf( j ,k,l,i) + &
                                        !             pp2_thetaL*qL_prim_rs${XYZ}$_vf( j ,k,l,i)
                                        ! qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i) = (1-pp2_thetaR)*q_prim_rs${XYZ}$_vf( j+1 ,k,l,i) + &
                                        !             pp2_thetaR*qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i)
                                    enddo
                                    !!! ----------------------------------------------------------------------------
        
        
                                    !!! check pressure ===========================================================
                                    if (q_prim_rs${XYZ}$_vf(j,k,l,E_idx) < 0d0  &
                                            .or. ieee_is_nan(q_prim_rs${XYZ}$_vf(j,k,l,E_idx)) &
                                            )then
                                        pp4_thetaL = 0d0
                                    else
                                        if (qL_prim_rs${XYZ}$_vf(j,k,l,E_idx) < 0d0 &
                                            .or. ieee_is_nan(qL_prim_rs${XYZ}$_vf(j,k,l,E_idx)) &
                                            )then
                                            ! pp4_thetaL = (q_prim_rs${XYZ}$_vf( j ,k,l,E_idx)-0d0) / &
                                            !     (q_prim_rs${XYZ}$_vf( j ,k,l,E_idx) - qL_prim_rs${XYZ}$_vf( j ,k,l,E_idx))
        
                                            pp4_thetaL = 0d0
                                        endif

                                        if (wave_propagation) then
                                            if (qR_prim_rs${XYZ}$_vf(j,k,l,E_idx) < 0d0 &
                                                .or. ieee_is_nan(qR_prim_rs${XYZ}$_vf(j,k,l,E_idx)) &
                                                )then
                                                pp4_theta_mL = 0d0
                                            endif
                                        end if
                                    endif
        
                                    if (q_prim_rs${XYZ}$_vf( j+1 ,k,l,E_idx) <0d0 &
                                            .or. ieee_is_nan(q_prim_rs${XYZ}$_vf( j+1 ,k,l,E_idx)) &
                                            )then
                                        pp4_thetaR = 0d0 
                                    else
                                        if (qR_prim_rs${XYZ}$_vf(j+1,k,l,E_idx) < 0d0 &
                                            .or. ieee_is_nan(qR_prim_rs${XYZ}$_vf(j+1,k,l,E_idx)) &
                                            )then
                                            ! pp4_thetaR = (q_prim_rs${XYZ}$_vf( j+1 ,k,l,E_idx) -0d0)/ &
                                            !     (q_prim_rs${XYZ}$_vf( j+1 ,k,l,E_idx) - qR_prim_rs${XYZ}$_vf( j+1 ,k,l,E_idx))
        
                                            pp4_thetaR = 0d0
                                        endif
                                    endif
                                    !!! ----------------------------------------------------------------------------
        
                                    !!! check rho_c2 ===============================================================
        
                                    lo_rho_c2_L = ((lo_gamma_L + 1d0)*lo_pres_L + lo_pi_inf_L)/lo_gamma_L
                                    lo_rho_c2_R = ((lo_gamma_R + 1d0)*lo_pres_R + lo_pi_inf_R)/lo_gamma_R
        
                                    rho_c2_L = ((gamma_L + 1d0)*pres_L + pi_inf_L)/gamma_L
                                    rho_c2_R = ((gamma_R + 1d0)*pres_R + pi_inf_R)/gamma_R

                                    rho_c2_L_cell = ((gamma_L_cell + 1d0)*pres_L_cell + pi_inf_L_cell)/gamma_L_cell
                                    ! print *, 'rho_c2_L_cell', rho_c2_L_cell, 'rho_c2_L', rho_c2_L, 'lo_rho_c2_L', lo_rho_c2_L, &
                                    !     'rho_c2_R', rho_c2_R, 'lo_rho_c2_R', lo_rho_c2_R

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        phasic_rho_c2_L(i) = ((gammas(i) + 1d0)*pres_L + pi_infs(i))/gammas(i)
                                        phasic_rho_c2_R(i) = ((gammas(i) + 1d0)*pres_R + pi_infs(i))/gammas(i)
                                        lo_phasic_rho_c2_L(i) = ((gammas(i) + 1d0)*lo_pres_L + pi_infs(i))/gammas(i)
                                        lo_phasic_rho_c2_R(i) = ((gammas(i) + 1d0)*lo_pres_R + pi_infs(i))/gammas(i)
                                    end do
        
                                    if (lo_rho_c2_L < 1d-9 &
                                            .or. ieee_is_nan(lo_rho_c2_L) &
                                            .or. lo_phasic_rho_c2_L(1) < 1d-9 &
                                            .or. lo_phasic_rho_c2_L(2) < 1d-9 &
                                            ) then
                                        pp3_thetaL = 0d0
                                    else
                                        if(rho_c2_L < 1d-9 &
                                            .or. ieee_is_nan(rho_c2_L) &
                                            .or. phasic_rho_c2_L(1) < 1d-9 &
                                            .or. phasic_rho_c2_L(2) < 1d-9 &
                                            ) then
                                            ! pp3_thetaL = (lo_rho_c2_L-1d-9)/(lo_rho_c2_L - rho_c2_L)
                                            pp3_thetaL = 0d0
                                        endif

                                        if (wave_propagation) then
                                            if(rho_c2_L_cell < 1d-9 &
                                                .or. ieee_is_nan(rho_c2_L_cell) &
                                                ) then
                                                pp3_theta_mL = 0d0
                                            endif
                                        end if
                                    endif
        
                                    if (lo_rho_c2_R < 1d-9 &
                                            .or. ieee_is_nan(lo_rho_c2_R) &
                                            .or. lo_phasic_rho_c2_R(1) < 1d-9 &
                                            .or. lo_phasic_rho_c2_R(2) < 1d-9 &
                                            ) then
                                        pp3_thetaR = 0d0
                                    else
                                        if(rho_c2_R < 1d-9 &
                                            .or. ieee_is_nan(rho_c2_R) &
                                            .or. phasic_rho_c2_R(1) < 1d-9 &
                                            .or. phasic_rho_c2_R(2) < 1d-9 &
                                            ) then
                                            ! pp3_thetaR = (lo_rho_c2_R-1d-9)/(lo_rho_c2_R - rho_c2_R)
                                            pp3_thetaR = 0d0
                                        endif
                                    endif
                                    !!! --------------------------------------------------------------------------------------
        
                                    !!! balance between high order and low order scheme ======================================
                                    pp_thetaL = pp1_thetaL*pp2_thetaL*pp3_thetaL*pp4_thetaL
                                    pp_thetaR = pp1_thetaR*pp2_thetaR*pp3_thetaR*pp4_thetaR
                                    pp_theta_mL = pp1_theta_mL*pp2_theta_mL*pp3_theta_mL*pp4_theta_mL
                                    !!! --------------------------------------------------------------------------------------
        
                                    ! ! make sure this limiter is not active in 1D cases
                                    ! if(n<=0)then
                                    !     pp_thetaL = 1d0
                                    !     pp_thetaR = 1d0
                                    ! endif
        
                                    do i = 1, sys_size
                                        if(pp_thetaL <1d0) then
                                            qL_prim_rs${XYZ}$_vf( j ,k,l,i) = q_prim_rs${XYZ}$_vf( j ,k,l,i)
                                        endif
                                        if(pp_thetaR <1d0) then
                                            qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i) = q_prim_rs${XYZ}$_vf( j+1 ,k,l,i)
                                        endif
                                        if(pp_theta_mL <1d0) then
                                            qR_prim_rs${XYZ}$_vf( j ,k,l,i) = q_prim_rs${XYZ}$_vf( j ,k,l,i)
                                        endif

                                        ! if (pp_theta_mL < 1d0 .or. pp_thetaL < 1d0 .or. pp_thetaR < 1d0) then
                                        !     qL_prim_rs${XYZ}$_vf( j ,k,l,i) = q_prim_rs${XYZ}$_vf( j ,k,l,i)
                                        !     qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i) = q_prim_rs${XYZ}$_vf( j+1 ,k,l,i)
                                        !     qR_prim_rs${XYZ}$_vf( j ,k,l,i) = q_prim_rs${XYZ}$_vf( j ,k,l,i)
                                        ! endif
                                        ! qL_prim_rs${XYZ}$_vf( j ,k,l,i) = (1d0 - pp_thetaL)*q_prim_rs${XYZ}$_vf( j ,k,l,i) + &
                                        !     pp_thetaL*qL_prim_rs${XYZ}$_vf( j ,k,l,i)
                                        ! qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i) = (1d0 - pp_thetaR)*q_prim_rs${XYZ}$_vf( j+1 ,k,l,i) + &
                                        !     pp_thetaR*qR_prim_rs${XYZ}$_vf( j+1 ,k,l,i)
                                    enddo
                                endif

                                vel_L_rms = 0d0; vel_R_rms = 0d0

                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2d0
                                    vel_R_rms = vel_R_rms + vel_R(i)**2d0
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                                rho_L = 0d0
                                gamma_L = 0d0
                                pi_inf_L = 0d0
                                qv_L = 0d0

                                rho_R = 0d0
                                gamma_R = 0d0
                                pi_inf_R = 0d0
                                qv_R = 0d0

                                alpha_L_sum = 0d0
                                alpha_R_sum = 0d0

                                if (mpp_lim) then
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, i) = max(0d0, qL_prim_rs${XYZ}$_vf(j, k, l, i))
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = min(max(0d0, qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)), 1d0)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)/max(alpha_L_sum, sgm_eps)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) = max(0d0, qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = min(max(0d0, qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)), 1d0)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                    pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                    qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)

                                    rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                    pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                    qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)

                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, advxb + i - 1)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, advxb + i - 1)

                                    alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxb + i - 1)
                                    alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxb + i - 1)

                                    rho_K_L(i) = alpha_rho_L(i)/max(alpha_L(i), sgm_eps)
                                    rho_K_R(i) = alpha_rho_R(i)/max(alpha_R(i), sgm_eps)
                                    kappa_PC(i) = 1d0/gammas(i)
                                    Pi_PC(i) = pi_infs(i)/gammas(i)
                                    Ek_L(i) = gammas(i)*pres_L + pi_infs(i) + 5d-1*rho_K_L(i)*vel_L_rms + &
                                                    rho_K_L(i)*qvs(i)
                                    Ek_R(i) = gammas(i)*pres_R + pi_infs(i) + 5d-1*rho_K_R(i)*vel_R_rms + &
                                                    rho_K_R(i)*qvs(i)
                                end do
                                !$acc loop seq
                                do i = 1, num_fluids
                                    Y_PC_L(i) = alpha_rho_L(i)/max(rho_L, sgm_eps)
                                    Y_PC_R(i) = alpha_rho_R(i)/max(rho_R, sgm_eps)
                                end do

                                if (any(Re_size > 0)) then
                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_L(i) = dflt_real

                                        if (Re_size(i) > 0) Re_L(i) = 0d0

                                        !$acc loop seq
                                        do q = 1, Re_size(i)
                                            Re_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_L(i)
                                        end do

                                        Re_L(i) = 1d0/max(Re_L(i), sgm_eps)

                                    end do

                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_R(i) = dflt_real

                                        if (Re_size(i) > 0) Re_R(i) = 0d0

                                        !$acc loop seq
                                        do q = 1, Re_size(i)
                                            Re_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_R(i)
                                        end do

                                        Re_R(i) = 1d0/max(Re_R(i), sgm_eps)
                                    end do
                                end if

                                E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*vel_L_rms + qv_L

                                E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*vel_R_rms + qv_R

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                              vel_L_rms, c_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                              vel_R_rms, c_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.

                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                              vel_avg_rms, c_avg)

                                if (any(Re_size > 0)) then
                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
                                    end do
                                end if

                                if (wave_speeds == 1) then
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                           (s_L - vel_L(dir_idx(1))) - &
                                           rho_R*vel_R(dir_idx(1))* &
                                           (s_R - vel_R(dir_idx(1)))) &
                                          /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                                            rho_R*(s_R - vel_R(dir_idx(1))))
                                elseif (wave_speeds == 2) then
                                    pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
                                                    (vel_L(dir_idx(1)) - &
                                                     vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
                                                         (pres_SL/pres_L - 1d0)*pres_L/ &
                                                         ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
                                    Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
                                                         (pres_SR/pres_R - 1d0)*pres_R/ &
                                                         ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                                (pres_L - pres_R)/ &
                                                (rho_avg*c_avg))
                                end if

                                ! cell L and R variables for the total fluctuation in wave-propagation method
                                if (wave_propagation) then

                                    vel_L_rms_cell = 0d0; vel_R_rms_cell = 0d0

                                    !$acc loop seq
                                    do i = 1, num_dims
                                        vel_L_cell(i) = qR_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                        vel_R_cell(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                        vel_L_rms_cell = vel_L_rms_cell + vel_L_cell(i)**2d0
                                        vel_R_rms_cell = vel_R_rms_cell + vel_R_cell(i)**2d0
                                    end do

                                    pres_L_cell = qR_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                    pres_R_cell = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)

                                    rho_L_cell = 0d0
                                    gamma_L_cell = 0d0
                                    pi_inf_L_cell = 0d0
                                    qv_L_cell = 0d0

                                    rho_R_cell = 0d0
                                    gamma_R_cell = 0d0
                                    pi_inf_R_cell = 0d0
                                    qv_R_cell = 0d0

                                    alpha_L_sum_cell = 0d0

                                    if (mpp_lim) then
                                    
                                        !$acc loop seq
                                        do i = 1, num_fluids
                                            qR_prim_rs${XYZ}$_vf(j, k, l, i) = max(0d0, qR_prim_rs${XYZ}$_vf(j, k, l, i))
                                            qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = &
                                                            min(max(0d0, qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)), 1d0)
                                            alpha_L_sum_cell = alpha_L_sum_cell + qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                        end do
    
                                        !$acc loop seq
                                        do i = 1, num_fluids
                                            qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = &
                                                            qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)/max(alpha_L_sum_cell, sgm_eps)
                                        end do
    
                                    end if

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        rho_L_cell = rho_L_cell + qR_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_L_cell = gamma_L_cell + qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                        pi_inf_L_cell = pi_inf_L_cell + qR_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                        qv_L_cell = qv_L_cell + qR_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)

                                        rho_R_cell = rho_R_cell + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_R_cell = gamma_R_cell + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                        pi_inf_R_cell = pi_inf_R_cell + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                        qv_R_cell = qv_R_cell + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)

                                        alpha_L_cell(i) = qR_prim_rs${XYZ}$_vf(j, k, l, advxb + i - 1)
                                        alpha_R_cell(i) = qL_prim_rs${XYZ}$_vf(j, k, l, advxb + i - 1)

                                        alpha_rho_L_cell(i) = qR_prim_rs${XYZ}$_vf(j, k, l, contxb + i - 1)
                                        alpha_rho_R_cell(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxb + i - 1)

                                        rho_K_L_cell(i) = alpha_rho_L_cell(i)/max(alpha_L_cell(i), sgm_eps)
                                        rho_K_R_cell(i) = alpha_rho_R_cell(i)/max(alpha_R_cell(i), sgm_eps)
                                        Ek_L_cell(i) = gammas(i)*pres_L_cell + pi_infs(i) + 5d-1*rho_K_L_cell(i)*vel_L_rms_cell + &
                                                    rho_K_L_cell(i)*qvs(i)
                                        Ek_R_cell(i) = gammas(i)*pres_R_cell + pi_infs(i) + 5d-1*rho_K_R_cell(i)*vel_R_rms_cell + &
                                                    rho_K_R_cell(i)*qvs(i)
                                    end do

                                    E_L_cell = gamma_L_cell*pres_L_cell + pi_inf_L_cell + 5d-1*rho_L_cell*vel_L_rms_cell + qv_L_cell

                                    E_R_cell = gamma_R_cell*pres_R_cell + pi_inf_R_cell + 5d-1*rho_R_cell*vel_R_rms_cell + qv_R_cell

                                    H_L_cell = (E_L_cell + pres_L_cell)/rho_L_cell
                                    H_R_cell = (E_R_cell + pres_R_cell)/rho_R_cell

                                    call s_compute_speed_of_sound(pres_L_cell, rho_L_cell, gamma_L_cell, pi_inf_L_cell, H_L_cell, &
                                                                alpha_L_cell, vel_L_rms_cell, c_L_cell)

                                    call s_compute_speed_of_sound(pres_R_cell, rho_R_cell, gamma_R_cell, pi_inf_R_cell, H_R_cell, &
                                                                alpha_R_cell, vel_R_rms_cell, c_R_cell)

                                    s_L_cell = min(vel_L_cell(dir_idx(1)) - c_L_cell, vel_R_cell(dir_idx(1)) - c_R_cell)
                                    s_R_cell = max(vel_R_cell(dir_idx(1)) + c_R_cell, vel_L_cell(dir_idx(1)) + c_L_cell)

                                    s_S_cell = (pres_R_cell - pres_L_cell + rho_L_cell*vel_L_cell(dir_idx(1))* &
                                        (s_L_cell - vel_L_cell(dir_idx(1))) - &
                                        rho_R_cell*vel_R_cell(dir_idx(1))* &
                                        (s_R_cell - vel_R_cell(dir_idx(1)))) &
                                        /(rho_L_cell*(s_L_cell - vel_L_cell(dir_idx(1))) - &
                                            rho_R_cell*(s_R_cell - vel_R_cell(dir_idx(1))))

                                end if

                                if (wave_propagation) then
                                    consider_nonconservative = .false.
                                    if (consider_nonconservative) then
                                        ! path-conservative HLLC Riemann solver
                                        ! Variables precomputed for the path-conservative HLLC RS
                                        m_PC = 2
                                        
                                        !> Compute u_hat at the left shock wave
                                        xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                        rho_Star = rho_L*xi_L
                                        rhou_p = rho_Star * s_S
                                        rhou_m = rho_L * vel_L(dir_idx(1))
                                        rho_p = rho_Star
                                        rho_m = rho_L
                                        call s_compute_PS_u_hat(u_hat_1_L, 1, m_PC, rhou_p, rhou_m, rho_p, rho_m)
                                        call s_compute_PS_u_hat(u_hat_2_L, 2, m_PC, rhou_p, rhou_m, rho_p, rho_m)
                                        call s_compute_PS_u_hat(u_hat_3_L, 3, m_PC, rhou_p, rhou_m, rho_p, rho_m)
                                        ! print *, 'u_hat_1_L', u_hat_1_L
                                        ! print *, 'u_hat_2_L', u_hat_2_L
                                        ! print *, 'u_hat_3_L', u_hat_3_L
                                        !> Compute Sigma_234 at the left shock wave
                                        Sigma_234_L = (Y_PC_L(2)*(Y_PC_L(1)*kappa_PC(2) - kappa_PC(1)*(Y_PC_L(1) - 0.5d0))*&
                                                        (xi_L - 1d0)*alpha_rho_L(1) - &
                                                        Y_PC_L(1)*(Y_PC_L(2)*kappa_PC(1) - kappa_PC(2)*(Y_PC_L(2) - 0.5d0))*&
                                                        (xi_L - 1d0)*alpha_rho_L(2))*u_hat_3_L + &
                                                        Y_PC_L(1)*Y_PC_L(2)*(kappa_PC(1) - kappa_PC(2))*(rhou_p - rhou_m)*u_hat_2_L +&
                                                        (Y_PC_L(2)*kappa_PC(1)*qvs(1)*(xi_L - 1d0)*alpha_rho_L(1) - &
                                                        Y_PC_L(1)*kappa_PC(2)*qvs(2)*(xi_L - 1d0)*alpha_rho_L(2))*u_hat_1_L
                                        a_PC_L = (Y_PC_L(1)*kappa_PC(2)*alpha_L(2)*u_hat_1_L)/(max(alpha_L(1),1d-8)*&
                                                    (s_L - (1d0 + kappa_PC(1))*s_S + Y_PC_L(2)*kappa_PC(1)*u_hat_1_L))

                                        b_PC_L = (Ek_L(1)*(s_L - vel_L(dir_idx(1)) + Y_PC_L(2)*kappa_PC(1)*u_hat_1_L) - &
                                                    s_S*(Pi_PC(1) + kappa_PC(1)*rho_K_L(1)*xi_L*(5d-1*s_S**2 + qvs(1))) - &
                                                    pres_L*vel_L(dir_idx(1)) - Ek_L(2)*Y_PC_L(1)*kappa_PC(2)*u_hat_1_L*&
                                                    alpha_L(2)/max(alpha_L(1),1d-8) + Sigma_234_L/max(alpha_L(1),1d-8))/&
                                                    (s_L - (1d0 + kappa_PC(1))*s_S + Y_PC_L(2)*kappa_PC(1)*u_hat_1_L)

                                        c_PC_L = (Y_PC_L(2)*kappa_PC(1)*alpha_L(1)*u_hat_1_L)/(max(alpha_L(2),1d-8)*&
                                                    (s_L - (1d0 + kappa_PC(2))*s_S + Y_PC_L(1)*kappa_PC(2)*u_hat_1_L))

                                        d_PC_L = (Ek_L(2)*(s_L - vel_L(dir_idx(1)) + Y_PC_L(1)*kappa_PC(2)*u_hat_1_L) - &
                                                    s_S*(Pi_PC(2) + kappa_PC(2)*rho_K_L(2)*xi_L*(5d-1*s_S**2 + qvs(2))) - &
                                                    pres_L*vel_L(dir_idx(1)) - Ek_L(1)*Y_PC_L(2)*kappa_PC(1)*u_hat_1_L*&
                                                    alpha_L(1)/max(alpha_L(2),1d-8) - Sigma_234_L/max(alpha_L(2),1d-8))/&
                                                    (s_L - (1d0 + kappa_PC(2))*s_S + Y_PC_L(1)*kappa_PC(2)*u_hat_1_L)

                                        Ek_PC_SL(1) = (a_PC_L*d_PC_L + b_PC_L)/(1d0 - a_PC_L*c_PC_L)
                                        Ek_PC_SL(2) = (c_PC_L*b_PC_L + d_PC_L)/(1d0 - a_PC_L*c_PC_L)
                                        E_PC_SL = alpha_L(1)*Ek_PC_SL(1) + alpha_L(2)*Ek_PC_SL(2)
                                        
                                        !> Compute u_hat at the right shock wave
                                        xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)
                                        rho_Star = rho_R*xi_R
                                        rhou_m = rho_Star * s_S
                                        rhou_p = rho_R * vel_R(dir_idx(1))
                                        rho_m = rho_Star
                                        rho_p = rho_R
                                        call s_compute_PS_u_hat(u_hat_1_R, 1, m_PC, rhou_p, rhou_m, rho_p, rho_m)
                                        call s_compute_PS_u_hat(u_hat_2_R, 2, m_PC, rhou_p, rhou_m, rho_p, rho_m)
                                        call s_compute_PS_u_hat(u_hat_3_R, 3, m_PC, rhou_p, rhou_m, rho_p, rho_m)
                                        !> Compute Sigma_234 at the right shock wave
                                        Sigma_234_R = (Y_PC_R(2)*(Y_PC_R(1)*kappa_PC(2) - kappa_PC(1)*(Y_PC_R(1) - 0.5d0))*&
                                                        (xi_R - 1d0)*alpha_rho_R(1) - &
                                                        Y_PC_R(1)*(Y_PC_R(2)*kappa_PC(1) - kappa_PC(2)*(Y_PC_R(2) - 0.5d0))*&
                                                        (xi_R - 1d0)*alpha_rho_R(2))*u_hat_3_R + &
                                                        Y_PC_R(1)*Y_PC_R(2)*(kappa_PC(1) - kappa_PC(2))*(rhou_m - rhou_p)*u_hat_2_R +&
                                                        (Y_PC_R(2)*kappa_PC(1)*qvs(1)*(xi_R - 1d0)*alpha_rho_R(1) - &
                                                        Y_PC_R(1)*kappa_PC(2)*qvs(2)*(xi_R - 1d0)*alpha_rho_R(2))*u_hat_1_R
                                        a_PC_R = (Y_PC_R(1)*kappa_PC(2)*alpha_R(2)*u_hat_1_R)/(max(alpha_R(1),1d-8)*&
                                                    (s_R - (1d0 + kappa_PC(1))*s_S + Y_PC_R(2)*kappa_PC(1)*u_hat_1_R))
                                        
                                        b_PC_R = (Ek_R(1)*(s_R - vel_R(dir_idx(1)) + Y_PC_R(2)*kappa_PC(1)*u_hat_1_R) - &
                                                    s_S*(Pi_PC(1) + kappa_PC(1)*rho_K_R(1)*xi_R*(5d-1*s_S**2 + qvs(1))) - &
                                                    pres_R*vel_R(dir_idx(1)) - Ek_R(2)*Y_PC_R(1)*kappa_PC(2)*u_hat_1_R*&
                                                    alpha_R(2)/max(alpha_R(1),1d-8) + Sigma_234_R/max(alpha_R(1),1d-8))/&
                                                    (s_R - (1d0 + kappa_PC(1))*s_S + Y_PC_R(2)*kappa_PC(1)*u_hat_1_R)

                                        c_PC_R = (Y_PC_R(2)*kappa_PC(1)*alpha_R(1)*u_hat_1_R)/(max(alpha_R(2),1d-8)*&
                                                    (s_R - (1d0 + kappa_PC(2))*s_S + Y_PC_R(1)*kappa_PC(2)*u_hat_1_R))

                                        d_PC_R = (Ek_R(2)*(s_R - vel_R(dir_idx(1)) + Y_PC_R(1)*kappa_PC(2)*u_hat_1_R) - &
                                                    s_S*(Pi_PC(2) + kappa_PC(2)*rho_K_R(2)*xi_R*(5d-1*s_S**2 + qvs(2))) - &
                                                    pres_R*vel_R(dir_idx(1)) - Ek_R(1)*Y_PC_R(2)*kappa_PC(1)*u_hat_1_R*&
                                                    alpha_R(1)/max(alpha_R(2),1d-8) - Sigma_234_R/max(alpha_R(2),1d-8))/&
                                                    (s_R - (1d0 + kappa_PC(2))*s_S + Y_PC_R(1)*kappa_PC(2)*u_hat_1_R)

                                        Ek_PC_SR(1) = (a_PC_R*d_PC_R + b_PC_R)/(1d0 - a_PC_R*c_PC_R)
                                        Ek_PC_SR(2) = (c_PC_R*b_PC_R + d_PC_R)/(1d0 - a_PC_R*c_PC_R)
                                        E_PC_SR = alpha_R(1)*Ek_PC_SR(1) + alpha_R(2)*Ek_PC_SR(2)

                                        ! compute the interfacial fluctuations
                                        do i = 1, cont_idx%end
                                            inter_fluc_rs${XYZ}$_m(j,k,l,i) = min(s_L,0d0)*alpha_rho_L(i)*(xi_L - 1d0) + &
                                                            min(s_R,0d0)*alpha_rho_R(i)*(1d0 - xi_R) + &
                                                            min(s_S,0d0)*(alpha_rho_R(i)*xi_R - alpha_rho_L(i)*xi_L)
                                            inter_fluc_rs${XYZ}$_p(j,k,l,i) = max(s_L,0d0)*alpha_rho_L(i)*(xi_L - 1d0) + &
                                                            max(s_R,0d0)*alpha_rho_R(i)*(1d0 - xi_R) + &
                                                            max(s_S,0d0)*(alpha_rho_R(i)*xi_R - alpha_rho_L(i)*xi_L)
                                        end do
                                        
                                        do i = mom_idx%beg, mom_idx%end
                                            q = i - mom_idx%beg + 1
                                            if (q == dir_idx(1)) then
                                                inter_fluc_rs${XYZ}$_m(j,k,l,i) = min(s_L,0d0)*(rho_L*xi_L*s_S - rho_L*vel_L(q)) + &
                                                                min(s_R,0d0)*(rho_R*vel_R(q) - rho_R*xi_R*s_S) + &
                                                                min(s_S,0d0)*(rho_R*xi_R - rho_L*xi_L)*s_S
                                                inter_fluc_rs${XYZ}$_p(j,k,l,i) = max(s_L,0d0)*(rho_L*xi_L*s_S - rho_L*vel_L(q)) + &
                                                                max(s_R,0d0)*(rho_R*vel_R(q) - rho_R*xi_R*s_S) + &
                                                                max(s_S,0d0)*(rho_R*xi_R - rho_L*xi_L)*s_S
                                            else
                                                inter_fluc_rs${XYZ}$_m(j,k,l,i) = min(s_L,0d0)*rho_L*vel_L(q)*(xi_L - 1d0) + &
                                                                min(s_R,0d0)*rho_R*vel_R(q)*(1d0 - xi_R) + &
                                                                min(s_S,0d0)*(rho_R*vel_R(q)*xi_R - rho_L*vel_L(q)*xi_L)
                                                inter_fluc_rs${XYZ}$_p(j,k,l,i) = max(s_L,0d0)*rho_L*vel_L(q)*(xi_L - 1d0) + &
                                                                max(s_R,0d0)*rho_R*vel_R(q)*(1d0 - xi_R) + &
                                                                max(s_S,0d0)*(rho_R*vel_R(q)*xi_R - rho_L*vel_L(q)*xi_L)
                                            end if
                                        end do

                                        inter_fluc_rs${XYZ}$_m(j,k,l,E_idx) = min(s_L,0d0)*(E_PC_SL - E_L) + &
                                                            min(s_R,0d0)*(E_R - E_PC_SR) + &
                                                            min(s_S,0d0)*(E_PC_SR - E_PC_SL)
                                        inter_fluc_rs${XYZ}$_p(j,k,l,E_idx) = max(s_L,0d0)*(E_PC_SL - E_L) + &
                                                            max(s_R,0d0)*(E_R - E_PC_SR) + &
                                                            max(s_S,0d0)*(E_PC_SR - E_PC_SL)

                                        do i = advxb, advxe
                                            inter_fluc_rs${XYZ}$_m(j,k,l,i) = min(s_S,0d0)*(alpha_R(i - advxb + 1) - alpha_L(i - advxb + 1))
                                            inter_fluc_rs${XYZ}$_p(j,k,l,i) = max(s_S,0d0)*(alpha_R(i - advxb + 1) - alpha_L(i - advxb + 1))
                                        end do

                                        do i = intxb, intxe
                                            inter_fluc_rs${XYZ}$_m(j,k,l,i) = min(s_L,0d0)*(Ek_PC_SL(i - intxb + 1) - Ek_L(i - intxb + 1))*&
                                                                    alpha_L(i - intxb + 1) + &
                                                            min(s_R,0d0)*(Ek_R(i - intxb + 1) - Ek_PC_SR(i - intxb + 1))*&
                                                                    alpha_R(i - intxb + 1) + &
                                                            min(s_S,0d0)*(Ek_PC_SR(i - intxb + 1)*alpha_R(i - intxb + 1) -&
                                                                        Ek_PC_SL(i - intxb + 1)*alpha_L(i - intxb + 1))
                                            inter_fluc_rs${XYZ}$_p(j,k,l,i) = max(s_L,0d0)*(Ek_PC_SL(i - intxb + 1) - Ek_L(i - intxb + 1))*&
                                                                    alpha_L(i - intxb + 1) + &
                                                            max(s_R,0d0)*(Ek_R(i - intxb + 1) - Ek_PC_SR(i - intxb + 1))*&
                                                                    alpha_R(i - intxb + 1) + &
                                                            max(s_S,0d0)*(Ek_PC_SR(i - intxb + 1)*alpha_R(i - intxb + 1) -&
                                                                        Ek_PC_SL(i - intxb + 1)*alpha_L(i - intxb + 1))
                                        end do
                                    else
                                        ! ! for the wave-propagation method without considering non-conservative variables
                                        ! ! A+dU = F(i+1) - F_HLLC, A-dU = F_HLLC - F(i)
                                        ! ! lo_F_R = F(i+1), lo_F_L = F(i)
                                        ! ! Here we first calculate lo_F_R and lo_F_L
                                        ! do i = 1, num_fluids
                                        !     ! mass
                                        !     lo_F_R(i) = lo_alpha_rho_R(i)*lo_vel_R(dir_idx(1))
                                        !     lo_F_L(i) = lo_alpha_rho_L(i)*lo_vel_L(dir_idx(1))
                                        !     ! volume fraction
                                        !     lo_F_R(i + advxb - 1) = lo_alpha_R(i)*lo_vel_R(dir_idx(1))
                                        !     lo_F_L(i + advxb - 1) = lo_alpha_L(i)*lo_vel_L(dir_idx(1))
                                        !     ! partial total energy
                                        !     lo_F_R(i + intxb - 1) = &
                                        !             lo_alpha_R(i)*(lo_Ek_R(i)+lo_pres_R)*lo_vel_R(dir_idx(1))
                                        !     lo_F_L(i + intxb - 1) = &
                                        !             lo_alpha_L(i)*(lo_Ek_L(i)+lo_pres_L)*lo_vel_L(dir_idx(1))
                                        ! end do
                                        ! ! momentum
                                        ! do i = 1, num_dims
                                        !     lo_F_R(momxb - 1 + dir_idx(i)) = &
                                        !         lo_rho_R*lo_vel_R(dir_idx(1))*lo_vel_R(dir_idx(i)) + &
                                        !         dir_flg(dir_idx(i))*lo_pres_R
                                        !     lo_F_L(momxb - 1 + dir_idx(i)) = &
                                        !         lo_rho_L*lo_vel_L(dir_idx(1))*lo_vel_L(dir_idx(i)) + &
                                        !         dir_flg(dir_idx(i))*lo_pres_L
                                        ! end do
                                        ! ! energy
                                        ! lo_F_R(E_idx) = (lo_E_R + lo_pres_R)*lo_vel_R(dir_idx(1))
                                        ! lo_F_L(E_idx) = (lo_E_L + lo_pres_L)*lo_vel_L(dir_idx(1))
                                        ! ! compute F_HLLC
                                        ! if (s_L >= 0d0) then
                                            
                                        !     do i = 1, num_fluids
                                        !         ! mass
                                        !         F_HLLC(i) = alpha_rho_L(i)*vel_L(dir_idx(1))
                                        !         ! volume fraction
                                        !         F_HLLC(i + advxb - 1) = alpha_L(i)*vel_L(dir_idx(1))
                                        !         ! partial total energy
                                        !         F_HLLC(i + intxb - 1) = alpha_L(i)*(Ek_L(i)+pres_L)*vel_L(dir_idx(1))
                                        !     end do
                                        !     ! momentum
                                        !     do i = 1, num_dims
                                        !         F_HLLC(momxb - 1 + dir_idx(i)) = &
                                        !             rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + &
                                        !             dir_flg(dir_idx(i))*pres_L
                                        !     end do
                                        !     ! energy
                                        !     F_HLLC(E_idx) = (E_L + pres_L)*vel_L(dir_idx(1))

                                        ! else if (s_R <= 0d0) then
                                            
                                        !     do i = 1, num_fluids
                                        !         ! mass
                                        !         F_HLLC(i) = alpha_rho_R(i)*vel_R(dir_idx(1))
                                        !         ! volume fraction
                                        !         F_HLLC(i + advxb - 1) = alpha_R(i)*vel_R(dir_idx(1))
                                        !         ! partial total energy
                                        !         F_HLLC(i + intxb - 1) = alpha_R(i)*(Ek_R(i)+pres_R)*vel_R(dir_idx(1))
                                        !     end do
                                        !     ! momentum
                                        !     do i = 1, num_dims
                                        !         F_HLLC(momxb - 1 + dir_idx(i)) = &
                                        !             rho_R*vel_R(dir_idx(1))*vel_R(dir_idx(i)) + &
                                        !             dir_flg(dir_idx(i))*pres_R
                                        !     end do
                                        !     ! energy
                                        !     F_HLLC(E_idx) = (E_R + pres_R)*vel_R(dir_idx(1))

                                        ! else if (s_S >= 0d0) then
                                        !     xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                            
                                        !     p_Star = rho_L*(s_L - vel_L(dir_idx(1)))*(s_S - vel_L(dir_idx(1))) + pres_L
                                        !     E_Star = 0d0
                                        !     do i = 1, num_fluids
                                        !         p_K_Star = pres_L + rho_K_L(i)*(s_L - vel_L(dir_idx(1)))*(s_S - vel_L(dir_idx(1)))
                                        !         Ek_SL(i) = xi_L*(Ek_L(i) + &
                                        !                 (s_S - vel_L(dir_idx(1)))*(rho_K_L(i)*s_S + &
                                        !                 pres_L/(s_L - vel_L(dir_idx(1)))))
                                        !         E_Star = E_Star + Ek_SL(i)*alpha_L(i)
                                        !         ! mass
                                        !         F_HLLC(i) = alpha_rho_L(i)*xi_L*s_S
                                        !         ! volume fraction
                                        !         F_HLLC(i + advxb - 1) = alpha_L(i)*s_S
                                        !         ! partial total energy
                                        !         F_HLLC(i + intxb - 1) = (Ek_SL(i)+p_K_Star)*s_S*alpha_L(i)
                                        !     end do
                                        !     ! momentum
                                        !     do i = 1, num_dims
                                        !         F_HLLC(momxb - 1 + dir_idx(i)) = &
                                        !                 rho_L*xi_L*s_S*(s_S*dir_flg(dir_idx(i)) + vel_L(dir_idx(i))* &
                                        !                     (1d0 - dir_flg(dir_idx(i)))) + dir_flg(dir_idx(i))*p_Star
                                        !     end do
                                        !     ! energy
                                        !     F_HLLC(E_idx) = (E_Star + p_Star)*s_S

                                        ! else
                                        !     xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)
                                        !     p_Star = rho_R*(s_R - vel_R(dir_idx(1)))*(s_S - vel_R(dir_idx(1))) + pres_R
                                        !     E_Star = 0d0
                                        !     do i = 1, num_fluids
                                        !         p_K_Star = pres_R + rho_K_R(i)*(s_R - vel_R(dir_idx(1)))*(s_S - vel_R(dir_idx(1)))
                                        !         Ek_SR(i) = xi_R*(Ek_R(i) + &
                                        !                 (s_S - vel_R(dir_idx(1)))*(rho_K_R(i)*s_S + &
                                        !                 pres_R/(s_R - vel_R(dir_idx(1)))))
                                        !         E_Star = E_Star + Ek_SR(i)*alpha_R(i)
                                        !         ! mass
                                        !         F_HLLC(i) = alpha_rho_R(i)*xi_R*s_S
                                        !         ! volume fraction
                                        !         F_HLLC(i + advxb - 1) = alpha_R(i)*s_S
                                        !         ! partial total energy
                                        !         F_HLLC(i + intxb - 1) = (Ek_SR(i)+p_K_Star)*s_S*alpha_R(i)
                                        !     end do
                                        !     ! momentum
                                        !     do i = 1, num_dims
                                        !         F_HLLC(momxb - 1 + dir_idx(i)) = &
                                        !                 rho_R*xi_R*s_S*(s_S*dir_flg(dir_idx(i)) + vel_R(dir_idx(i))* &
                                        !                     (1d0 - dir_flg(dir_idx(i)))) + dir_flg(dir_idx(i))*p_Star
                                        !     end do
                                        !     ! energy
                                        !     F_HLLC(E_idx) = (E_Star + p_Star)*s_S

                                        ! end if

                                        ! ! compute the interfacial fluctuations
                                        ! do i = 1, sys_size
                                        !     inter_fluc_rs${XYZ}$_m(j,k,l,i) = F_HLLC(i) - lo_F_L(i)
                                        !     inter_fluc_rs${XYZ}$_p(j,k,l,i) = lo_F_R(i) - F_HLLC(i)
                                        !     cell_fluc_rs${XYZ}$(j,k,l,i) = 0d0
                                        ! end do

                                        !> If we directly use the formula to calculte fluctuations, we can get:

                                        ! compute the interfacial fluctuations
                                        xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                        xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)
                                        do i = 1, cont_idx%end
                                            inter_fluc_rs${XYZ}$_m(j,k,l,i) = &
                                                            min(s_L,0d0)*alpha_rho_L(i)*(xi_L - 1d0) + &
                                                            min(s_R,0d0)*alpha_rho_R(i)*(1d0 - xi_R) + &
                                                            min(s_S,0d0)*(alpha_rho_R(i)*xi_R - alpha_rho_L(i)*xi_L)
                                            inter_fluc_rs${XYZ}$_p(j,k,l,i) = &
                                                            max(s_L,0d0)*alpha_rho_L(i)*(xi_L - 1d0) + &
                                                            max(s_R,0d0)*alpha_rho_R(i)*(1d0 - xi_R) + &
                                                            max(s_S,0d0)*(alpha_rho_R(i)*xi_R - alpha_rho_L(i)*xi_L)
                                        end do

                                        do i = mom_idx%beg, mom_idx%end
                                            q = i - mom_idx%beg + 1
                                            if (q == dir_idx(1)) then
                                                inter_fluc_rs${XYZ}$_m(j,k,l,i) = &
                                                                min(s_L,0d0)*(rho_L*xi_L*s_S - rho_L*vel_L(q)) + &
                                                                min(s_R,0d0)*(rho_R*vel_R(q) - rho_R*xi_R*s_S) + &
                                                                min(s_S,0d0)*(rho_R*xi_R - rho_L*xi_L)*s_S
                                                inter_fluc_rs${XYZ}$_p(j,k,l,i) = &
                                                                max(s_L,0d0)*(rho_L*xi_L*s_S - rho_L*vel_L(q)) + &
                                                                max(s_R,0d0)*(rho_R*vel_R(q) - rho_R*xi_R*s_S) + &
                                                                max(s_S,0d0)*(rho_R*xi_R - rho_L*xi_L)*s_S
                                            else
                                                inter_fluc_rs${XYZ}$_m(j,k,l,i) = &
                                                                min(s_L,0d0)*rho_L*vel_L(q)*(xi_L - 1d0) + &
                                                                min(s_R,0d0)*rho_R*vel_R(q)*(1d0 - xi_R) + &
                                                                min(s_S,0d0)*(rho_R*vel_R(q)*xi_R - rho_L*vel_L(q)*xi_L)
                                                inter_fluc_rs${XYZ}$_p(j,k,l,i) = &
                                                                max(s_L,0d0)*rho_L*vel_L(q)*(xi_L - 1d0) + &
                                                                max(s_R,0d0)*rho_R*vel_R(q)*(1d0 - xi_R) + &
                                                                max(s_S,0d0)*(rho_R*vel_R(q)*xi_R - rho_L*vel_L(q)*xi_L)
                                            end if
                                        end do
                                        
                                        E_PC_SR = 0d0
                                        E_PC_SL = 0d0
                                        do i = 1, num_fluids
                                            Ek_SL(i) = xi_L*(Ek_L(i) + &
                                                    (s_S - vel_L(dir_idx(1)))*(rho_K_L(i)*s_S + &
                                                    pres_L/(s_L - vel_L(dir_idx(1)))))
                                            E_PC_SL = E_PC_SL + Ek_SL(i)*alpha_L(i)
                                            Ek_SR(i) = xi_R*(Ek_R(i) + &
                                                    (s_S - vel_R(dir_idx(1)))*(rho_K_R(i)*s_S + &
                                                    pres_R/(s_R - vel_R(dir_idx(1)))))
                                            E_PC_SR = E_PC_SR + Ek_SR(i)*alpha_R(i)
                                        end do

                                        inter_fluc_rs${XYZ}$_m(j,k,l,E_idx) = min(s_L,0d0)*(E_PC_SL - E_L) + &
                                                            min(s_R,0d0)*(E_R - E_PC_SR) + &
                                                            min(s_S,0d0)*(E_PC_SR - E_PC_SL)
                                        inter_fluc_rs${XYZ}$_p(j,k,l,E_idx) = max(s_L,0d0)*(E_PC_SL - E_L) + &
                                                            max(s_R,0d0)*(E_R - E_PC_SR) + &
                                                            max(s_S,0d0)*(E_PC_SR - E_PC_SL)

                                        do i = advxb, advxe
                                            inter_fluc_rs${XYZ}$_m(j,k,l,i) = min(s_S,0d0)*(alpha_R(i - advxb + 1) - alpha_L(i - advxb + 1))
                                            inter_fluc_rs${XYZ}$_p(j,k,l,i) = max(s_S,0d0)*(alpha_R(i - advxb + 1) - alpha_L(i - advxb + 1))
                                        end do

                                        do i = intxb, intxe
                                            inter_fluc_rs${XYZ}$_m(j,k,l,i) = &
                                                            min(s_L,0d0)*(Ek_SL(i - intxb + 1) - Ek_L(i - intxb + 1))*&
                                                                    alpha_L(i - intxb + 1) + &
                                                            min(s_R,0d0)*(Ek_R(i - intxb + 1) - Ek_SR(i - intxb + 1))*&
                                                                    alpha_R(i - intxb + 1) + &
                                                            min(s_S,0d0)*(Ek_SR(i - intxb + 1)*alpha_R(i - intxb + 1) -&
                                                                        Ek_SL(i - intxb + 1)*alpha_L(i - intxb + 1))
                                            inter_fluc_rs${XYZ}$_p(j,k,l,i) = &
                                                            max(s_L,0d0)*(Ek_SL(i - intxb + 1) - Ek_L(i - intxb + 1))*&
                                                                    alpha_L(i - intxb + 1) + &
                                                            max(s_R,0d0)*(Ek_R(i - intxb + 1) - Ek_SR(i - intxb + 1))*&
                                                                    alpha_R(i - intxb + 1) + &
                                                            max(s_S,0d0)*(Ek_SR(i - intxb + 1)*alpha_R(i - intxb + 1) -&
                                                                        Ek_SL(i - intxb + 1)*alpha_L(i - intxb + 1))
                                        end do

                                        ! Compute the cell total fluctuations
                                        xi_L_cell = (s_L_cell - vel_L_cell(dir_idx(1)))/(s_L_cell - s_S_cell)
                                        xi_R_cell = (s_R_cell - vel_R_cell(dir_idx(1)))/(s_R_cell - s_S_cell)
                                        ! mass
                                        do i = 1, num_fluids
                                            cell_fluc_rs${XYZ}$(j,k,l,i) = s_L_cell*(alpha_rho_L_cell(i)*(xi_L_cell - 1d0)) + &
                                                            s_R_cell*(alpha_rho_R_cell(i)*(1d0 - xi_R_cell)) + &
                                                            s_S_cell*(alpha_rho_R_cell(i)*xi_R_cell - alpha_rho_L_cell(i)*xi_L_cell)
                                        end do
                                        ! momentum
                                        do i = mom_idx%beg, mom_idx%end
                                            q = i - mom_idx%beg + 1
                                            if (q == dir_idx(1)) then
                                                cell_fluc_rs${XYZ}$(j,k,l,i) = &
                                                                s_L_cell*(rho_L_cell*xi_L_cell*s_S_cell - rho_L_cell*vel_L_cell(q)) + &
                                                                s_R_cell*(rho_R_cell*vel_R_cell(q) - rho_R_cell*xi_R_cell*s_S_cell) + &
                                                                s_S_cell*(rho_R_cell*xi_R_cell - rho_L_cell*xi_L_cell)*s_S_cell
                                            else
                                                cell_fluc_rs${XYZ}$(j,k,l,i) = &
                                                                s_L_cell*rho_L_cell*vel_L_cell(q)*(xi_L_cell - 1d0) + &
                                                                s_R_cell*rho_R_cell*vel_R_cell(q)*(1d0 - xi_R_cell) + &
                                                                s_S_cell*(rho_R_cell*vel_R_cell(q)*xi_R_cell - rho_L_cell*vel_L_cell(q)*xi_L_cell)
                                            end if
                                        end do
                                        E_SR_cell = 0d0
                                        E_SL_cell = 0d0
                                        do i = 1, num_fluids
                                            Ek_SL_cell(i) = xi_L_cell*(Ek_L_cell(i) + &
                                                    (s_S_cell - vel_L_cell(dir_idx(1)))*(rho_K_L_cell(i)*s_S_cell + &
                                                    pres_L_cell/(s_L_cell - vel_L_cell(dir_idx(1)))))
                                            E_SL_cell = E_SL_cell + Ek_SL_cell(i)*alpha_L_cell(i)
                                            Ek_SR_cell(i) = xi_R_cell*(Ek_R_cell(i) + &
                                                    (s_S_cell - vel_R_cell(dir_idx(1)))*(rho_K_R_cell(i)*s_S_cell + &
                                                    pres_R_cell/(s_R_cell - vel_R_cell(dir_idx(1)))))
                                            E_SR_cell = E_SR_cell + Ek_SR_cell(i)*alpha_R_cell(i)
                                        end do
                                        ! energy
                                        cell_fluc_rs${XYZ}$(j,k,l,E_idx) = &
                                                            s_L_cell*(E_SL_cell - E_L_cell) + &
                                                            s_R_cell*(E_R_cell - E_SR_cell) + &
                                                            s_S_cell*(E_SR_cell - E_SL_cell)

                                        ! volume fraction
                                        do i = advxb, advxe
                                            cell_fluc_rs${XYZ}$(j,k,l,i) = s_S_cell*(alpha_R_cell(i - advxb + 1) - alpha_L_cell(i - advxb + 1))
                                        end do
                                        ! partial total energy
                                        do i = intxb, intxe
                                            cell_fluc_rs${XYZ}$(j,k,l,i) = &
                                                                    s_L_cell*(Ek_SL_cell(i - intxb + 1) - Ek_L_cell(i - intxb + 1))*&
                                                                    alpha_L_cell(i - intxb + 1) + &
                                                                    s_R_cell*(Ek_R_cell(i - intxb + 1) - Ek_SR_cell(i - intxb + 1))*&
                                                                    alpha_R_cell(i - intxb + 1) + &
                                                                    s_S_cell*(Ek_SR_cell(i - intxb + 1)*alpha_R_cell(i - intxb + 1) -&
                                                                        Ek_SL_cell(i - intxb + 1)*alpha_L_cell(i - intxb + 1))
                                        end do

                                    end if
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = 0.0d0
                                    
                                else

                                    ! print *, "s_L", s_L, "s_S", s_S, "s_R", s_R
                                    if (s_L >= 0d0) then
                                        p_Star = pres_L ! Only useful to recalculate the radial momentum geometric source flux
                                        !$acc loop seq
                                        do i = 1, num_fluids
                                            ! volume fraction flux alpha*u
                                            flux_rs${XYZ}$_vf(j, k, l, i + advxb - 1) = &
                                                qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)*vel_L(dir_idx(1))
                                            ! phasic mass flux alpharho*u
                                            flux_rs${XYZ}$_vf(j, k, l, i + contxb - 1) = &
                                                qL_prim_rs${XYZ}$_vf(j, k, l, i + contxb - 1)*vel_L(dir_idx(1))
                                            ! phaisc total energy flux alpha*(E+P)*u
                                            flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                                qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)* &
                                                (Ek_L(i) + pres_L)*vel_L(dir_idx(1))
                                            
                                            ! ! phasic internal energy flux
                                            ! flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                            !     (qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)* &
                                            !     (gammas(i)*pres_L + pi_infs(i)) + &
                                            !     qL_prim_rs${XYZ}$_vf(j, k, l, i + contxb - 1)* &
                                            !     qvs(i))*vel_L(dir_idx(1))
                                            ! phasic pressure source term
                                            pha_pres_src_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                                qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)*pres_L
                                        end do
                                        ! momentum flux
                                        !$acc loop seq
                                        do i = 1, num_dims
                                            flux_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(i)) = &
                                                rho_L*vel_L(dir_idx(1))*vel_L(dir_idx(i)) + dir_flg(dir_idx(i))*pres_L
                                            ! Compute the star velocities for the non-conservative terms
                                            vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = vel_L(dir_idx(i))
                                        end do
                                        ! energy flux
                                        flux_rs${XYZ}$_vf(j, k, l, E_idx) = (E_L + pres_L)*vel_L(dir_idx(1))

                                        ! Compute right solution state
                                    else if (s_R <= 0d0) then
                                        p_Star = pres_R
                                        ! Only useful to recalculate the radial momentum geometric source flux
                                        !$acc loop seq
                                        do i = 1, num_fluids
                                            ! volume fraction flux alpha*u
                                            flux_rs${XYZ}$_vf(j, k, l, i + advxb - 1) = &
                                                qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)*vel_R(dir_idx(1))
                                            ! phasic mass flux alpharho*u
                                            flux_rs${XYZ}$_vf(j, k, l, i + contxb - 1) = &
                                                qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + contxb - 1)*vel_R(dir_idx(1))
                                            ! phaisc total energy flux alpha*(E+P)*u
                                            flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                                qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)* &
                                                (Ek_R(i) + pres_R)*vel_R(dir_idx(1))
                                            ! ! phasic internal energy flux
                                            ! flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                            !     (qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)* &
                                            !     (gammas(i)*pres_R + pi_infs(i)) + &
                                            !     qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + contxb - 1)* &
                                            !     qvs(i))*vel_R(dir_idx(1))
                                            ! phasic pressure source term
                                            pha_pres_src_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                                qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)*pres_R
                                            ! print *, "right",pha_pres_src_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! momentum flux
                                        !$acc loop seq
                                        do i = 1, num_dims
                                            flux_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(i)) = &
                                                rho_R*vel_R(dir_idx(1))*vel_R(dir_idx(i)) + dir_flg(dir_idx(i))*pres_R
                                            ! Compute the star velocities for the non-conservative terms
                                            vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = vel_R(dir_idx(i))
                                            
                                        end do
                                        ! energy flux
                                        flux_rs${XYZ}$_vf(j, k, l, E_idx) = (E_R + pres_R)*vel_R(dir_idx(1))

                                        ! Compute left star solution state
                                    else if (s_S >= 0d0) then
                                        xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                        rho_Star = rho_L*xi_L
                                        p_Star = rho_L*(s_L - vel_L(dir_idx(1)))*(s_S - vel_L(dir_idx(1))) + pres_L
                                        E_Star = 0d0
                                        !$acc loop seq
                                        do i = 1, num_fluids
                                            ! ! p_star in the original code
                                            ! p_K_Star = (pres_L + pi_infs(i)/(1d0 + gammas(i)))* &
                                            !         xi_L**(1d0/gammas(i) + 1d0) - pi_infs(i)/(1d0 + gammas(i))

                                            ! !! p_star based on Saurel (2009)
                                            ! rho_K_Star = rho_K_L(i)*xi_L
                                            ! p_K_Star = (pres_L + pi_infs(i)/(1d0 + gammas(i)))* &
                                            !            (1d0/gammas(i)*rho_K_L(i) - (1d0/gammas(i) + 2d0)*rho_K_Star)/ &
                                            !            (1d0/gammas(i)*rho_K_Star - (1d0/gammas(i) + 2d0)*rho_K_L(i))- &
                                            !            pi_infs(i)/(1d0 + gammas(i))

                                            ! volume fraction flux alpha*s_S
                                            flux_rs${XYZ}$_vf(j, k, l, i + advxb - 1) = &
                                                qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)*s_S
                                            ! phasic mass flux alpharho*s_S*xi_L
                                            flux_rs${XYZ}$_vf(j, k, l, i + contxb - 1) = &
                                                qL_prim_rs${XYZ}$_vf(j, k, l, i + contxb - 1)*xi_L*s_S

                                            ! phaisc total energy flux alphak*(Ek_S+Pk_S)*s_S
                                            p_K_Star = pres_L + rho_K_L(i)*(s_L - vel_L(dir_idx(1)))*(s_S - vel_L(dir_idx(1)))
                                            Ek_SL(i) = xi_L*(Ek_L(i) + &
                                                (s_S - vel_L(dir_idx(1)))*(rho_K_L(i)*s_S + &
                                                pres_L/(s_L - vel_L(dir_idx(1)))))
                                            E_Star = E_Star + Ek_SL(i)*alpha_L(i)
                                            flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                                qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)* &
                                                (Ek_SL(i) + p_K_Star)*s_S
                                            ! !! try different approach based on Saurel (2009)
                                            ! ! phasic internal energy flux
                                            ! flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                            !     (qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)* &
                                            !     (gammas(i)*p_K_Star + pi_infs(i)) + &
                                            !     qL_prim_rs${XYZ}$_vf(j, k, l, i + contxb - 1)* &
                                            !     xi_L*qvs(i))*s_S
                                            ! phasic pressure source term
                                            pha_pres_src_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                                qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)*p_K_Star
                                            ! print *, "star_left",pha_pres_src_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! momentum flux
                                        !$acc loop seq
                                        do i = 1, num_dims
                                            flux_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(i)) = &
                                                rho_Star*s_S*(s_S*dir_flg(dir_idx(i)) + vel_L(dir_idx(i))* &
                                                            (1d0 - dir_flg(dir_idx(i)))) + dir_flg(dir_idx(i))*p_Star

                                            vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = vel_L(dir_idx(i)) + &
                                                                                        dir_flg(dir_idx(i))*(s_S - vel_L(dir_idx(i)))
                                            ! Compute the star velocities for the non-conservative terms
                                        end do
                                        ! energy flux
                                        flux_rs${XYZ}$_vf(j, k, l, E_idx) = (E_Star + p_Star)*s_S

                                        ! Compute right star solution state
                                    else
                                        xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                                        rho_Star = rho_R*xi_R

                                        p_Star = rho_R*(s_R - vel_R(dir_idx(1)))*(s_S - vel_R(dir_idx(1))) + pres_R
                                        E_Star = 0d0
                                        !$acc loop seq
                                        do i = 1, num_fluids
                                            ! ! p_star in the original code
                                            ! p_K_Star = (pres_R + pi_infs(i)/(1d0 + gammas(i)))* &
                                            !         xi_R**(1d0/gammas(i) + 1d0) - pi_infs(i)/(1d0 + gammas(i))

                                            ! !! p_star based on Saurel (2009)
                                            ! rho_K_Star = rho_K_R(i)*xi_R
                                            ! p_K_Star = (pres_R + pi_infs(i)/(1d0 + gammas(i)))* &
                                            !            (1d0/gammas(i)*rho_K_R(i) - (1d0/gammas(i) + 2d0)*rho_K_Star)/ &
                                            !            (1d0/gammas(i)*rho_K_Star - (1d0/gammas(i) + 2d0)*rho_K_R(i))- &
                                            !            pi_infs(i)/(1d0 + gammas(i))
                                            
                                            ! volume fraction flux alpha*s_S
                                            flux_rs${XYZ}$_vf(j, k, l, i + advxb - 1) = &
                                                qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)*s_S
                                            ! phasic mass flux alpharho*s_S*xi_R
                                            flux_rs${XYZ}$_vf(j, k, l, i + contxb - 1) = &
                                                qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + contxb - 1)*xi_R*s_S
                                            ! phaisc total energy flux alphak*(Ek_S+Pk_S)*s_S
                                            p_K_Star = pres_R + rho_K_R(i)*(s_R - vel_R(dir_idx(1)))*(s_S - vel_R(dir_idx(1)))
                                            Ek_SR(i) = xi_R*(Ek_R(i) + &
                                                    (s_S - vel_R(dir_idx(1)))*(rho_K_R(i)*s_S + &
                                                    pres_R/(s_R - vel_R(dir_idx(1)))))
                                            E_Star = E_Star + Ek_SR(i)*alpha_R(i)
                                            flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                                qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)* &
                                                (Ek_SR(i) + p_K_Star)*s_S
                                            ! ! phasic internal energy flux
                                            ! !! try different approach based on Saurel (2009)
                                            ! flux_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                            !     (qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)* &
                                            !     (gammas(i)*p_K_Star + pi_infs(i)) + &
                                            !     qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + contxb - 1)* &
                                            !     xi_R*qvs(i))*s_S
                                            ! phasic pressure source term
                                            pha_pres_src_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                                qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)*p_K_Star
                                            ! print *, "star_right",pha_pres_src_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! momentum flux
                                        !$acc loop seq
                                        do i = 1, num_dims
                                            flux_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(i)) = &
                                                rho_Star*s_S*(s_S*dir_flg(dir_idx(i)) + vel_R(dir_idx(i))*(1d0 - dir_flg(dir_idx(i)))) + &
                                                                                                dir_flg(dir_idx(i))*p_Star

                                            vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = vel_R(dir_idx(i)) + &
                                                                                        dir_flg(dir_idx(i))*(s_S - vel_R(dir_idx(i)))
                                            ! Compute the star velocities for the non-conservative terms
                                        end do
                                        ! energy flux
                                        flux_rs${XYZ}$_vf(j, k, l, E_idx) = (E_Star + p_Star)*s_S

                                    end if   
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                end if

                                ! compute the volume fraction and partial total energy source term
                                if (s_L >= 0d0) then
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_L(dir_idx(1))
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        pha_pres_src_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                            qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)*pres_L
                                    end do
                                else if (s_R <= 0d0) then
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_R(dir_idx(1))
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        pha_pres_src_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                            qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)*pres_R
                                    end do
                                else if (s_S >= 0d0) then
                                    xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                    ! flux_src_rs${XYZ}$_vf(j, k, l, advxb) = s_S*xi_L
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = s_S
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        p_K_Star = pres_L + rho_K_L(i)*(s_L - vel_L(dir_idx(1)))*(s_S - vel_L(dir_idx(1)))
                                        pha_pres_src_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                            qL_prim_rs${XYZ}$_vf(j, k, l, i + advxb - 1)*p_K_Star
                                    end do
                                else
                                    xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)
                                    ! flux_src_rs${XYZ}$_vf(j, k, l, advxb) = s_S*xi_R
                                    flux_src_rs${XYZ}$_vf(j, k, l, advxb) = s_S
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        p_K_Star = pres_R + rho_K_R(i)*(s_R - vel_R(dir_idx(1)))*(s_S - vel_R(dir_idx(1)))
                                        pha_pres_src_rs${XYZ}$_vf(j, k, l, i + intxb - 1) = &
                                            qR_prim_rs${XYZ}$_vf(j + 1, k, l, i + advxb - 1)*p_K_Star
                                    end do
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                if (cyl_coord .and. norm_dir == 2) then
                                    ! Substituting the advective flux into the inviscid geometrical source flux
                                    !$acc loop seq
                                    do i = 1, E_idx
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                    !$acc loop seq
                                    do i = intxb, intxe
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                    end do
                                    ! Recalculating the radial momentum geometric source flux (subtracting the pressure part)
                                    flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(1)) = &
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb - 1 + dir_idx(1)) - p_Star
                                    ! Geometrical source of the void fraction(s) is zero
                                    !$acc loop seq
                                    do i = advxb, advxe
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0d0
                                    end do
                                end if

                            end do
                        end do
                    end do
                elseif (model_eqns == 4) then
                    !ME4
                    !$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, vel_avg, &
                    !$acc rho_avg, h_avg, gamma_avg, s_L, s_R, s_S, vel_avg_rms, nbub_L, nbub_R, ptilde_L, ptilde_R)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                !$acc loop seq
                                do i = 1, contxe
                                    alpha_rho_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    alpha_rho_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                end do

                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                end do

                                vel_L_rms = 0d0; vel_R_rms = 0d0
                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L_rms = vel_L_rms + vel_L(i)**2d0
                                    vel_R_rms = vel_R_rms + vel_R(i)**2d0
                                end do

                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                                rho_L = 0d0
                                gamma_L = 0d0
                                pi_inf_L = 0d0
                                qv_L = 0d0
                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_L = rho_L + alpha_rho_L(i)
                                    gamma_L = gamma_L + alpha_L(i)*gammas(i)
                                    pi_inf_L = pi_inf_L + alpha_L(i)*pi_infs(i)
                                    qv_L = qv_L + alpha_rho_L(i)*qvs(i)
                                end do

                                rho_R = 0d0
                                gamma_R = 0d0
                                pi_inf_R = 0d0
                                qv_R = 0d0
                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_R = rho_R + alpha_rho_R(i)
                                    gamma_R = gamma_R + alpha_R(i)*gammas(i)
                                    pi_inf_R = pi_inf_R + alpha_R(i)*pi_infs(i)
                                    qv_R = qv_R + alpha_rho_R(i)*qvs(i)
                                end do

                                E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*vel_L_rms + qv_L

                                E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*vel_R_rms + qv_R

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                              vel_L_rms, c_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                              vel_R_rms, c_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.

                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                              vel_avg_rms, c_avg)

                                if (wave_speeds == 1) then
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                           (s_L - vel_L(dir_idx(1))) - &
                                           rho_R*vel_R(dir_idx(1))* &
                                           (s_R - vel_R(dir_idx(1)))) &
                                          /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                                            rho_R*(s_R - vel_R(dir_idx(1))))
                                elseif (wave_speeds == 2) then
                                    pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
                                                    (vel_L(dir_idx(1)) - &
                                                     vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
                                                         (pres_SL/pres_L - 1d0)*pres_L/ &
                                                         ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
                                    Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
                                                         (pres_SR/pres_R - 1d0)*pres_R/ &
                                                         ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                                (pres_L - pres_R)/ &
                                                (rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M = min(0d0, s_L); s_P = max(0d0, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                                ! goes with numerical velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5d-1 + sign(5d-1, s_S))
                                xi_P = (5d-1 - sign(5d-1, s_S))

                                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*alpha_rho_L(i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                        + xi_P*alpha_rho_R(i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                end do

                                ! Momentum flux.
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                     vel_L(dir_idx(i)) + &
                                                     s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + &
                                                                (1d0 - dir_flg(dir_idx(i)))* &
                                                                vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) + &
                                              dir_flg(dir_idx(i))*pres_L) &
                                        + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                       vel_R(dir_idx(i)) + &
                                                       s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + &
                                                                  (1d0 - dir_flg(dir_idx(i)))* &
                                                                  vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) + &
                                                dir_flg(dir_idx(i))*pres_R)
                                end do

                                if (bubbles) then
                                    ! Put p_tilde in
                                    !$acc loop seq
                                    do i = 1, num_dims
                                        flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                            flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) + &
                                            xi_M*(dir_flg(dir_idx(i))*(-1d0*ptilde_L)) &
                                            + xi_P*(dir_flg(dir_idx(i))*(-1d0*ptilde_R))
                                    end do
                                end if

                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = 0.d0

                                !$acc loop seq
                                do i = alf_idx, alf_idx !only advect the void fraction
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                end do

                                ! Source for volume fraction advection equation
                                !$acc loop seq
                                do i = 1, num_dims

                                    vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = 0d0
                                    !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0d0
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                ! Add advection flux for bubble variables
                                if (bubbles) then
                                    !$acc loop seq
                                    do i = bubxb, bubxe
                                        flux_rs${XYZ}$_vf(j, k, l, i) = &
                                            xi_M*nbub_L*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                            *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                            + xi_P*nbub_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                            *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                    end do
                                end if

                                ! Geometrical source flux for cylindrical coordinates

                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        !$acc loop seq
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = &
                                            xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                         vel_L(dir_idx(1)) + &
                                                         s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + &
                                                                    (1d0 - dir_flg(dir_idx(1)))* &
                                                                    vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                            + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                           vel_R(dir_idx(1)) + &
                                                           s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + &
                                                                      (1d0 - dir_flg(dir_idx(1)))* &
                                                                      vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        !$acc loop seq
                                        do i = advxb, advxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0d0
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        !$acc loop seq
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0d0
                                        end do
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb + 1) = &
                                            -xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                          vel_L(dir_idx(1)) + &
                                                          s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + &
                                                                     (1d0 - dir_flg(dir_idx(1)))* &
                                                                     vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                            - xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                           vel_R(dir_idx(1)) + &
                                                           s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + &
                                                                      (1d0 - dir_flg(dir_idx(1)))* &
                                                                      vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxe) = flux_rs${XYZ}$_vf(j, k, l, momxb + 1)
                                    end if
                                #:endif
                            end do
                        end do
                    end do

                elseif (model_eqns == 2 .and. bubbles) then
                    !$acc parallel loop collapse(3) gang vector default(present) private(R0_L, R0_R, V0_L, V0_R, P0_L, P0_R, pbw_L, pbw_R, vel_L, vel_R, &
                    !$acc rho_avg, alpha_L, alpha_R, h_avg, gamma_avg, s_L, s_R, s_S, nbub_L, nbub_R, ptilde_L, ptilde_R, vel_avg_rms, Re_L, Re_R)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                end do

                                vel_L_rms = 0d0; vel_R_rms = 0d0

                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2d0
                                    vel_R_rms = vel_R_rms + vel_R(i)**2d0
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                                rho_L = 0d0
                                gamma_L = 0d0
                                pi_inf_L = 0d0
                                qv_L = 0d0

                                if (mpp_lim .and. (num_fluids > 2)) then
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                        pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                        qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)
                                    end do
                                else if (num_fluids > 2) then
                                    !$acc loop seq
                                    do i = 1, num_fluids - 1
                                        rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                        gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                        pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                        qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)
                                    end do
                                else
                                    rho_L = qL_prim_rs${XYZ}$_vf(j, k, l, 1)
                                    gamma_L = gammas(1)
                                    pi_inf_L = pi_infs(1)
                                    qv_L = qvs(1)
                                end if

                                rho_R = 0d0
                                gamma_R = 0d0
                                pi_inf_R = 0d0
                                qv_R = 0d0

                                if (mpp_lim .and. (num_fluids > 2)) then
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                        gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                        pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                        qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                    end do
                                else if (num_fluids > 2) then
                                    !$acc loop seq
                                    do i = 1, num_fluids - 1
                                        rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                        gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                        pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                        qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                    end do
                                else
                                    rho_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, 1)
                                    gamma_R = gammas(1)
                                    pi_inf_R = pi_infs(1)
                                    qv_R = qvs(1)
                                end if

                                if (any(Re_size > 0)) then
                                    if (num_fluids == 1) then ! Need to consider case with num_fluids >= 2
                                        !$acc loop seq
                                        do i = 1, 2
                                            Re_L(i) = dflt_real

                                            if (Re_size(i) > 0) Re_L(i) = 0d0

                                            !$acc loop seq
                                            do q = 1, Re_size(i)
                                                Re_L(i) = (1d0 - qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q)))/Res(i, q) &
                                                          + Re_L(i)
                                            end do

                                            Re_L(i) = 1d0/max(Re_L(i), sgm_eps)

                                        end do

                                        !$acc loop seq
                                        do i = 1, 2
                                            Re_R(i) = dflt_real

                                            if (Re_size(i) > 0) Re_R(i) = 0d0

                                            !$acc loop seq
                                            do q = 1, Re_size(i)
                                                Re_R(i) = (1d0 - qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q)))/Res(i, q) &
                                                          + Re_R(i)
                                            end do

                                            Re_R(i) = 1d0/max(Re_R(i), sgm_eps)
                                        end do
                                    end if
                                end if

                                E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*vel_L_rms

                                E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*vel_R_rms

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R
                                if (avg_state == 2) then
                                    !$acc loop seq
                                    do i = 1, nb
                                        R0_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, rs(i))
                                        R0_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, rs(i))

                                        V0_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, vs(i))
                                        V0_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, vs(i))
                                        if (.not. polytropic .and. .not. qbmm) then
                                            P0_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, ps(i))
                                            P0_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, ps(i))
                                        end if
                                    end do

                                    if (.not. qbmm) then
                                        if (adv_n) then
                                            nbub_L = qL_prim_rs${XYZ}$_vf(j, k, l, n_idx)
                                            nbub_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, n_idx)
                                        else
                                            nbub_L_denom = 0d0
                                            nbub_R_denom = 0d0
                                            !$acc loop seq
                                            do i = 1, nb
                                                nbub_L_denom = nbub_L_denom + (R0_L(i)**3d0)*weight(i)
                                                nbub_R_denom = nbub_R_denom + (R0_R(i)**3d0)*weight(i)
                                            end do
                                            nbub_L = (3.d0/(4.d0*pi))*qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + num_fluids)/nbub_L_denom
                                            nbub_R = (3.d0/(4.d0*pi))*qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + num_fluids)/nbub_R_denom
                                        end if
                                    else
                                        !nb stored in 0th moment of first R0 bin in variable conversion module
                                        nbub_L = qL_prim_rs${XYZ}$_vf(j, k, l, bubxb)
                                        nbub_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, bubxb)
                                    end if

                                    !$acc loop seq
                                    do i = 1, nb
                                        if (.not. qbmm) then
                                            if (polytropic) then
                                                pbw_L(i) = f_cpbw_KM(R0(i), R0_L(i), V0_L(i), 0d0)
                                                pbw_R(i) = f_cpbw_KM(R0(i), R0_R(i), V0_R(i), 0d0)
                                            else
                                                pbw_L(i) = f_cpbw_KM(R0(i), R0_L(i), V0_L(i), P0_L(i))
                                                pbw_R(i) = f_cpbw_KM(R0(i), R0_R(i), V0_R(i), P0_R(i))
                                            end if
                                        end if
                                    end do

                                    if (qbmm) then
                                        PbwR3Lbar = mom_sp_rs${XYZ}$_vf(j, k, l, 4)
                                        PbwR3Rbar = mom_sp_rs${XYZ}$_vf(j + 1, k, l, 4)

                                        R3Lbar = mom_sp_rs${XYZ}$_vf(j, k, l, 1)
                                        R3Rbar = mom_sp_rs${XYZ}$_vf(j + 1, k, l, 1)

                                        R3V2Lbar = mom_sp_rs${XYZ}$_vf(j, k, l, 3)
                                        R3V2Rbar = mom_sp_rs${XYZ}$_vf(j + 1, k, l, 3)
                                    else

                                        PbwR3Lbar = 0d0
                                        PbwR3Rbar = 0d0

                                        R3Lbar = 0d0
                                        R3Rbar = 0d0

                                        R3V2Lbar = 0d0
                                        R3V2Rbar = 0d0

                                        !$acc loop seq
                                        do i = 1, nb
                                            PbwR3Lbar = PbwR3Lbar + pbw_L(i)*(R0_L(i)**3.d0)*weight(i)
                                            PbwR3Rbar = PbwR3Rbar + pbw_R(i)*(R0_R(i)**3.d0)*weight(i)

                                            R3Lbar = R3Lbar + (R0_L(i)**3.d0)*weight(i)
                                            R3Rbar = R3Rbar + (R0_R(i)**3.d0)*weight(i)

                                            R3V2Lbar = R3V2Lbar + (R0_L(i)**3.d0)*(V0_L(i)**2.d0)*weight(i)
                                            R3V2Rbar = R3V2Rbar + (R0_R(i)**3.d0)*(V0_R(i)**2.d0)*weight(i)
                                        end do
                                    end if

                                    if (qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + num_fluids) < small_alf .or. R3Lbar < small_alf) then
                                        ptilde_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + num_fluids)*pres_L
                                    else
                                        ptilde_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + num_fluids)*(pres_L - PbwR3Lbar/R3Lbar - &
                                                                                                      rho_L*R3V2Lbar/R3Lbar)
                                    end if

                                    if (qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + num_fluids) < small_alf .or. R3Rbar < small_alf) then
                                        ptilde_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + num_fluids)*pres_R
                                    else
                                        ptilde_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + num_fluids)*(pres_R - PbwR3Rbar/R3Rbar - &
                                                                                                          rho_R*R3V2Rbar/R3Rbar)
                                    end if

                                    if ((ptilde_L /= ptilde_L) .or. (ptilde_R /= ptilde_R)) then
                                    end if

                                    rho_avg = 5d-1*(rho_L + rho_R)
                                    H_avg = 5d-1*(H_L + H_R)
                                    gamma_avg = 5d-1*(gamma_L + gamma_R)
                                    vel_avg_rms = 0d0

                                    !$acc loop seq
                                    do i = 1, num_dims
                                        vel_avg_rms = vel_avg_rms + (5d-1*(vel_L(i) + vel_R(i)))**2d0
                                    end do

                                end if

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                              vel_L_rms, c_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                              vel_R_rms, c_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.

                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                              vel_avg_rms, c_avg)

                                if (any(Re_size > 0)) then
                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
                                    end do
                                end if

                                if (wave_speeds == 1) then
                                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                                           (s_L - vel_L(dir_idx(1))) - &
                                           rho_R*vel_R(dir_idx(1))* &
                                           (s_R - vel_R(dir_idx(1)))) &
                                          /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                                            rho_R*(s_R - vel_R(dir_idx(1))))
                                elseif (wave_speeds == 2) then
                                    pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
                                                    (vel_L(dir_idx(1)) - &
                                                     vel_R(dir_idx(1))))

                                    pres_SR = pres_SL

                                    Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
                                                         (pres_SL/pres_L - 1d0)*pres_L/ &
                                                         ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
                                    Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
                                                         (pres_SR/pres_R - 1d0)*pres_R/ &
                                                         ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

                                    s_L = vel_L(dir_idx(1)) - c_L*Ms_L
                                    s_R = vel_R(dir_idx(1)) + c_R*Ms_R

                                    s_S = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                                                (pres_L - pres_R)/ &
                                                (rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M = min(0d0, s_L); s_P = max(0d0, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                                xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                                ! goes with numerical velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5d-1 + sign(5d-1, s_S))
                                xi_P = (5d-1 - sign(5d-1, s_S))

                                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                end do

                                if (bubbles .and. (num_fluids > 1)) then
                                    ! Kill mass transport @ gas density
                                    flux_rs${XYZ}$_vf(j, k, l, contxe) = 0.d0
                                end if

                                ! Momentum flux.
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)

                                ! Include p_tilde

                                !$acc loop seq
                                do i = 1, num_dims
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(i)) = &
                                        xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                     vel_L(dir_idx(i)) + &
                                                     s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + &
                                                                (1d0 - dir_flg(dir_idx(i)))* &
                                                                vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) + &
                                              dir_flg(dir_idx(i))*(pres_L - ptilde_L)) &
                                        + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                       vel_R(dir_idx(i)) + &
                                                       s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + &
                                                                  (1d0 - dir_flg(dir_idx(i)))* &
                                                                  vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) + &
                                                dir_flg(dir_idx(i))*(pres_R - ptilde_R))
                                    ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                                end do

                                ! Energy flux.
                                ! f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))

                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    xi_M*(vel_L(dir_idx(1))*(E_L + pres_L - ptilde_L) + &
                                          s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx(1)))* &
                                                     (rho_L*s_S + (pres_L - ptilde_L)/ &
                                                      (s_L - vel_L(dir_idx(1))))) - E_L)) &
                                    + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R - ptilde_R) + &
                                            s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx(1)))* &
                                                       (rho_R*s_S + (pres_R - ptilde_R)/ &
                                                        (s_R - vel_R(dir_idx(1))))) - E_R))

                                ! Volume fraction flux

                                !$acc loop seq
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                end do

                                ! Source for volume fraction advection equation
                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(i)) = &
                                        xi_M*(vel_L(dir_idx(i)) + &
                                              dir_flg(dir_idx(i))* &
                                              s_M*(xi_L - 1d0)) &
                                        + xi_P*(vel_R(dir_idx(i)) + &
                                                dir_flg(dir_idx(i))* &
                                                s_P*(xi_R - 1d0))

                                    !IF ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0d0
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, dir_idx(1))

                                ! Add advection flux for bubble variables
                                !$acc loop seq
                                do i = bubxb, bubxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*nbub_L*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                        + xi_P*nbub_R*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                end do

                                if (qbmm) then
                                    flux_rs${XYZ}$_vf(j, k, l, bubxb) = &
                                        xi_M*nbub_L &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                        + xi_P*nbub_R &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                end if

                                if (adv_n) then
                                    flux_rs${XYZ}$_vf(j, k, l, n_idx) = &
                                        xi_M*nbub_L &
                                        *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                                        + xi_P*nbub_R &
                                        *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                                end if

                                ! Geometrical source flux for cylindrical coordinates
                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        ! Substituting the advective flux into the inviscid geometrical source flux
                                        !$acc loop seq
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + dir_idx(1)) = &
                                            xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                         vel_L(dir_idx(1)) + &
                                                         s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + &
                                                                    (1d0 - dir_flg(dir_idx(1)))* &
                                                                    vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                            + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                           vel_R(dir_idx(1)) + &
                                                           s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + &
                                                                      (1d0 - dir_flg(dir_idx(1)))* &
                                                                      vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        !$acc loop seq
                                        do i = advxb, advxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0d0
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        !$acc loop seq
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0d0
                                        end do

                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb + 1) = &
                                            -xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                                          vel_L(dir_idx(1)) + &
                                                          s_M*(xi_L*(dir_flg(dir_idx(1))*s_S + &
                                                                     (1d0 - dir_flg(dir_idx(1)))* &
                                                                     vel_L(dir_idx(1))) - vel_L(dir_idx(1))))) &
                                            - xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                                           vel_R(dir_idx(1)) + &
                                                           s_P*(xi_R*(dir_flg(dir_idx(1))*s_S + &
                                                                      (1d0 - dir_flg(dir_idx(1)))* &
                                                                      vel_R(dir_idx(1))) - vel_R(dir_idx(1)))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxe) = flux_rs${XYZ}$_vf(j, k, l, momxb + 1)

                                    end if
                                #:endif
                            end do
                        end do
                    end do
                    !$acc end parallel loop
                else
                    !$acc parallel loop collapse(3) gang vector default(present) private(vel_L, vel_R, Re_L, Re_R, &
                    !$acc rho_avg, h_avg, gamma_avg, alpha_L, alpha_R, s_L, s_R, s_S, vel_avg_rms) copyin(is1,is2,is3)
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end

                                idx1 = 1; if (dir_idx(1) == 2) idx1 = 2; if (dir_idx(1) == 3) idx1 = 3

                                !$acc loop seq
                                do i = 1, num_fluids
                                    alpha_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    alpha_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                end do

                                vel_L_rms = 0d0; vel_R_rms = 0d0
                                !$acc loop seq
                                do i = 1, num_dims
                                    vel_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, contxe + i)
                                    vel_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, contxe + i)
                                    vel_L_rms = vel_L_rms + vel_L(i)**2d0
                                    vel_R_rms = vel_R_rms + vel_R(i)**2d0
                                end do

                                pres_L = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx)
                                pres_R = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx)

                                rho_L = 0d0
                                gamma_L = 0d0
                                pi_inf_L = 0d0
                                qv_L = 0d0

                                rho_R = 0d0
                                gamma_R = 0d0
                                pi_inf_R = 0d0
                                qv_R = 0d0

                                alpha_L_sum = 0d0
                                alpha_R_sum = 0d0

                                if (mpp_lim) then
                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, i) = max(0d0, qL_prim_rs${XYZ}$_vf(j, k, l, i))
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = min(max(0d0, qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)), 1d0)
                                        alpha_L_sum = alpha_L_sum + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)/max(alpha_L_sum, sgm_eps)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) = max(0d0, qR_prim_rs${XYZ}$_vf(j + 1, k, l, i))
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = min(max(0d0, qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)), 1d0)
                                        alpha_R_sum = alpha_R_sum + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)
                                    end do

                                    !$acc loop seq
                                    do i = 1, num_fluids
                                        qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)/max(alpha_R_sum, sgm_eps)
                                    end do
                                end if

                                !$acc loop seq
                                do i = 1, num_fluids
                                    rho_L = rho_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)
                                    gamma_L = gamma_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*gammas(i)
                                    pi_inf_L = pi_inf_L + qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + i)*pi_infs(i)
                                    qv_L = qv_L + qL_prim_rs${XYZ}$_vf(j, k, l, i)*qvs(i)

                                    rho_R = rho_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)
                                    gamma_R = gamma_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*gammas(i)
                                    pi_inf_R = pi_inf_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + i)*pi_infs(i)
                                    qv_R = qv_R + qR_prim_rs${XYZ}$_vf(j + 1, k, l, i)*qvs(i)
                                end do

                                if (any(Re_size > 0)) then
                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_L(i) = dflt_real

                                        if (Re_size(i) > 0) Re_L(i) = 0d0

                                        !$acc loop seq
                                        do q = 1, Re_size(i)
                                            Re_L(i) = qL_prim_rs${XYZ}$_vf(j, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_L(i)
                                        end do

                                        Re_L(i) = 1d0/max(Re_L(i), sgm_eps)

                                    end do

                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_R(i) = dflt_real

                                        if (Re_size(i) > 0) Re_R(i) = 0d0

                                        !$acc loop seq
                                        do q = 1, Re_size(i)
                                            Re_R(i) = qR_prim_rs${XYZ}$_vf(j + 1, k, l, E_idx + Re_idx(i, q))/Res(i, q) &
                                                      + Re_R(i)
                                        end do

                                        Re_R(i) = 1d0/max(Re_R(i), sgm_eps)
                                    end do
                                end if

                                E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*vel_L_rms + qv_L

                                E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*vel_R_rms + qv_R

                                H_L = (E_L + pres_L)/rho_L
                                H_R = (E_R + pres_R)/rho_R

                                @:compute_average_state()

                                call s_compute_speed_of_sound(pres_L, rho_L, gamma_L, pi_inf_L, H_L, alpha_L, &
                                                              vel_L_rms, c_L)

                                call s_compute_speed_of_sound(pres_R, rho_R, gamma_R, pi_inf_R, H_R, alpha_R, &
                                                              vel_R_rms, c_R)

                                !> The computation of c_avg does not require all the variables, and therefore the non '_avg'
                                ! variables are placeholders to call the subroutine.

                                call s_compute_speed_of_sound(pres_R, rho_avg, gamma_avg, pi_inf_R, H_avg, alpha_R, &
                                                              vel_avg_rms, c_avg)

                                if (any(Re_size > 0)) then
                                    !$acc loop seq
                                    do i = 1, 2
                                        Re_avg_rs${XYZ}$_vf(j, k, l, i) = 2d0/(1d0/Re_L(i) + 1d0/Re_R(i))
                                    end do
                                end if

                                if (wave_speeds == 1) then
                                    s_L = min(vel_L(idx1) - c_L, vel_R(idx1) - c_R)
                                    s_R = max(vel_R(idx1) + c_R, vel_L(idx1) + c_L)

                                    s_S = (pres_R - pres_L + rho_L*vel_L(idx1)* &
                                           (s_L - vel_L(idx1)) - &
                                           rho_R*vel_R(idx1)* &
                                           (s_R - vel_R(idx1))) &
                                          /(rho_L*(s_L - vel_L(idx1)) - &
                                            rho_R*(s_R - vel_R(idx1)))

                                elseif (wave_speeds == 2) then
                                    pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
                                                    (vel_L(idx1) - &
                                                     vel_R(idx1)))

                                    pres_SR = pres_SL

                                    Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
                                                         (pres_SL/pres_L - 1d0)*pres_L/ &
                                                         ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
                                    Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
                                                         (pres_SR/pres_R - 1d0)*pres_R/ &
                                                         ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

                                    s_L = vel_L(idx1) - c_L*Ms_L
                                    s_R = vel_R(idx1) + c_R*Ms_R

                                    s_S = 5d-1*((vel_L(idx1) + vel_R(idx1)) + &
                                                (pres_L - pres_R)/ &
                                                (rho_avg*c_avg))
                                end if

                                ! follows Einfeldt et al.
                                ! s_M/P = min/max(0.,s_L/R)
                                s_M = min(0d0, s_L); s_P = max(0d0, s_R)

                                ! goes with q_star_L/R = xi_L/R * (variable)
                                ! xi_L/R = ( ( s_L/R - u_L/R )/(s_L/R - s_star) )
                                xi_L = (s_L - vel_L(idx1))/(s_L - s_S)
                                xi_R = (s_R - vel_R(idx1))/(s_R - s_S)

                                ! goes with numerical velocity in x/y/z directions
                                ! xi_P/M = 0.5 +/m sgn(0.5,s_star)
                                xi_M = (5d-1 + sign(5d-1, s_S))
                                xi_P = (5d-1 - sign(5d-1, s_S))

                                !$acc loop seq
                                do i = 1, contxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1d0)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1d0))
                                end do

                                ! Momentum flux.
                                ! f = \rho u u + p I, q = \rho u, q_star = \xi * \rho*(s_star, v, w)
                                !$acc loop seq
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    flux_rs${XYZ}$_vf(j, k, l, contxe + idxi) = &
                                        xi_M*(rho_L*(vel_L(idx1)* &
                                                     vel_L(idxi) + &
                                                     s_M*(xi_L*(dir_flg(idxi)*s_S + &
                                                                (1d0 - dir_flg(idxi))* &
                                                                vel_L(idxi)) - vel_L(idxi))) + &
                                              dir_flg(idxi)*(pres_L)) &
                                        + xi_P*(rho_R*(vel_R(idx1)* &
                                                       vel_R(idxi) + &
                                                       s_P*(xi_R*(dir_flg(idxi)*s_S + &
                                                                  (1d0 - dir_flg(idxi))* &
                                                                  vel_R(idxi)) - vel_R(idxi))) + &
                                                dir_flg(idxi)*(pres_R))
                                    ! if (j==0) print*, 'flux_rs_vf', flux_rs_vf(cont_idx%end+dir_idx(i))%sf(j,k,l)
                                end do

                                ! Energy flux.
                                ! f = u*(E+p), q = E, q_star = \xi*E+(s-u)(\rho s_star + p/(s-u))
                                flux_rs${XYZ}$_vf(j, k, l, E_idx) = &
                                    xi_M*(vel_L(idx1)*(E_L + pres_L) + &
                                          s_M*(xi_L*(E_L + (s_S - vel_L(idx1))* &
                                                     (rho_L*s_S + pres_L/ &
                                                      (s_L - vel_L(idx1)))) - E_L)) &
                                    + xi_P*(vel_R(idx1)*(E_R + pres_R) + &
                                            s_P*(xi_R*(E_R + (s_S - vel_R(idx1))* &
                                                       (rho_R*s_S + pres_R/ &
                                                        (s_R - vel_R(idx1)))) - E_R))

                                ! Volume fraction flux
                                !$acc loop seq
                                do i = advxb, advxe
                                    flux_rs${XYZ}$_vf(j, k, l, i) = &
                                        xi_M*qL_prim_rs${XYZ}$_vf(j, k, l, i) &
                                        *(vel_L(idx1) + s_M*(xi_L - 1d0)) &
                                        + xi_P*qR_prim_rs${XYZ}$_vf(j + 1, k, l, i) &
                                        *(vel_R(idx1) + s_P*(xi_R - 1d0))
                                end do

                                ! Source for volume fraction advection equation
                                !$acc loop seq
                                do i = 1, num_dims
                                    idxi = dir_idx(i)
                                    vel_src_rs${XYZ}$_vf(j, k, l, idxi) = &
                                        xi_M*(vel_L(idxi) + &
                                              dir_flg(idxi)* &
                                              s_M*(xi_L - 1d0)) &
                                        + xi_P*(vel_R(idxi) + &
                                                dir_flg(idxi)* &
                                                s_P*(xi_R - 1d0))

                                    !if ( (model_eqns == 4) .or. (num_fluids==1) ) vel_src_rs_vf(dir_idx(i))%sf(j,k,l) = 0d0
                                end do

                                flux_src_rs${XYZ}$_vf(j, k, l, advxb) = vel_src_rs${XYZ}$_vf(j, k, l, idx1)

                                ! Geometrical source flux for cylindrical coordinates

                                #:if (NORM_DIR == 2)
                                    if (cyl_coord) then
                                        !Substituting the advective flux into the inviscid geometrical source flux
                                        !$acc loop seq
                                        do i = 1, E_idx
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = flux_rs${XYZ}$_vf(j, k, l, i)
                                        end do
                                        ! Recalculating the radial momentum geometric source flux
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, contxe + idx1) = &
                                            xi_M*(rho_L*(vel_L(idx1)* &
                                                         vel_L(idx1) + &
                                                         s_M*(xi_L*(dir_flg(idx1)*s_S + &
                                                                    (1d0 - dir_flg(idx1))* &
                                                                    vel_L(idx1)) - vel_L(idx1)))) &
                                            + xi_P*(rho_R*(vel_R(idx1)* &
                                                           vel_R(idx1) + &
                                                           s_P*(xi_R*(dir_flg(idx1)*s_S + &
                                                                      (1d0 - dir_flg(idx1))* &
                                                                      vel_R(idx1)) - vel_R(idx1))))
                                        ! Geometrical source of the void fraction(s) is zero
                                        !$acc loop seq
                                        do i = advxb, advxe
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0d0
                                        end do
                                    end if
                                #:endif
                                #:if (NORM_DIR == 3)
                                    if (grid_geometry == 3) then
                                        !$acc loop seq
                                        do i = 1, sys_size
                                            flux_gsrc_rs${XYZ}$_vf(j, k, l, i) = 0d0
                                        end do

                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxb + 1) = &
                                            -xi_M*(rho_L*(vel_L(idx1)* &
                                                          vel_L(idx1) + &
                                                          s_M*(xi_L*(dir_flg(idx1)*s_S + &
                                                                     (1d0 - dir_flg(idx1))* &
                                                                     vel_L(idx1)) - vel_L(idx1)))) &
                                            - xi_P*(rho_R*(vel_R(idx1)* &
                                                           vel_R(idx1) + &
                                                           s_P*(xi_R*(dir_flg(idx1)*s_S + &
                                                                      (1d0 - dir_flg(idx1))* &
                                                                      vel_R(idx1)) - vel_R(idx1))))
                                        flux_gsrc_rs${XYZ}$_vf(j, k, l, momxe) = flux_rs${XYZ}$_vf(j, k, l, momxb + 1)

                                    end if
                                #:endif

                            end do
                        end do
                    end do
                end if
            end if
        #:endfor
        ! Computing HLLC flux and source flux for Euler system of equations

        if (any(Re_size > 0)) then
            if (weno_Re_flux) then
                call s_compute_viscous_source_flux( &
                    qL_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    dqL_prim_dz_vf(momxb:momxe), &
                    qR_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    dqR_prim_dz_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy, iz)
            else
                call s_compute_viscous_source_flux( &
                    q_prim_vf(momxb:momxe), &
                    dqL_prim_dx_vf(momxb:momxe), &
                    dqL_prim_dy_vf(momxb:momxe), &
                    dqL_prim_dz_vf(momxb:momxe), &
                    q_prim_vf(momxb:momxe), &
                    dqR_prim_dx_vf(momxb:momxe), &
                    dqR_prim_dy_vf(momxb:momxe), &
                    dqR_prim_dz_vf(momxb:momxe), &
                    flux_src_vf, norm_dir, ix, iy, iz)
            end if
        end if

        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz, &
                                       inter_fluc_m_vf, inter_fluc_p_vf, &
                                       cell_fluc_vf, &
                                       phasic_pres_src_vf, &
                                       wave_propagation)

    end subroutine s_hllc_riemann_solver



    !>  This subroutine is utilized to compute u_hat in Path-conservative HLLC
        !!   Riemann solver.
    subroutine s_compute_PS_u_hat(u_hat, order, m, rhou_p, rhou_m, rho_p, rho_m) ! ---

        real(kind(0.0d0)), intent(inout) :: u_hat
        real(kind(0.0d0)), intent(in)    :: rhou_p, rhou_m, rho_p, rho_m
        integer, intent(in)              :: order, m

        ! local variables
        real(kind(0.0d0)) :: tolerance
        integer :: max_pts, init_pts
        real(kind(0.0d0)), allocatable :: x(:), w(:), s_nodes(:), weights(:)
        real(kind(0.0d0)) :: prev_u_hat, current_u_hat, integral, s, numerator, denominator
        integer :: n, i, iter_points, n_used
        logical :: converged = .false.

        ! Setting the tolerance and maximum number of points
        tolerance = 1d-6
        max_pts = 100
        init_pts = 5

        prev_u_hat = 0d0
        current_u_hat = 0d0
        n_used = init_pts
        converged = .false.
        do while (n_used <= max_pts .and. .not. converged)
            ! Setting the number of points
            n = n_used

            ! Allocating the nodes and weights
            allocate(x(n), w(n), s_nodes(n), weights(n))

            ! Setting the nodes and weights
            call s_compute_gauss_legendre_nodes_and_weights(n, x, w)

            ! Setting the nodes and weights for the s variable
            s_nodes = 0.5d0*(x + 1d0)
            weights = 0.5d0*w

            ! Computing the integral
            integral = 0d0
            do i = 1, n
                s = s_nodes(i)
                numerator = (rhou_m + s**m*(rhou_p - rhou_m))
                denominator = (rho_m + s**m*(rho_p - rho_m))
                integral = integral + weights(i)*(numerator/denominator)**order*s**(m-1)
            end do

            ! Computing the current u_hat
            current_u_hat = m*integral

            ! Checking for convergence
            if (abs(current_u_hat - prev_u_hat) < tolerance) then
                converged = .true.
            else
                prev_u_hat = current_u_hat
                n_used = n_used + 5
            end if

            ! Deallocating the nodes and weights
            deallocate(x, w, s_nodes, weights)

            u_hat = current_u_hat

        end do

    end subroutine s_compute_PS_u_hat

    subroutine s_compute_gauss_legendre_nodes_and_weights(n, x, w) ! ----------------

        integer, intent(in) :: n
        real(kind(0.0d0)), intent(out) :: x(n), w(n)

        ! local variables
        real(kind(0.0d0))   ::  pi = 4.0d0*atan(1.0d0)
        integer             ::  i, j
        real(kind(0.0d0))   ::  z, pp, p_pre

        do i = 1, n
            z = cos(pi*(i - 0.25d0)/(n + 0.5d0))

            do j = 1, 100
                call s_compute_legendre_polynomial(n, z, pp)
                if (abs(pp) < 1d-15) exit
                call s_compute_legendre_polynomial(n - 1, z, p_pre)
                z = z - pp/(n*(p_pre - z*pp)/(1.0d0 - z**2))
            end do
            x(n - i + 1) = z
            call s_compute_legendre_polynomial(n - 1, z, p_pre)
            w(n - i + 1) = 2d0*(1d0 - z**2)/(n**2*p_pre**2)
        end do

    end subroutine s_compute_gauss_legendre_nodes_and_weights

    subroutine s_compute_legendre_polynomial(n, x, p_now) ! ----------------

        integer, intent(in) :: n
        real(kind(0.0d0)), intent(in) :: x
        real(kind(0.0d0)), intent(out) :: p_now

        ! local variables
        real(kind(0.0d0)) :: p0, p1, p2
        integer :: i

        p0 = 1d0
        p1 = x
        if (n == 0) then
            p_now = p0
        else if (n == 1) then
            p_now = p1
        else
            do i = 2, n
                p2 = ((2d0*i - 1d0)*x*p1 - (i - 1d0)*p0)/i
                p0 = p1
                p1 = p2
            end do
            p_now = p2
        end if

    end subroutine s_compute_legendre_polynomial

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_riemann_solvers_module() ! ---------------------

        ! Allocating the variables that will be utilized to formulate the
        ! left, right, and average states of the Riemann problem, as well
        ! the Riemann problem solution
        integer :: i, j
        INTEGER :: xbeg,xend,ybeg,yend,zbeg,zend

        ! Setting up special bounds for cell-average values
        xbeg = -buff_size; ybeg = 0; zbeg = 0
        IF (n > 0) ybeg = -buff_size; IF (p > 0) zbeg = -buff_size
        xend = m - xbeg; yend = n - ybeg; zend = p - zbeg

        @:ALLOCATE_GLOBAL(Gs(1:num_fluids))

        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        !$acc update device(Gs)

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

        !$acc enter data copyin(is1, is2, is3, isx, isy, isz)

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution of a given Riemann problem
        if (riemann_solver == 1) then
            s_riemann_solver => s_hll_riemann_solver
        elseif (riemann_solver == 2) then
            s_riemann_solver => s_hllc_riemann_solver
        end if

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to compute the viscous source flux
        if (grid_geometry == 3) then
            s_compute_viscous_source_flux => s_compute_cylindrical_viscous_source_flux
        else
            s_compute_viscous_source_flux => s_compute_cartesian_viscous_source_flux
        end if

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = m; is2%end = n; is3%end = p

        @:ALLOCATE_GLOBAL(flux_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(flux_gsrc_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(flux_src_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, advxb:sys_size))
        @:ALLOCATE_GLOBAL(vel_src_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:num_dims))
        @:ALLOCATE_GLOBAL(cell_fluc_rsx(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(inter_fluc_rsx_m(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(inter_fluc_rsx_p(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(pha_pres_src_rsx_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, intxb:sys_size))
        @:ALLOCATE_GLOBAL(q_prim_rsx_vf(xbeg:xend, &
            ybeg:yend, &
            zbeg:zend, 1:sys_size))
        if (qbmm) then
            @:ALLOCATE_GLOBAL(mom_sp_rsx_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (any(Re_size > 0)) then
            @:ALLOCATE_GLOBAL(Re_avg_rsx_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:2))
        end if

        if (n == 0) return

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = n; is2%end = m; is3%end = p

        @:ALLOCATE_GLOBAL(flux_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(flux_gsrc_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(flux_src_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, advxb:sys_size))
        @:ALLOCATE_GLOBAL(vel_src_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:num_dims))
        @:ALLOCATE_GLOBAL(cell_fluc_rsy(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(inter_fluc_rsy_m(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(inter_fluc_rsy_p(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(pha_pres_src_rsy_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, intxb:sys_size))
        @:ALLOCATE_GLOBAL(q_prim_rsy_vf(ybeg:yend, &
            xbeg:xend, &
            zbeg:zend, 1:sys_size))

        if (qbmm) then
            @:ALLOCATE_GLOBAL(mom_sp_rsy_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (any(Re_size > 0)) then
            @:ALLOCATE_GLOBAL(Re_avg_rsy_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:2))
        end if

        if (p == 0) return

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = p; is2%end = n; is3%end = m

        @:ALLOCATE_GLOBAL(flux_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(flux_gsrc_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(flux_src_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, advxb:sys_size))
        @:ALLOCATE_GLOBAL(vel_src_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:num_dims))
        @:ALLOCATE_GLOBAL(cell_fluc_rsz(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(inter_fluc_rsz_m(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(inter_fluc_rsz_p(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))
        @:ALLOCATE_GLOBAL(pha_pres_src_rsz_vf(is1%beg:is1%end, &
            is2%beg:is2%end, &
            is3%beg:is3%end, intxb:sys_size))
        @:ALLOCATE_GLOBAL(q_prim_rsz_vf(zbeg:zend, &
            ybeg:yend, &
            xbeg:xend, 1:sys_size))

        if (qbmm) then
            @:ALLOCATE_GLOBAL(mom_sp_rsz_vf(is1%beg:is1%end + 1, is2%beg:is2%end, is3%beg:is3%end, 1:4))
        end if

        if (any(Re_size > 0)) then
            @:ALLOCATE_GLOBAL(Re_avg_rsz_vf(is1%beg:is1%end, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:2))
        end if

    end subroutine s_initialize_riemann_solvers_module ! -------------------

    !>  The purpose of this subroutine is to populate the buffers
        !!      of the left and right Riemann states variables, depending
        !!      on the boundary conditions.
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param dqL_prim_dx_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqL_prim_dy_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqL_prim_dz_vf The  left WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order x-dir spatial derivatives
        !!  @param dqR_prim_dy_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order y-dir spatial derivatives
        !!  @param dqR_prim_dz_vf The right WENO-reconstructed cell-boundary values of the
        !!      first-order z-dir spatial derivatives
        !!  @param gm_alphaL_vf  Left averaged gradient magnitude
        !!  @param gm_alphaR_vf Right averaged gradient magnitude
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
    subroutine s_populate_riemann_states_variables_buffers( & ! ------------
        qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, dqL_prim_dx_vf, &
        dqL_prim_dy_vf, &
        dqL_prim_dz_vf, &
        qL_prim_vf, &
        qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf, dqR_prim_dx_vf, &
        dqR_prim_dy_vf, &
        dqR_prim_dz_vf, &
        qR_prim_vf, &
        norm_dir, ix, iy, iz)

        real(kind(0d0)), dimension(startx:, starty:, startz:, 1:), intent(INOUT) :: qL_prim_rsx_vf, qL_prim_rsy_vf, qL_prim_rsz_vf, qR_prim_rsx_vf, qR_prim_rsy_vf, qR_prim_rsz_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             dqL_prim_dy_vf, dqR_prim_dy_vf, &
                             dqL_prim_dz_vf, dqR_prim_dz_vf, &
                             qL_prim_vf, qR_prim_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k, l !< Generic loop iterator

        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1d0, 0d0, 0d0/)
        elseif (norm_dir == 2) then
            is1 = iy; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0d0, 1d0, 0d0/)
        else
            is1 = iz; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0d0, 0d0, 1d0/)
        end if

        !$acc update device(is1, is2, is3)

        if (hypoelasticity) then
            if (norm_dir == 1) then
                dir_idx_tau = (/1, 2, 4/)
            else if (norm_dir == 2) then
                dir_idx_tau = (/3, 2, 5/)
            else
                dir_idx_tau = (/6, 4, 5/)
            end if
        end if

        isx = ix; isy = iy; isz = iz
        !$acc update device(isx, isy, isz) ! for stuff in the same module
        !$acc update device(dir_idx, dir_flg,  dir_idx_tau) ! for stuff in different modules

        ! Population of Buffers in x-direction =============================
        if (norm_dir == 1) then

            if (bc_x%beg == -4) then    ! Riemann state extrap. BC at beginning
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsx_vf(-1, k, l, i) = &
                                qR_prim_rsx_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end

                                dqL_prim_dx_vf(i)%sf(-1, k, l) = &
                                    dqR_prim_dx_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do

                    if (n > 0) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end

                                    dqL_prim_dy_vf(i)%sf(-1, k, l) = &
                                        dqR_prim_dy_vf(i)%sf(0, k, l)
                                end do
                            end do
                        end do

                        if (p > 0) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = momxb, momxe
                                do l = isz%beg, isz%end
                                    do k = isy%beg, isy%end

                                        dqL_prim_dz_vf(i)%sf(-1, k, l) = &
                                            dqR_prim_dz_vf(i)%sf(0, k, l)
                                    end do
                                end do
                            end do
                        end if

                    end if

                end if

            end if

            if (bc_x%end == -4) then    ! Riemann state extrap. BC at end

                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsx_vf(m + 1, k, l, i) = &
                                qL_prim_rsx_vf(m, k, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do k = isy%beg, isy%end

                                dqR_prim_dx_vf(i)%sf(m + 1, k, l) = &
                                    dqL_prim_dx_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do

                    if (n > 0) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do k = isy%beg, isy%end

                                    dqR_prim_dy_vf(i)%sf(m + 1, k, l) = &
                                        dqL_prim_dy_vf(i)%sf(m, k, l)
                                end do
                            end do
                        end do

                        if (p > 0) then
                            !$acc parallel loop collapse(3) gang vector default(present)
                            do i = momxb, momxe
                                do l = isz%beg, isz%end
                                    do k = isy%beg, isy%end

                                        dqR_prim_dz_vf(i)%sf(m + 1, k, l) = &
                                            dqL_prim_dz_vf(i)%sf(m, k, l)
                                    end do
                                end do
                            end do
                        end if

                    end if

                end if

            end if
            ! END: Population of Buffers in x-direction ========================

            ! Population of Buffers in y-direction =============================
        elseif (norm_dir == 2) then

            if (bc_y%beg == -4) then    ! Riemann state extrap. BC at beginning
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsy_vf(-1, k, l, i) = &
                                qR_prim_rsy_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, -1, l) = &
                                    dqR_prim_dx_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, -1, l) = &
                                    dqR_prim_dy_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do

                    if (p > 0) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do j = isx%beg, isx%end
                                    dqL_prim_dz_vf(i)%sf(j, -1, l) = &
                                        dqR_prim_dz_vf(i)%sf(j, 0, l)
                                end do
                            end do
                        end do
                    end if

                end if

            end if

            if (bc_y%end == -4) then    ! Riemann state extrap. BC at end

                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsy_vf(n + 1, k, l, i) = &
                                qL_prim_rsy_vf(n, k, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, n + 1, l) = &
                                    dqL_prim_dx_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do l = isz%beg, isz%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, n + 1, l) = &
                                    dqL_prim_dy_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do

                    if (p > 0) then
                        !$acc parallel loop collapse(3) gang vector default(present)
                        do i = momxb, momxe
                            do l = isz%beg, isz%end
                                do j = isx%beg, isx%end
                                    dqR_prim_dz_vf(i)%sf(j, n + 1, l) = &
                                        dqL_prim_dz_vf(i)%sf(j, n, l)
                                end do
                            end do
                        end do
                    end if

                end if

            end if
            ! END: Population of Buffers in y-direction ========================

            ! Population of Buffers in z-direction =============================
        else

            if (bc_z%beg == -4) then    ! Riemann state extrap. BC at beginning
                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qL_prim_rsz_vf(-1, k, l, i) = &
                                qR_prim_rsz_vf(0, k, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dx_vf(i)%sf(j, k, -1) = &
                                    dqR_prim_dx_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dy_vf(i)%sf(j, k, -1) = &
                                    dqR_prim_dy_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqL_prim_dz_vf(i)%sf(j, k, -1) = &
                                    dqR_prim_dz_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                end if

            end if

            if (bc_z%end == -4) then    ! Riemann state extrap. BC at end

                !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            qR_prim_rsz_vf(p + 1, k, l, i) = &
                                qL_prim_rsz_vf(p, k, l, i)
                        end do
                    end do
                end do

                if (any(Re_size > 0)) then
                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dx_vf(i)%sf(j, k, p + 1) = &
                                    dqL_prim_dx_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dy_vf(i)%sf(j, k, p + 1) = &
                                    dqL_prim_dy_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do

                    !$acc parallel loop collapse(3) gang vector default(present)
                    do i = momxb, momxe
                        do k = isy%beg, isy%end
                            do j = isx%beg, isx%end
                                dqR_prim_dz_vf(i)%sf(j, k, p + 1) = &
                                    dqL_prim_dz_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                end if

            end if

        end if
        ! END: Population of Buffers in z-direction ========================

    end subroutine s_populate_riemann_states_variables_buffers ! -----------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures needed to configure the chosen Riemann
        !!      solver algorithm.
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_initialize_riemann_solver( &
        q_prim_vf, &
        flux_vf, flux_src_vf, &
        flux_gsrc_vf, &
        norm_dir, ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k, l ! Generic loop iterators
        INTEGER :: xbeg,xend,ybeg,yend,zbeg,zend
        INTEGER :: s1beg,s1end,s2beg,s2end,s3beg,s3end

        ! Setting up special bounds for cell-average values
        xbeg = -buff_size; ybeg = 0; zbeg = 0
        IF (n > 0) ybeg = -buff_size; IF (p > 0) zbeg = -buff_size
        xend = m - xbeg; yend = n - ybeg; zend = p - zbeg

        ! Configuring the coordinate direction indexes
        IF (norm_dir == 1) THEN
            s1beg = xbeg; s1end = xend; s2beg = ybeg; s2end = yend; s3beg = zbeg; s3end = zend
        ELSEIF (norm_dir == 2) THEN
            s1beg = ybeg; s1end = yend; s2beg = xbeg; s2end = xend; s3beg = zbeg; s3end = zend
        ELSE
            s1beg = zbeg; s1end = zend; s2beg = ybeg; s2end = yend; s3beg = xbeg; s3end = xend
        END IF

        ! Reshaping Inputted Data in x-direction ===========================

        if (norm_dir == 1) then

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = s3beg, s3end
                    do k = s2beg, s2end
                        do j = s1beg, s1end
                            q_prim_rsx_vf(j, k, l, i) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do

            if (any(Re_size > 0)) then

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = momxb, E_idx
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k, l) = 0d0
                            end do
                        end do
                    end do
                end do
            end if

            if (qbmm) then

                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsx_vf(j, k, l, i) = mom_sp(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
            end if

            ! ==================================================================

            ! Reshaping Inputted Data in y-direction ===========================
        elseif (norm_dir == 2) then

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = s3beg, s3end
                    do k = s2beg, s2end
                        do j = s1beg, s1end
                            q_prim_rsy_vf(j, k, l, i) = q_prim_vf(i)%sf(k, j, l)
                        end do
                    end do
                end do
            end do

            if (any(Re_size > 0)) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = momxb, E_idx
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j, l) = 0d0
                            end do
                        end do
                    end do
                end do
            end if

            if (qbmm) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsy_vf(j, k, l, i) = mom_sp(i)%sf(k, j, l)
                            end do
                        end do
                    end do
                end do
            end if

            ! ==================================================================

            ! Reshaping Inputted Data in z-direction ===========================
        else

            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = s3beg, s3end
                    do k = s2beg, s2end
                        do j = s1beg, s1end
                            q_prim_rsz_vf(j, k, l, i) = q_prim_vf(i)%sf(l, k, j)
                        end do
                    end do
                end do
            end do

            if (any(Re_size > 0)) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = momxb, E_idx
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_src_vf(i)%sf(l, k, j) = 0d0
                            end do
                        end do
                    end do
                end do
            end if

            if (qbmm) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, 4
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end + 1
                                mom_sp_rsz_vf(j, k, l, i) = mom_sp(i)%sf(l, k, j)
                            end do
                        end do
                    end do
                end do
            end if

        end if

        ! ==================================================================

    end subroutine s_initialize_riemann_solver ! ---------------------------

    !>  The goal of this subroutine is to evaluate and account
        !!      for the contribution of viscous stresses in the source
        !!      flux for the momentum and energy.
        !!  @param velL_vf  Left, WENO reconstructed, cell-boundary values of the velocity
        !!  @param velR_vf Right, WENO reconstructed, cell-boundary values of the velocity
        !!  @param dvelL_dx_vf  Left, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelL_dy_vf  Left, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelL_dz_vf  Left, WENO reconstructed cell-avg. z-dir derivative of the velocity
        !!  @param dvelR_dx_vf Right, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelR_dy_vf Right, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelR_dz_vf Right, WENO reconstructed cell-avg. z-dir derivative of the velocity
        !!  @param flux_src_vf Intercell flux
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix Index bounds in  first coordinate direction
        !!  @param iy Index bounds in second coordinate direction
        !!  @param iz Index bounds in  third coordinate direction
    subroutine s_compute_cylindrical_viscous_source_flux(velL_vf, & ! -------------
                                                         dvelL_dx_vf, &
                                                         dvelL_dy_vf, &
                                                         dvelL_dz_vf, &
                                                         velR_vf, &
                                                         dvelR_dx_vf, &
                                                         dvelR_dy_vf, &
                                                         dvelR_dz_vf, &
                                                         flux_src_vf, &
                                                         norm_dir, &
                                                         ix, iy, iz)

        type(scalar_field), &
            dimension(num_dims), &
            intent(IN) :: velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_src_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy, iz

        ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
        ! boundary values of cell-average first-order spatial derivatives
        ! of velocity
        real(kind(0d0)), dimension(num_dims) :: avg_vel
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dx
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dy
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dz

        ! Viscous stress tensor
        real(kind(0d0)), dimension(num_dims, num_dims) :: tau_Re

        ! Generic loop iterators
        integer :: i, j, k, l

        ! Viscous Stresses in z-direction ==================================
        if (norm_dir == 1) then
            if (Re_size(1) > 0) then ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (4d0/3d0)*dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (n == 0) return

            if (Re_size(1) > 0) then ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j + 1, k, l))

                            !$acc loop seq
                            do i = 1, 2
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(2) = 5d-1*(dvelL_dx_vf(2)%sf(j, k, l) &
                                                   + dvelR_dx_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2d0/3d0)*(dvel_avg_dy(2) + &
                                                       avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            tau_Re(1, 2) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            !$acc loop seq
                            do i = 1, 2
                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(1, i)
                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsx_vf(j, k, l, i)* &
                                    tau_Re(1, i)
                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel,  dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j + 1, k, l))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (dvel_avg_dy(2) + &
                                            avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (Re_size(1) > 0) then ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 3, 2
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(3) = 5d-1*(dvelL_dx_vf(3)%sf(j, k, l) &
                                                   + dvelR_dx_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2d0/3d0)*dvel_avg_dz(3)/y_cc(k)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            tau_Re(1, 3) = (dvel_avg_dz(1)/y_cc(k) + dvel_avg_dx(3))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            !$acc loop seq
                            do i = 1, 3, 2

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsx_vf(j, k, l, i)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( avg_vel, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dz(3)/y_cc(k)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in z-direction =============================

            ! Viscous Stresses in r-direction ==================================
        elseif (norm_dir == 2) then

            if (Re_size(1) > 0) then ! Shear stresses

                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j, k + 1, l))

                            !$acc loop seq
                            do i = 1, 2

                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k, l) &
                                          + dvelR_dx_vf(i)%sf(j, k + 1, l))

                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j, k + 1, l))

                            end do

                            tau_Re(2, 1) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            tau_Re(2, 2) = (4d0*dvel_avg_dy(2) &
                                            - 2d0*dvel_avg_dx(1) &
                                            - 2d0*avg_vel(2)/y_cb(k))/ &
                                           (3d0*Re_avg_rsy_vf(k, j, l, 1))

                            !$acc loop seq
                            do i = 1, 2

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsy_vf(k, j, l, i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j, k + 1, l))

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j, k + 1, l))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j, k + 1, l))

                            tau_Re(2, 2) = (dvel_avg_dx(1) + dvel_avg_dy(2) + &
                                            avg_vel(2)/y_cb(k))/ &
                                           Re_avg_rsy_vf(k, j, l, 2)

                            flux_src_vf(momxb + 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsy_vf(k, j, l, 2)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel,  dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(3) = 5d-1*(velL_vf(3)%sf(j, k, l) &
                                               + velR_vf(3)%sf(j, k + 1, l))

                            !$acc loop seq
                            do i = 2, 3
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j, k + 1, l))
                            end do

                            dvel_avg_dy(3) = 5d-1*(dvelL_dy_vf(3)%sf(j, k, l) &
                                                   + dvelR_dy_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = -(2d0/3d0)*dvel_avg_dz(3)/y_cb(k)/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            tau_Re(2, 3) = ((dvel_avg_dz(2) - avg_vel(3))/ &
                                            y_cb(k) + dvel_avg_dy(3))/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            !$acc loop seq
                            do i = 2, 3

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsy_vf(k, j, l, i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel,  dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = dvel_avg_dz(3)/y_cb(k)/ &
                                           Re_avg_rsy_vf(k, j, l, 2)

                            flux_src_vf(momxb + 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsy_vf(k, j, l, 2)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in r-direction =============================

            ! Viscous Stresses in theta-direction ==================================
        else

            if (Re_size(1) > 0) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 2, 3
                                avg_vel(i) = 5d-1*(velL_vf(i)%sf(j, k, l) &
                                                   + velR_vf(i)%sf(j, k, l + 1))
                            end do

                            !$acc loop seq
                            do i = 1, 3, 2
                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k, l) &
                                          + dvelR_dx_vf(i)%sf(j, k, l + 1))
                            end do

                            do i = 2, 3
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j, k, l + 1))
                            end do

                            !$acc loop seq
                            do i = 1, 3
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j, k, l + 1))
                            end do

                            tau_Re(3, 1) = (dvel_avg_dz(1)/y_cc(k) + dvel_avg_dx(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 1)/ &
                                           y_cc(k)

                            tau_Re(3, 2) = ((dvel_avg_dz(2) - avg_vel(3))/ &
                                            y_cc(k) + dvel_avg_dy(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 1)/ &
                                           y_cc(k)

                            tau_Re(3, 3) = (4d0*dvel_avg_dz(3)/y_cc(k) &
                                            - 2d0*dvel_avg_dx(1) &
                                            - 2d0*dvel_avg_dy(2) &
                                            + 4d0*avg_vel(2)/y_cc(k))/ &
                                           (3d0*Re_avg_rsz_vf(l, k, j, 1))/ &
                                           y_cc(k)

                            !$acc loop seq
                            do i = 1, 3
                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(3, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsz_vf(l, k, j, i)* &
                                    tau_Re(3, i)
                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(avg_vel, dvel_avg_dx, dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            avg_vel(2) = 5d-1*(velL_vf(2)%sf(j, k, l) &
                                               + velR_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j, k, l + 1))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j, k, l + 1))

                            tau_Re(3, 3) = (dvel_avg_dx(1) &
                                            + dvel_avg_dy(2) &
                                            + dvel_avg_dz(3)/y_cc(k) &
                                            + avg_vel(2)/y_cc(k))/ &
                                           Re_avg_rsz_vf(l, k, j, 2)/ &
                                           y_cc(k)

                            flux_src_vf(momxe)%sf(j, k, l) = &
                                flux_src_vf(momxe)%sf(j, k, l) - &
                                tau_Re(3, 3)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsz_vf(l, k, j, 3)* &
                                tau_Re(3, 3)

                        end do
                    end do
                end do
            end if

        end if
        ! END: Viscous Stresses in theta-direction =============================

    end subroutine s_compute_cylindrical_viscous_source_flux ! -------------------------

    !>  The goal of this subroutine is to evaluate and account
        !!      for the contribution of viscous stresses in the source
        !!      flux for the momentum and energy.
        !!  @param velL_vf  Left, WENO reconstructed, cell-boundary values of the velocity
        !!  @param velR_vf Right, WENO reconstructed, cell-boundary values of the velocity
        !!  @param dvelL_dx_vf  Left, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelL_dy_vf  Left, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelL_dz_vf  Left, WENO reconstructed cell-avg. z-dir derivative of the velocity
        !!  @param dvelR_dx_vf Right, WENO reconstructed cell-avg. x-dir derivative of the velocity
        !!  @param dvelR_dy_vf Right, WENO reconstructed cell-avg. y-dir derivative of the velocity
        !!  @param dvelR_dz_vf Right, WENO reconstructed cell-avg. z-dir derivative of the velocity
        !!  @param flux_src_vf Intercell flux
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix Index bounds in  first coordinate direction
        !!  @param iy Index bounds in second coordinate direction
        !!  @param iz Index bounds in  third coordinate direction
    subroutine s_compute_cartesian_viscous_source_flux(velL_vf, & ! -------------
                                                       dvelL_dx_vf, &
                                                       dvelL_dy_vf, &
                                                       dvelL_dz_vf, &
                                                       velR_vf, &
                                                       dvelR_dx_vf, &
                                                       dvelR_dy_vf, &
                                                       dvelR_dz_vf, &
                                                       flux_src_vf, &
                                                       norm_dir, &
                                                       ix, iy, iz)

        type(scalar_field), &
            dimension(num_dims), &
            intent(IN) :: velL_vf, velR_vf, &
                          dvelL_dx_vf, dvelR_dx_vf, &
                          dvelL_dy_vf, dvelR_dy_vf, &
                          dvelL_dz_vf, dvelR_dz_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_src_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy, iz

        ! Arithmetic mean of the left and right, WENO-reconstructed, cell-
        ! boundary values of cell-average first-order spatial derivatives
        ! of velocity
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dx
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dy
        real(kind(0d0)), dimension(num_dims) :: dvel_avg_dz

        real(kind(0d0)), dimension(num_dims, num_dims) :: tau_Re !< Viscous stress tensor

        integer :: i, j, k, l !< Generic loop iterators

        ! Viscous Stresses in x-direction ==================================
        if (norm_dir == 1) then

            if (Re_size(1) > 0) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = (4d0/3d0)*dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dx(1)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (n == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 2
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(2) = 5d-1*(dvelL_dx_vf(2)%sf(j, k, l) &
                                                   + dvelR_dx_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2d0/3d0)*dvel_avg_dy(2)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            tau_Re(1, 2) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            !$acc loop seq
                            do i = 1, 2

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsx_vf(j, k, l, i)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dy(2)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 3, 2
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j + 1, k, l))
                            end do

                            dvel_avg_dx(3) = 5d-1*(dvelL_dx_vf(3)%sf(j, k, l) &
                                                   + dvelR_dx_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = -(2d0/3d0)*dvel_avg_dz(3)/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            tau_Re(1, 3) = (dvel_avg_dz(1) + dvel_avg_dx(3))/ &
                                           Re_avg_rsx_vf(j, k, l, 1)

                            !$acc loop seq
                            do i = 1, 3, 2
                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(1, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsx_vf(j, k, l, i)* &
                                    tau_Re(1, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j + 1, k, l))

                            tau_Re(1, 1) = dvel_avg_dz(3)/ &
                                           Re_avg_rsx_vf(j, k, l, 2)

                            flux_src_vf(momxb)%sf(j, k, l) = &
                                flux_src_vf(momxb)%sf(j, k, l) - &
                                tau_Re(1, 1)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsx_vf(j, k, l, 1)* &
                                tau_Re(1, 1)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in x-direction =============================

            ! Viscous Stresses in y-direction ==================================
        elseif (norm_dir == 2) then

            if (Re_size(1) > 0) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 2

                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k, l) &
                                          + dvelR_dx_vf(i)%sf(j, k + 1, l))

                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j, k + 1, l))

                            end do

                            tau_Re(2, 1) = (dvel_avg_dy(1) + dvel_avg_dx(2))/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            tau_Re(2, 2) = (4d0*dvel_avg_dy(2) &
                                            - 2d0*dvel_avg_dx(1))/ &
                                           (3d0*Re_avg_rsy_vf(k, j, l, 1))

                            !$acc loop seq
                            do i = 1, 2

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsy_vf(k, j, l, i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j, k + 1, l))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j, k + 1, l))

                            tau_Re(2, 2) = (dvel_avg_dx(1) + dvel_avg_dy(2))/ &
                                           Re_avg_rsy_vf(k, j, l, 2)

                            flux_src_vf(momxb + 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsy_vf(k, j, l, 2)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if

            if (p == 0) return

            if (Re_size(1) > 0) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private(  dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 2, 3
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j, k + 1, l))
                            end do

                            dvel_avg_dy(3) = 5d-1*(dvelL_dy_vf(3)%sf(j, k, l) &
                                                   + dvelR_dy_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = -(2d0/3d0)*dvel_avg_dz(3)/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            tau_Re(2, 3) = (dvel_avg_dz(2) + dvel_avg_dy(3))/ &
                                           Re_avg_rsy_vf(k, j, l, 1)

                            !$acc loop seq
                            do i = 2, 3

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(2, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsy_vf(k, j, l, i)* &
                                    tau_Re(2, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j, k + 1, l))

                            tau_Re(2, 2) = dvel_avg_dz(3)/ &
                                           Re_avg_rsy_vf(k, j, l, 2)

                            flux_src_vf(momxb + 1)%sf(j, k, l) = &
                                flux_src_vf(momxb + 1)%sf(j, k, l) - &
                                tau_Re(2, 2)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsy_vf(k, j, l, 2)* &
                                tau_Re(2, 2)

                        end do
                    end do
                end do
            end if
            ! END: Viscous Stresses in y-direction =============================

            ! Viscous Stresses in z-direction ==================================
        else

            if (Re_size(1) > 0) then              ! Shear stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            !$acc loop seq
                            do i = 1, 3, 2
                                dvel_avg_dx(i) = &
                                    5d-1*(dvelL_dx_vf(i)%sf(j, k, l) &
                                          + dvelR_dx_vf(i)%sf(j, k, l + 1))
                            end do

                            !$acc loop seq
                            do i = 2, 3
                                dvel_avg_dy(i) = &
                                    5d-1*(dvelL_dy_vf(i)%sf(j, k, l) &
                                          + dvelR_dy_vf(i)%sf(j, k, l + 1))
                            end do

                            !$acc loop seq
                            do i = 1, 3
                                dvel_avg_dz(i) = &
                                    5d-1*(dvelL_dz_vf(i)%sf(j, k, l) &
                                          + dvelR_dz_vf(i)%sf(j, k, l + 1))
                            end do

                            tau_Re(3, 1) = (dvel_avg_dz(1) + dvel_avg_dx(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 1)

                            tau_Re(3, 2) = (dvel_avg_dz(2) + dvel_avg_dy(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 1)

                            tau_Re(3, 3) = (4d0*dvel_avg_dz(3) &
                                            - 2d0*dvel_avg_dx(1) &
                                            - 2d0*dvel_avg_dy(2))/ &
                                           (3d0*Re_avg_rsz_vf(l, k, j, 1))

                            !$acc loop seq
                            do i = 1, 3

                                flux_src_vf(contxe + i)%sf(j, k, l) = &
                                    flux_src_vf(contxe + i)%sf(j, k, l) - &
                                    tau_Re(3, i)

                                flux_src_vf(E_idx)%sf(j, k, l) = &
                                    flux_src_vf(E_idx)%sf(j, k, l) - &
                                    vel_src_rsz_vf(l, k, j, i)* &
                                    tau_Re(3, i)

                            end do

                        end do
                    end do
                end do
            end if

            if (Re_size(2) > 0) then              ! Bulk stresses
                !$acc parallel loop collapse(3) gang vector default(present) private( dvel_avg_dx, dvel_avg_dy, dvel_avg_dz, tau_Re)
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end

                            dvel_avg_dx(1) = 5d-1*(dvelL_dx_vf(1)%sf(j, k, l) &
                                                   + dvelR_dx_vf(1)%sf(j, k, l + 1))

                            dvel_avg_dy(2) = 5d-1*(dvelL_dy_vf(2)%sf(j, k, l) &
                                                   + dvelR_dy_vf(2)%sf(j, k, l + 1))

                            dvel_avg_dz(3) = 5d-1*(dvelL_dz_vf(3)%sf(j, k, l) &
                                                   + dvelR_dz_vf(3)%sf(j, k, l + 1))

                            tau_Re(3, 3) = (dvel_avg_dx(1) &
                                            + dvel_avg_dy(2) &
                                            + dvel_avg_dz(3))/ &
                                           Re_avg_rsz_vf(l, k, j, 2)

                            flux_src_vf(momxe)%sf(j, k, l) = &
                                flux_src_vf(momxe)%sf(j, k, l) - &
                                tau_Re(3, 3)

                            flux_src_vf(E_idx)%sf(j, k, l) = &
                                flux_src_vf(E_idx)%sf(j, k, l) - &
                                vel_src_rsz_vf(l, k, j, 3)* &
                                tau_Re(3, 3)

                        end do
                    end do
                end do
            end if

        end if
        ! END: Viscous Stresses in z-direction =============================

    end subroutine s_compute_cartesian_viscous_source_flux ! -------------------------

    !>  Deallocation and/or disassociation procedures that are
        !!      needed to finalize the selected Riemann problem solver
        !!  @param flux_vf       Intercell fluxes
        !!  @param flux_src_vf   Intercell source fluxes
        !!  @param flux_gsrc_vf  Intercell geometric source fluxes
        !!  @param norm_dir Dimensional splitting coordinate direction
        !!  @param ix   Index bounds in  first coordinate direction
        !!  @param iy   Index bounds in second coordinate direction
        !!  @param iz   Index bounds in  third coordinate direction
    subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, & ! --------
                                         flux_gsrc_vf, &
                                         norm_dir, ix, iy, iz, &
                                         inter_fluc_m_vf, inter_fluc_p_vf, &
                                         cell_fluc_vf, &
                                         phasic_pres_src_vf, &
                                         wave_propagation)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf, inter_fluc_m_vf, inter_fluc_p_vf, cell_fluc_vf, phasic_pres_src_vf

        integer, intent(IN) :: norm_dir

        type(int_bounds_info), intent(IN) :: ix, iy, iz
        logical, intent(IN) :: wave_propagation

        integer :: i, j, k, l !< Generic loop iterators

        ! Reshaping Outputted Data in y-direction ==========================
        if (norm_dir == 2) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            if (.not. wave_propagation) then
                                flux_vf(i)%sf(k, j, l) = &
                                    flux_rsy_vf(j, k, l, i)
                            else
                                cell_fluc_vf(i)%sf(k, j, l) = &
                                    cell_fluc_rsy(j, k, l, i)
                                inter_fluc_m_vf(i)%sf(k, j, l) = &
                                    inter_fluc_rsy_m(j, k, l, i)
                                inter_fluc_p_vf(i)%sf(k, j, l) = &
                                    inter_fluc_rsy_p(j, k, l, i)
                            end if
                            ! assign phasic pressure source terms in phasic totoal energy equation
                            if (i >= intxb) then
                                phasic_pres_src_vf(i)%sf(k, j, l) = &
                                    pha_pres_src_rsy_vf(j, k, l, i)
                            end if
                        end do
                    end do
                end do
            end do

            if (cyl_coord) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_gsrc_vf(i)%sf(k, j, l) = &
                                    flux_gsrc_rsy_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
            end if

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3%beg, is3%end
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        flux_src_vf(advxb)%sf(k, j, l) = &
                            flux_src_rsy_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                    do l = is3%beg, is3%end
                        do j = is1%beg, is1%end
                            do k = is2%beg, is2%end
                                flux_src_vf(i)%sf(k, j, l) = &
                                    flux_src_rsy_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do

            end if
            ! ==================================================================
            ! Reshaping Outputted Data in z-direction ==========================
        elseif (norm_dir == 3) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
                            if (.not. wave_propagation) then
                                flux_vf(i)%sf(l, k, j) = &
                                    flux_rsz_vf(j, k, l, i)
                            else
                                cell_fluc_vf(i)%sf(l, k, j) = &
                                    cell_fluc_rsz(j, k, l, i)
                                inter_fluc_m_vf(i)%sf(l, k, j) = &
                                    inter_fluc_rsz_m(j, k, l, i)
                                inter_fluc_p_vf(i)%sf(l, k, j) = &
                                    inter_fluc_rsz_p(j, k, l, i)
                            end if
                            ! assign phasic pressure source terms in phasic totoal energy equation
                            if (i >= intxb) then
                                phasic_pres_src_vf(i)%sf(l, k, j) = &
                                    pha_pres_src_rsz_vf(j, k, l, i)
                            end if
                        end do
                    end do
                end do
            end do
            if (grid_geometry == 3) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, sys_size
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end

                                flux_gsrc_vf(i)%sf(l, k, j) = &
                                    flux_gsrc_rsz_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
            end if

            !$acc parallel loop collapse(3) gang vector default(present)
            do j = is1%beg, is1%end
                do k = is2%beg, is2%end
                    do l = is3%beg, is3%end
                        flux_src_vf(advxb)%sf(l, k, j) = &
                            flux_src_rsz_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            do l = is3%beg, is3%end
                                flux_src_vf(i)%sf(l, k, j) = &
                                    flux_src_rsz_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do

            end if
        elseif (norm_dir == 1) then
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            if (.not. wave_propagation) then
                                flux_vf(i)%sf(j, k, l) = &
                                    flux_rsx_vf(j, k, l, i)
                            else
                                cell_fluc_vf(i)%sf(j, k, l) = &
                                    cell_fluc_rsx(j, k, l, i)
                                inter_fluc_m_vf(i)%sf(j, k, l) = &
                                    inter_fluc_rsx_m(j, k, l, i)
                                inter_fluc_p_vf(i)%sf(j, k, l) = &
                                    inter_fluc_rsx_p(j, k, l, i)
                            end if
                            ! assign phasic pressure source terms in phasic totoal energy equation
                            if (i >= intxb) then
                                phasic_pres_src_vf(i)%sf(j, k, l) = &
                                    pha_pres_src_rsx_vf(j, k, l, i)
                            end if
                        end do
                    end do
                end do
            end do

            !$acc parallel loop collapse(3) gang vector default(present)
            do l = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = is1%beg, is1%end
                        flux_src_vf(advxb)%sf(j, k, l) = &
                            flux_src_rsx_vf(j, k, l, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                !$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                    do l = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = is1%beg, is1%end
                                flux_src_vf(i)%sf(j, k, l) = &
                                    flux_src_rsx_vf(j, k, l, i)
                            end do
                        end do
                    end do
                end do
            end if
        end if

        ! ==================================================================

        ! ==================================================================

    end subroutine s_finalize_riemann_solver ! -----------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_riemann_solvers_module() ! -----------------------

        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the solution of a given Riemann problem
        s_riemann_solver => null()

        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the viscous source flux
        s_compute_viscous_source_flux => null()

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        ! s_convert_to_mixture_variables => null()

        if (Re_size(1) > 0) then
            @:DEALLOCATE_GLOBAL(Re_avg_rsx_vf)
        end if
        @:DEALLOCATE_GLOBAL(vel_src_rsx_vf)
        @:DEALLOCATE_GLOBAL(cell_fluc_rsx)
        @:DEALLOCATE_GLOBAL(inter_fluc_rsx_m)
        @:DEALLOCATE_GLOBAL(inter_fluc_rsx_p)
        @:DEALLOCATE_GLOBAL(pha_pres_src_rsx_vf)
        @:DEALLOCATE_GLOBAL(q_prim_rsx_vf)
        @:DEALLOCATE_GLOBAL(flux_rsx_vf)
        @:DEALLOCATE_GLOBAL(flux_src_rsx_vf)
        @:DEALLOCATE_GLOBAL(flux_gsrc_rsx_vf)
        if (qbmm) then
            @:DEALLOCATE_GLOBAL(mom_sp_rsx_vf)
        end if

        if (n == 0) return

        if (Re_size(1) > 0) then
            @:DEALLOCATE_GLOBAL(Re_avg_rsy_vf)
        end if
        @:DEALLOCATE_GLOBAL(vel_src_rsy_vf)
        @:DEALLOCATE_GLOBAL(cell_fluc_rsy)
        @:DEALLOCATE_GLOBAL(inter_fluc_rsy_m)
        @:DEALLOCATE_GLOBAL(inter_fluc_rsy_p)
        @:DEALLOCATE_GLOBAL(pha_pres_src_rsy_vf)
        @:DEALLOCATE_GLOBAL(q_prim_rsy_vf)
        @:DEALLOCATE_GLOBAL(flux_rsy_vf)
        @:DEALLOCATE_GLOBAL(flux_src_rsy_vf)
        @:DEALLOCATE_GLOBAL(flux_gsrc_rsy_vf)
        if (qbmm) then
            @:DEALLOCATE_GLOBAL(mom_sp_rsy_vf)
        end if

        if (p == 0) return

        if (Re_size(1) > 0) then
            @:DEALLOCATE_GLOBAL(Re_avg_rsz_vf)
        end if
        @:DEALLOCATE_GLOBAL(vel_src_rsz_vf)
        @:DEALLOCATE_GLOBAL(cell_fluc_rsz)
        @:DEALLOCATE_GLOBAL(inter_fluc_rsz_m)
        @:DEALLOCATE_GLOBAL(inter_fluc_rsz_p)
        @:DEALLOCATE_GLOBAL(pha_pres_src_rsz_vf)
        @:DEALLOCATE_GLOBAL(q_prim_rsz_vf)
        @:DEALLOCATE_GLOBAL(flux_rsz_vf)
        @:DEALLOCATE_GLOBAL(flux_src_rsz_vf)
        @:DEALLOCATE_GLOBAL(flux_gsrc_rsz_vf)
        if (qbmm) then
            @:DEALLOCATE_GLOBAL(mom_sp_rsz_vf)
        end if

    end subroutine s_finalize_riemann_solvers_module ! ---------------------

end module m_riemann_solvers
