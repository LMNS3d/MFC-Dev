!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers

#:include 'macros.fpp'

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles              !< Bubble dynamics routines

    use m_ibm

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_boundary_conditions

    use m_helper

    use m_fftw

    use m_nvtx

    use ieee_arithmetic
    ! ==========================================================================

    implicit none

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), q_cons_ts)
    !! Cell-average conservative variables at each time-stage (TS)

    @:CRAY_DECLARE_GLOBAL(type(scalar_field), dimension(:), q_prim_vf)
    !! Cell-average primitive variables at the current time-stage

    @:CRAY_DECLARE_GLOBAL(type(scalar_field), dimension(:), rhs_vf)
    !! Cell-average RHS variables at the current time-stage

    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), lo_rhs_vf)
    !! Cell-average RHS variables with low-order scheme at the current time-stage

    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), q_prim_ts)
    !! Cell-average primitive variables at consecutive TIMESTEPS

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :, :), rhs_pb)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :, :, :), rhs_mv)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), MOOD_indicator)

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

    !$acc declare link(q_cons_ts,q_prim_vf,rhs_vf,lo_rhs_vf,q_prim_ts, rhs_mv, rhs_pb, MOOD_indicator)
#else
    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: lo_rhs_vf !<
    !! Cell-average RHS variables with low-order scheme at the current time-stage

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: rhs_pb

    real(kind(0d0)), allocatable, dimension(:, :, :, :, :) :: rhs_mv

    real(kind(0d0)), allocatable, dimension(:, :, :) :: MOOD_indicator

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

    !$acc declare create(q_cons_ts,q_prim_vf,rhs_vf,lo_rhs_vf,q_prim_ts, rhs_mv, rhs_pb,MOOD_indicator)
#endif

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module() ! -----------------------

        type(int_bounds_info) :: ix_t, iy_t, iz_t !<
            !! Indical bounds in the x-, y- and z-directions

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3/))) then
            num_ts = 2
        end if

        ! Setting the indical bounds in the x-, y- and z-directions
        ix_t%beg = -buff_size; ix_t%end = m + buff_size

        if (n > 0) then
            iy_t%beg = -buff_size; iy_t%end = n + buff_size

            if (p > 0) then
                iz_t%beg = -buff_size; iz_t%end = p + buff_size
            else
                iz_t%beg = 0; iz_t%end = 0
            end if
        else
            iy_t%beg = 0; iy_t%end = 0
            iz_t%beg = 0; iz_t%end = 0
        end if

        ! Allocating the cell-average conservative variables
        @:ALLOCATE_GLOBAL(q_cons_ts(1:num_ts))

        do i = 1, num_ts
            @:ALLOCATE(q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
                @:ALLOCATE(q_cons_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
            end do
            @:ACC_SETUP_VFs(q_cons_ts(i))
        end do

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
            @:ALLOCATE_GLOBAL(q_prim_ts(0:3))

            do i = 0, 3
                @:ALLOCATE(q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
                    @:ALLOCATE(q_prim_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, &
                        iy_t%beg:iy_t%end, &
                        iz_t%beg:iz_t%end))
                end do
            end do

            do i = 0, 3
                @:ACC_SETUP_VFs(q_prim_ts(i))
            end do
        end if

        ! Allocating the cell-average primitive variables
        @:ALLOCATE_GLOBAL(q_prim_vf(1:sys_size))

        do i = 1, adv_idx%end
            @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end))
            @:ACC_SETUP_SFs(q_prim_vf(i))
        end do

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
            if (adv_n) then
                @:ALLOCATE(q_prim_vf(n_idx)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(n_idx))
            end if
        end if

        if (hypoelasticity) then

            do i = stress_idx%beg, stress_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                    iy_t%beg:iy_t%end, &
                    iz_t%beg:iz_t%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do
        end if

        @:ALLOCATE_GLOBAL(pb_ts(1:2))
        !Initialize bubble variables pb and mv at all quadrature nodes for all R0 bins
        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(pb_ts(1)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE_GLOBAL(rhs_pb(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(pb_ts(1)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE_GLOBAL(rhs_pb(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
        end if

        @:ALLOCATE_GLOBAL(mv_ts(1:2))

        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(mv_ts(1)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE_GLOBAL(rhs_mv(ix_t%beg:ix_t%end, &
                iy_t%beg:iy_t%end, &
                iz_t%beg:iz_t%end, 1:nnode, 1:nb))

        else if (qbmm .and. polytropic) then
            @:ALLOCATE(mv_ts(1)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE_GLOBAL(rhs_mv(ix_t%beg:ix_t%beg + 1, &
                iy_t%beg:iy_t%beg + 1, &
                iz_t%beg:iz_t%beg + 1, 1:nnode, 1:nb))
        end if

        ! Allocating the cell-average RHS variables
        @:ALLOCATE_GLOBAL(rhs_vf(1:sys_size))

        do i = 1, sys_size
            @:ALLOCATE(rhs_vf(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(rhs_vf(i))
        end do

        ! Allocating the cell-average RHS variables with low-order scheme
        @:ALLOCATE_GLOBAL(lo_rhs_vf(1:sys_size))
        do i = 1, sys_size
            @:ALLOCATE(lo_rhs_vf(i)%sf(0:m, 0:n, 0:p))
            @:ACC_SETUP_SFs(lo_rhs_vf(i))
        end do

        ! Allocating the MOOD indicator variables
        @:ALLOCATE_GLOBAL(MOOD_indicator(-3:m+3, -3:n+3, -3:p+3))

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

    end subroutine s_initialize_time_steppers_module ! ---------------------

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: start, finish
        real(kind(0d0)) :: nR3bar
        logical :: MOOD_button

        MOOD_button = .true.

        ! Stage 1 of 1 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, &
                                t_step, MOOD_button, lo_rhs_vf)

#ifdef DEBUG
        print *, 'got rhs'
#endif

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

#ifdef DEBUG
        print *, 'wrote runtime info'
#endif

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_Wakimura2024(q_cons_ts(1)%vf)

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_1st_order_tvd_rk ! ------------------------------------

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l, q!< Generic loop iterator
        real(kind(0d0)) :: start, finish
        real(kind(0d0)) :: nR3bar
        logical :: MOOD_button

        MOOD_button = .true.

        ! Stage 1 of 2 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, &
                                t_step, MOOD_button, lo_rhs_vf)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3) then
            call s_pressure_relaxation_Wakimura2024(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if
        ! ==================================================================

        ! Stage 2 of 2 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, &
                                t_step, MOOD_button, lo_rhs_vf)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/2d0
                    end do
                end do
            end do
        end do

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_pb(j, k, l, q, i))/2d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt*rhs_mv(j, k, l, q, i))/2d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) then
            call s_pressure_relaxation_Wakimura2024(q_cons_ts(1)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_2nd_order_tvd_rk ! ------------------------------------

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step, time_avg, dt_in) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg
        real(kind(0d0)), intent(IN) :: dt_in

        integer :: i, j, k, l, q !< Generic loop iterator
        real(kind(0d0)) :: ts_error, denom, error_fraction, time_step_factor !< Generic loop iterator
        real(kind(0d0)) :: start, finish
        real(kind(0d0)) :: nR3bar
        logical :: MOOD_button
        real(kind(0d0)) :: MOOD_eps1, MOOD_eps2
        integer :: neighbors    !< Number of influenced neighbors of troubled cells
        real(kind(0d0)), dimension(1:sys_size) :: U_next !< Array of the next sub-step values

        CHARACTER(LEN = path_len + 2*name_len) :: t_step_dir !<
        !! Relative path to the current time-step directory


        CHARACTER(LEN = path_len + 3*name_len) :: file_path !<

        CHARACTER(LEN=15) :: FMT

        MOOD_button = .true.
        MOOD_eps1 = 1d-10
        MOOD_eps2 = 1d-8
        neighbors = 1

        ! Stage 1 of 3 =====================================================

        if (.not. adap_dt) then
            call cpu_time(start)
            call nvtxStartRange("Time_Step")
        end if

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, &
                                t_step, MOOD_button, lo_rhs_vf)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        ! variable for data restore
        FMT="(2F30.7)"
        WRITE(t_step_dir,'(A,I0,A,I0)') TRIM(case_dir) // '/D'

        ! NAD and PAD for MOOD in first sub-step ---------------------
        if (MOOD_button) then
            ! Initialize the MOOD indicator
            MOOD_indicator(:,:,:) = 0d0
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        ! Compute the next sub-step values with high-order scheme
                        do i = 1, sys_size
                            U_next(i) = &
                                q_cons_ts(1)%vf(i)%sf(j, k, l) &
                                + dt_in*rhs_vf(i)%sf(j, k, l)
                        end do
                        if (abs(rhs_vf(advxb)%sf(j, k, l) + rhs_vf(advxe)%sf(j, k, l))*dt_in > 1d-10) then
                            print *, 'volume fraction ', rhs_vf(advxb)%sf(j, k, l) + rhs_vf(advxe)%sf(j, k, l), &
                                        U_next(advxb) + U_next(advxe)
                        end if
                        ! Detect the troubled cells which violate the NAD and PAD criteria
                        call s_mood_nad_pad(j,k,l,U_next, MOOD_eps1, MOOD_eps2, neighbors)

                    end do
                end do
            end do

            ! do l = 0, p
            !     do k = 0, n
            !         do j = 0, m
            !             if (MOOD_indicator(j, k, l) > 0) then
            !                 ! print *, 'troubled cell in first sub-step: ', 'j =',j,'k =', k,'l =', l
            !                 i = advxb
            !                 U_next(i) = &
            !                     q_cons_ts(1)%vf(i)%sf(j, k, l) &
            !                     + dt_in*rhs_vf(i)%sf(j, k, l)
            !                 ! print *, 'rhs = ', rhs_vf(i)%sf(j, k, l)
            !                 ! print *, 'lo_rhs = ', lo_rhs_vf(i)%sf(j, k, l)
            !                 if (U_next(i) < MOOD_eps1 .or. U_next(i) > 1d0-MOOD_eps1) then
            !                     print *, 'PAD criteria violated by alpha '
                            
            !                     print *, 'U_now = ', q_cons_ts(1)%vf(i)%sf(j, k, l)
            !                     print *, 'U_original_next = ', U_next(i)
            !                     U_next(i) = &
            !                         q_cons_ts(1)%vf(i)%sf(j, k, l) &
            !                         + dt_in*(rhs_vf(i)%sf(j, k, l)*(1d0-MOOD_indicator(j, k, l)) &
            !                         + lo_rhs_vf(i)%sf(j, k, l)*MOOD_indicator(j, k, l))
            !                     print *, 'U_MOOD_next = ', U_next(i)
            !                 end if
            !             end if
            !         end do
            !     end do
            ! end do

            ! print *, sum(MOOD_indicator(0:m,0:n,0:p)), t_step+1, proc_rank

            IF(MOD(t_step+1-t_step_start, t_step_save) == 0) THEN

                IF (n==0 .AND. p==0) THEN

        WRITE(file_path,'(A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/mood1.', proc_rank, '.', t_step+1,'.dat'
                    open(2,file= trim(file_path))
                        do j=0,m
                            write(2,FMT) x_cb(j), MOOD_indicator(j,0,0)
                        end do
                    close(2)

                else

        WRITE(file_path,'(A,I2.2,A,I6.6,A)') TRIM(t_step_dir) // '/mood1.', proc_rank, '.', t_step+1,'.dat'
                    open(2,file= trim(file_path))
                        DO j=0,m
                        DO k=0,n
                            WRITE(2,FMT) x_cb(j),y_cb(k), MOOD_indicator(j,k,0)
                        END DO
                        WRITE(2,*)
                        END DO
                    close(2)
                endif

            END IF

            ! Compute the next sub-step values with the troubled cells recalculated by low-order scheme
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                            q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                            + dt_in*(rhs_vf(i)%sf*(1d0-MOOD_indicator(0:m,0:n,0:p)) &
                            + lo_rhs_vf(i)%sf*MOOD_indicator(0:m,0:n,0:p)) 
            END DO

            ! if (sum(MOOD_indicator(0:m,0:n,0:p))>0) then
            !     print *, 'total number of troubled cells in first sub-step = ', sum(MOOD_indicator(0:m,0:n,0:p))
            ! end if
            
            ! if (sum(MOOD_indicator(0:m,0:n,0:p))>3) then
            !     call s_mpi_abort('NAD and PAD criteria violated m')
            ! end if

        else
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                                q_cons_ts(1)%vf(i)%sf(j, k, l) &
                                + dt_in*rhs_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        end if
        ! end of NAD and PAD for MOOD in first sub-step ---------------------

        !Evolve pb and mv for non-polytropic qbmm
        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    pb_ts(1)%sf(j, k, l, q, i) &
                                    + dt_in*rhs_pb(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    mv_ts(1)%sf(j, k, l, q, i) &
                                    + dt_in*rhs_mv(j, k, l, q, i)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3) then
            call s_pressure_relaxation_Wakimura2024(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! ==================================================================

        ! Stage 2 of 3 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, &
                                t_step, MOOD_button, lo_rhs_vf)

        ! NAD and PAD for MOOD in second sub-step ---------------------
        if (MOOD_button) then
            ! Initialize the MOOD indicator
            MOOD_indicator(:,:,:) = 0d0
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        ! Compute the next sub-step values with high-order scheme
                        do i = 1, sys_size
                            U_next(i) = &
                                (3d0*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                                + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                                + dt_in*rhs_vf(i)%sf(j, k, l))/4d0
                        end do
                        ! Detect the troubled cells which debate the NAD and PAD criteria
                        call s_mood_nad_pad(j,k,l,U_next, MOOD_eps1, MOOD_eps2, neighbors)

                    end do
                end do
            end do

            ! Compute the next sub-step values with the troubled cells recalculated by low-order scheme
            DO i = 1, sys_size
                q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) = &
                                ( 3d0*q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                                +     q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) &
                                +     dt_in*(rhs_vf(i)%sf*(1d0-MOOD_indicator(0:m,0:n,0:p)) &
                                + lo_rhs_vf(i)%sf*MOOD_indicator(0:m,0:n,0:p)) ) / 4d0
            END DO
            ! if (sum(MOOD_indicator(0:m,0:n,0:p))>0) then
            !     print *, 'total number of troubled cells in second sub-step = ', sum(MOOD_indicator(0:m,0:n,0:p))
            ! end if

        else
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                                (3d0*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                                + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                                + dt_in*rhs_vf(i)%sf(j, k, l))/4d0
                        end do
                    end do
                end do
            end do
        end if
        ! end of NAD and PAD for MOOD in second sub-step ---------------------

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(2)%sf(j, k, l, q, i) = &
                                    (3d0*pb_ts(1)%sf(j, k, l, q, i) &
                                     + pb_ts(2)%sf(j, k, l, q, i) &
                                     + dt_in*rhs_pb(j, k, l, q, i))/4d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(2)%sf(j, k, l, q, i) = &
                                    (3d0*mv_ts(1)%sf(j, k, l, q, i) &
                                     + mv_ts(2)%sf(j, k, l, q, i) &
                                     + dt_in*rhs_mv(j, k, l, q, i))/4d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3) then
            call s_pressure_relaxation_Wakimura2024(q_cons_ts(2)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(2)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf, pb_ts(2)%sf, mv_ts(2)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(2)%vf, q_prim_vf)
            end if
        end if

        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, pb_ts(2)%sf, rhs_pb, mv_ts(2)%sf, rhs_mv, &
                                t_step, MOOD_button, lo_rhs_vf)

        ! NAD and PAD for MOOD in third sub-step ---------------------
        if (MOOD_button) then
            ! Initialize the MOOD indicator
            MOOD_indicator(:,:,:) = 0d0
            !$acc parallel loop collapse(3) gang vector default(present)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        ! Compute the next sub-step values with high-order scheme
                        do i = 1, sys_size
                            U_next(i) = &
                                (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                                + 2d0*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                                + 2d0*dt_in*rhs_vf(i)%sf(j, k, l))/3d0
                        end do
                        ! Detect the troubled cells which violate the NAD and PAD criteria
                        call s_mood_nad_pad(j,k,l,U_next, MOOD_eps1, MOOD_eps2, neighbors)

                    end do
                end do
            end do

            ! Compute the next sub-step values with the troubled cells recalculated by low-order scheme
            DO i = 1, sys_size
                q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) = &
                                (     q_cons_ts(1)%vf(i)%sf(0:m,0:n,0:p) &
                                + 2d0*q_cons_ts(2)%vf(i)%sf(0:m,0:n,0:p) &
                                + 2d0*dt_in*(rhs_vf(i)%sf*(1d0-MOOD_indicator(0:m,0:n,0:p)) &
                                + lo_rhs_vf(i)%sf*MOOD_indicator(0:m,0:n,0:p)) ) / 3d0
            END DO
            
            ! if (sum(MOOD_indicator(0:m,0:n,0:p))>0) then
            !     print *, 'total number of troubled cells in third sub-step = ', sum(MOOD_indicator(0:m,0:n,0:p))
            ! end if

        else
            !$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                                (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                                + 2d0*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                                + 2d0*dt_in*rhs_vf(i)%sf(j, k, l))/3d0
                        end do
                    end do
                end do
            end do
        end if
        ! end of NAD and PAD for MOOD in third sub-step ---------------------

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                pb_ts(1)%sf(j, k, l, q, i) = &
                                    (pb_ts(1)%sf(j, k, l, q, i) &
                                     + 2d0*pb_ts(2)%sf(j, k, l, q, i) &
                                     + 2d0*dt_in*rhs_pb(j, k, l, q, i))/3d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (qbmm .and. (.not. polytropic)) then
            !$acc parallel loop collapse(5) gang vector default(present)
            do i = 1, nb
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            do q = 1, nnode
                                mv_ts(1)%sf(j, k, l, q, i) = &
                                    (mv_ts(1)%sf(j, k, l, q, i) &
                                     + 2d0*mv_ts(2)%sf(j, k, l, q, i) &
                                     + 2d0*dt_in*rhs_mv(j, k, l, q, i))/3d0
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) then
            call s_pressure_relaxation_Wakimura2024(q_cons_ts(1)%vf)
        end if

        if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        if (ib) then
            if (qbmm .and. .not. polytropic) then
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            else
                call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
            end if
        end if

        if (.not. adap_dt) then
            call nvtxEndRange
            call cpu_time(finish)

            time = time + (finish - start)

            if (t_step >= 4) then
                time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
            else
                time_avg = 0d0
            end if
        end if
        ! ==================================================================

    end subroutine s_3rd_order_tvd_rk ! ------------------------------------

    !> This subroutine computes the next sub-step values with high-order scheme
        !!      and detects the troubled cells which violate the NAD and PAD criteria
        !! @param j,k,l Current cell index
        !! @param U_next Array of the next sub-step values with given rhs values
        !! @param MOOD_eps1,MOOD_eps2 Thresholds for the MOOD algorithm, 1 for volume fraction and phasic mass; 2 for pressure and sound speed
        !! @param neighbors Number of influenced neighbors of troubled cells
    subroutine s_mood_nad_pad(j,k,l,U_next, MOOD_eps1, MOOD_eps2, neighbors) ! ----------------

        integer, intent(IN) :: j, k, l
        real(kind(0d0)), intent(IN), dimension(1:sys_size) :: U_next
        real(kind(0d0)), intent(IN) :: MOOD_eps1, MOOD_eps2
        integer, intent(IN) :: neighbors

        integer :: i !< Generic loop iterator

        ! Variables for checking the MOOD criteria
        real(kind(0d0)), dimension(1:num_fluids) :: alpha_rho, alpha
        real(kind(0d0)) :: rho, rho_c2, total_Energy, Gamma_mix, Pi_mix, pressure, Qv_mix
        real(kind(0d0)) :: internal_Energy, dyn_pres
        real(kind(0d0)), dimension(1:num_dims) :: rho_u
        real(kind(0d0)), dimension(1:num_fluids) :: partial_total_Energy, partial_internal_energy
        real(kind(0d0)), dimension(1:num_fluids) :: rho_c2_k, pressure_k
        logical :: interface_indicator !< Interface indicator, pressure criterion is only active for interface cells

        ! NAD
        !$acc parallel loop gang vector default(present)
        do i = 1, sys_size
            if (ieee_is_nan(U_next(i)) &
                .or. (.not. ieee_is_finite(U_next(i)))) then

                rhs_vf(i)%sf(j,k,l) = 0d0

                ! print *, 'troubled cell caused by NAD: ', 'j =',j,'k =', k,'l =', l
                ! call s_mpi_abort('NAD criteria violated')

                MOOD_indicator(j-neighbors:j+neighbors, k, l) = 1d0
                MOOD_indicator(j, k-neighbors:k+neighbors, l) = 1d0
                MOOD_indicator(j, k, l-neighbors:l+neighbors) = 1d0

            endif
        enddo

        ! PAD
        !> Compute the pressure and sound speed and some other variables
        rho = 0d0; Gamma_mix = 0d0; Pi_mix = 0d0; Qv_mix = 0d0
        do i = 1, num_fluids
            alpha_rho(i) = U_next(i)
            alpha(i) = U_next(i + E_idx)
            rho = rho + alpha_rho(i)
            Gamma_mix = Gamma_mix + alpha(i)*fluid_pp(i)%gamma
            Pi_mix = Pi_mix + alpha(i)*fluid_pp(i)%pi_inf
            Qv_mix = Qv_mix + alpha_rho(i)*fluid_pp(i)%qv
        end do
        dyn_pres = 0d0
        do i = momxb, momxe
            dyn_pres = dyn_pres + 5d-1*U_next(i)* &
                        U_next(i)/max(rho, sgm_eps)
        end do
        ! calculate the phasic pressure and sound speed
        do i = 1, num_fluids
            partial_total_Energy(i) = U_next(i + adv_idx%end)
            partial_internal_energy(i) = partial_total_Energy(i) - &
                                    dyn_pres/max(rho, sgm_eps)*alpha_rho(i)
            pressure_k(i) = (partial_internal_energy(i) - alpha_rho(i)*fluid_pp(i)%qv &
                             - alpha(i)*fluid_pp(i)%pi_inf)/(alpha(i)*fluid_pp(i)%gamma)
            rho_c2_k(i) = ((fluid_pp(i)%gamma + 1d0)*pressure_k(i) + &
                            fluid_pp(i)%pi_inf)/(fluid_pp(i)%gamma)
        end do
        ! calculate the mixture pressure and sound speed
        total_Energy = U_next(E_idx)
        internal_Energy = total_Energy - dyn_pres
        pressure = (internal_Energy - Qv_mix - Pi_mix)/Gamma_mix
        rho_c2 = ((Gamma_mix+1d0)*pressure + Pi_mix)/Gamma_mix

        if (alpha(1) > 1d-6 .and. alpha(1) < 1d0 - 1d-6) then 
            interface_indicator = .TRUE.
        else
            interface_indicator = .FALSE.
        endif

        ! check the PAD MOOD criteria
        if (     alpha_rho(1) < MOOD_eps1 .or. alpha_rho(2) < MOOD_eps1 &
            .or. alpha(1) < MOOD_eps1 .or. alpha(1) > 1d0 - MOOD_eps1 &
            .or. alpha(2) < MOOD_eps1 .or. alpha(2) > 1d0 - MOOD_eps1 &
            ! .or. (internal_Energy < MOOD_eps2) &
            ! .or. (partial_internal_energy(1) < MOOD_eps2) &
            ! .or. (partial_internal_energy(2) < MOOD_eps2) &
            ! .or. (pressure < MOOD_eps2 .and. interface_indicator) &
            .or. (rho_c2 < MOOD_eps2) &
            .or. (rho_c2_k(1) < MOOD_eps2) &
            .or. (rho_c2_k(2) < MOOD_eps2) &
            ! .or. ieee_is_nan(rho_c2) .or. ieee_is_nan(internal_Energy) &
            ! .or. ieee_is_nan(pressure) &
            )THEN
            ! ! print * , 'alpha1, alpha2 ', alpha(1), alpha(2), sum(alpha)
            ! print *, 'troubled cell caused by PAD: ', 'j =',j,'k =', k,'l =', l
            ! print *, 'm =', m, 'n =', n, 'p =', p
            ! print *, 'alpha_rho1, alpha_rho2 ', alpha_rho(1), alpha_rho(2)
            ! print *, 'alpha1, alpha2 ', alpha(1), alpha(2)
            ! print *, 'rho_c2 ', rho_c2
            ! print *, 'internal_Energy ', internal_Energy
            ! print *, 'pressure ', pressure
            ! if (pressure < MOOD_eps2 .and. interface_indicator) then
            !     ! call s_mpi_abort('PAD criteria violated')
            !     print *, 'PAD criteria violated by pressure ', pressure
            ! end if
            ! if (alpha(1) < MOOD_eps1 .or. alpha(2) < MOOD_eps1 &
            !     .or. alpha(1) > 1d0 - MOOD_eps1 .or. alpha(2) > 1d0 - MOOD_eps1 &
            !     ) then
            !     ! call s_mpi_abort('NAD and PAD criteria violated')
            !     print *, 'PAD criteria violated by alpha ', alpha(1), alpha(2)
            ! end if
            ! if (internal_Energy < MOOD_eps2) then
            !     ! call s_mpi_abort('PAD criteria violated')
            !     print *, 'PAD criteria violated by internal_Energy ', internal_Energy
            ! end if
            ! if (partial_internal_Energy(1) < MOOD_eps2 &
            !     .or. partial_internal_Energy(2) < MOOD_eps2) then
            !     ! call s_mpi_abort('PAD criteria violated')
            !     print *, 'PAD criteria violated by partial_internal_Energy ', internal_Energy
            ! end if
            ! if (rho_c2 < MOOD_eps2) then
            !     ! call s_mpi_abort('PAD criteria violated')
            !     print *, 'PAD criteria violated by rho_c2 ', rho_c2
            ! end if
            ! if (rho_c2_k(1) < MOOD_eps2 &
            !     .or. rho_c2_k(2) < MOOD_eps2) then
            !     ! call s_mpi_abort('PAD criteria violated')
            !     print *, 'PAD criteria violated by rho_c2_k ', rho_c2_k
            ! end if
            MOOD_indicator(j-neighbors:j+neighbors, k, l) = 1d0
            MOOD_indicator(j, k-neighbors:k+neighbors, l) = 1d0
            MOOD_indicator(j, k, l-neighbors:l+neighbors) = 1d0

        endif

    end subroutine s_mood_nad_pad ! ----------------------------------------

    !> Strang splitting scheme with 3rd order TVD RK time-stepping algorithm for
        !!      the flux term and adaptive time stepping algorithm for
        !!      the source term
        !! @param t_step Current time-step
    subroutine s_strang_splitting(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        ! Stage 1 of 3 =====================================================
        call s_3rd_order_tvd_rk(t_step, time_avg, 0.5d0*dt)

        ! Stage 2 of 3 =====================================================
        call s_adaptive_dt_bubble(t_step)

        ! Stage 3 of 3 =====================================================
        call s_3rd_order_tvd_rk(t_step, time_avg, 0.5d0*dt)

        call nvtxEndRange

        call cpu_time(finish)

        time = time + (finish - start)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_strang_splitting ! ------------------------------------

    !> Bubble source part in Strang operator splitting scheme
        !! @param q_cons_vf conservative variables
    subroutine s_adaptive_dt_bubble(t_step) ! ------------------------

        integer, intent(IN) :: t_step

        type(int_bounds_info) :: ix, iy, iz
        type(vector_field) :: gm_alpha_qp

        integer :: i, j, k, l, q !< Generic loop iterator

        ix%beg = 0; iy%beg = 0; iz%beg = 0
        ix%end = m; iy%end = n; iz%end = p
        call s_convert_conservative_to_primitive_variables( &
            q_cons_ts(1)%vf, &
            q_prim_vf, &
            gm_alpha_qp%vf, &
            ix, iy, iz)

        call s_compute_bubble_source(q_cons_ts(1)%vf, q_prim_vf, t_step, rhs_vf)

    end subroutine s_adaptive_dt_bubble ! ------------------------------

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step) ! ----------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        do i = 1, sys_size
            !$acc update host(q_prim_vf(i)%sf)
        end do

        if (t_step == t_step_start) then
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 1) then
            do i = 1, sys_size
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 2) then
            do i = 1, sys_size
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 3) then
            do i = 1, sys_size
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        else ! All other timesteps
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_ts(2)%vf(i)%sf(:, :, :)
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_ts(1)%vf(i)%sf(:, :, :)
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_ts(0)%vf(i)%sf(:, :, :)
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        end if

    end subroutine s_time_step_cycling ! -----------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_time_steppers_module() ! -------------------------

        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
        do i = 1, num_ts

            do j = 1, sys_size
                @:DEALLOCATE(q_cons_ts(i)%vf(j)%sf)
            end do

            @:DEALLOCATE(q_cons_ts(i)%vf)

        end do

        @:DEALLOCATE_GLOBAL(q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size
                    @:DEALLOCATE(q_prim_ts(i)%vf(j)%sf)
                end do
                @:DEALLOCATE(q_prim_ts(i)%vf)
            end do
            @:DEALLOCATE_GLOBAL(q_prim_ts)
        end if

        ! Deallocating the cell-average primitive variables
        do i = 1, adv_idx%end
            @:DEALLOCATE(q_prim_vf(i)%sf)
        end do

        if (hypoelasticity) then
            do i = stress_idx%beg, stress_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        @:DEALLOCATE_GLOBAL(q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
            @:DEALLOCATE(rhs_vf(i)%sf)
        end do

        @:DEALLOCATE_GLOBAL(rhs_vf)

        ! Deallocating the cell-average RHS variables with low-order scheme
        do i = 1, sys_size
            @:DEALLOCATE(lo_rhs_vf(i)%sf)
        end do
        @:DEALLOCATE_GLOBAL(lo_rhs_vf)

        ! Deallocating the MOOD indicator variables
        @:DEALLOCATE_GLOBAL(MOOD_indicator)

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module ! -----------------------

end module m_time_steppers
