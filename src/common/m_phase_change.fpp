!>
!! @file m_phase_change.fpp
!! @brief Contains module m_phasechange

#:include 'macros.fpp'

!> @brief This module is used to relax the model equations (6-eqn model)
!> towards pressure and temperature (6-eqn to 4-eqn), and (if wanted) Gibbs free
!> energies (6-eqn to 4-eqn) equilibrium through an infinitely fast (algebraic)
!> procedure.
module m_phase_change

#ifndef MFC_POST_PROCESS

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use ieee_arithmetic

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_phasechange_module, &
 s_relaxation_solver, &
 s_infinite_relaxation_k, &
 s_finalize_relaxation_solver_module

    !> @name Abstract interface for creating function pointers
    !> @{
    abstract interface

        !> @name Abstract subroutine for the infinite relaxation solver
        !> @{
        subroutine s_abstract_relaxation_solver(q_cons_vf) ! -------
            import :: scalar_field, sys_size
            type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        end subroutine
        !> @}

    end interface
    !> @}

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter :: max_iter = 1e8        !< max # of iterations
    real(kind(0d0)), parameter :: pCr = 4.94d7   !< Critical water pressure
    real(kind(0d0)), parameter :: TCr = 385.05 + 273.15  !< Critical water temperature
    real(kind(0d0)), parameter :: pCr_dodecane = 1.82d6   !< Critical dodecane pressure
    real(kind(0d0)), parameter :: TCr_dodecane = 658.2   !< Critical dodecane temperature
    real(kind(0d0)), parameter :: mixM = 1.0d-8 !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter :: lp = 1    !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2    !< index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    real(kind(0d0)) :: A, B, C, D
    !> @}

    !$acc declare create(max_iter,pCr,TCr,mixM,lp,vp,A,B,C,D)

    procedure(s_abstract_relaxation_solver), pointer :: s_relaxation_solver => null()

contains

    !>  The purpose of this subroutine is to initialize the phase change module
        !!      by setting the parameters needed for phase change and
        !!      selecting the phase change module that will be used
        !!      (pT- or pTg-equilibrium)
    subroutine s_initialize_phasechange_module()
        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) &
             + qvps(vp) - qvps(lp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

        D = ((gs_min(lp) - 1.0d0)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

    end subroutine s_initialize_phasechange_module !-------------------------------

    !>  This subroutine is created to activate either the pT- (N fluids) or the
        !!      pTg-equilibrium (2 fluids for g-equilibrium)
        !!      model, also considering mass depletion, depending on the incoming
        !!      state conditions.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_infinite_relaxation_k(q_cons_vf) ! ----------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)) :: pS, pSOV, pSSL !< equilibrium pressure for mixture, overheated vapor, and subcooled liquid
        real(kind(0.0d0)) :: TS, TSOV, TSSL, TSatOV, TSatSL !< equilibrium temperature for mixture, overheated vapor, and subcooled liquid. Saturation Temperatures at overheated vapor and subcooled liquid
        real(kind(0.0d0)) :: rhoe, dynE, rhos !< total internal energy, kinetic energy, and total entropy
        real(kind(0.0d0)) :: rho, rM, m1, m2, MCT !< total density, total reacting mass, individual reacting masses
        real(kind(0.0d0)) :: TvF !< total volume fraction

        !$acc declare create(pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, dynE, rhos, rho, rM, m1, m2, MCT, TvF)

        real(kind(0d0)), dimension(num_fluids) :: p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok, Tk, alphak

        !< Generic loop iterators
        integer :: i, j, k, l, q
        integer :: relax_pT !< Relaxation procedure determination variable
        real(kind(0.0d0)) :: vof_threshold
        real(kind(0.0d0)) :: p_before_thermal, rhoq_before_thermal, Gamma_before_thermal, Pi_before_thermal
        real(kind(0.0d0)) :: p_before_chemical, T_before_chemical
        real(kind(0.0d0)), dimension(sys_size) :: q_before_chemical
        !$acc declare create(p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok, Tk)

        ! starting equilibrium solver
        !$acc parallel loop collapse(3) gang vector default(present) private(p_infOV, p_infpT, p_infSL, sk, hk, gk, ek, rhok,Tk,pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, dynE, rhos, rho, rM, m1, m2, MCT, TvF)
        do j = 0, m
            do k = 0, n
                do l = 0, p

                    rho = 0.0d0; TvF = 0.0d0
                    rhoq_before_thermal = 0.0d0
                    Gamma_before_thermal = 0.0d0
                    Pi_before_thermal = 0.0d0
                    !$acc loop seq
                    do i = 1, num_fluids

                        ! Mixture density
                        rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! Total Volume Fraction
                        TvF = TvF + q_cons_vf(i + advxb - 1)%sf(j, k, l)

                        rhoq_before_thermal = rhoq_before_thermal + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                        Gamma_before_thermal = Gamma_before_thermal + q_cons_vf(i + advxb - 1)%sf(j, k, l)*gammas(i)

                        Pi_before_thermal = Pi_before_thermal + q_cons_vf(i + advxb - 1)%sf(j, k, l)*pi_infs(i)

                    end do

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process.
                    rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! correcting negative (recating) mass fraction values in case they happen
                    call s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)

                    ! fixing m1 and m2 AFTER correcting the partial densities. Note that these values must be stored for the phase
                    ! change process that will happen a posteriori
                    m1 = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

                    m2 = q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! kinetic energy as an auxiliary variable to the calculation of the total internal energy
                    dynE = 0.0d0
                    !$acc loop seq
                    do i = momxb, momxe

                        dynE = dynE + 5.0d-1*q_cons_vf(i)%sf(j, k, l)**2/rho

                    end do

                    ! calculating the total energy that MUST be preserved throughout the pT- and pTg-relaxation procedures
                    ! at each of the cells. The internal energy is calculated as the total energy minus the kinetic
                    ! energy to preserved its value at sharp interfaces
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    vof_threshold = 1.0d-4
                    if (q_cons_vf(lp + advxb - 1)%sf(j, k, l) >= vof_threshold .and. &
                        q_cons_vf(lp + advxb - 1)%sf(j, k, l) <= 1.0d0 - vof_threshold) then

                        relax_pT = 1

                    else 
                        relax_pT = 0
                    endif

                    if (relax_pT == 1) then
                        ! Calling pT-equilibrium for either finishing phase-change module, or as an IC for the pTg-equilibrium
                        ! for this case, MFL cannot be either 0 or 1, so I chose it to be 2
                        ! call s_infinite_pt_relaxation_k(j, k, l, 2, pS, p_infpT, rM, q_cons_vf, rhoe, TS)
                        
                        p_before_thermal = (rhoe - rhoq_before_thermal - Pi_before_thermal)/Gamma_before_thermal
                        ! print *, 'p_before_thermal', p_before_thermal

                        call s_infinite_pt_Wakimura2024(j, k, l, pS, p_infpT, m1, m2, rM, q_cons_vf, rhoe, TS, alphak)
                        ! call s_infinite_pt_relaxation_k(j, k, l, 2, pS, p_infpT, rM, q_cons_vf, rhoe, TS)

                        ! call s_infinite_pt_anyrate(j, k, l, pS, p_infpT, m1, m2, rM, q_cons_vf, rhoe, TS, Tk)

                        ! print *, 'difference: ', abs(TS-Tk(1)) + abs(TS-Tk(2))

                        ! Calculations AFTER pT equilibrium
                        ! entropy
                        sk(1:num_fluids) = cvs(1:num_fluids)*DLOG((TS**gs_min(1:num_fluids)) &
                                /((pS + ps_inf(1:num_fluids))**(gs_min(1:num_fluids) - 1.0d0))) &
                                + qvps(1:num_fluids)

                        ! enthalpy
                        hk(1:num_fluids) = gs_min(1:num_fluids)*cvs(1:num_fluids)*TS &
                                            + qvs(1:num_fluids)

                        ! Gibbs-free energy
                        gk(1:num_fluids) = hk(1:num_fluids) - TS*sk(1:num_fluids)

                        ! print *, gk(1), gk(2), gk(1) - gk(2)
                        ! densities
                        rhok(1) = m1/alphak(1)
                        rhok(2) = m2/alphak(2)

                        ! internal energy
                        ek(1:num_fluids) = (pS + gs_min(1:num_fluids) &
                                            *ps_inf(1:num_fluids))/(gs_min(1:num_fluids) - 1d0) &
                                           /rhok(1:num_fluids) + qvs(1:num_fluids)

                        ! calculating volume fractions, internal energies
                        rhos = 0.0d0
                        !$acc loop seq
                        do i = 1, num_fluids

                            ! volume fractions
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = alphak(i)

                            ! alpha*E
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)*ek(i) &
                                + dynE/rho*q_cons_vf(i + contxb - 1)%sf(j, k, l)

                            ! Total entropy
                            rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)

                        end do
                        ! q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(intxb)%sf(j, k, l) + q_cons_vf(intxe)%sf(j, k, l)

                        do i = 1, sys_size
                            if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l)) .or. &
                                ieee_is_nan(pS) .or. ieee_is_nan(TS)) then
                                print *, 'troubled cell caused by NAD in thermal relaxation step', i
                                print *, 'pS after thermal relaxation', pS
                                print *, 'TS after thermal relaxation', TS
                                do q = 1, sys_size
                                    print *, 'q_cons_vf', q_cons_vf(q)%sf(j, k, l)
                                end do
                                call s_mpi_abort('NAD criteria violated')
                            end if
                        end do
                    else
                        
                        ! calculate the current pressure and temperature without heat transfer
                        call s_compute_p_T(pS, TS, q_cons_vf, rhoe, j, k, l)

                    endif

                    ! ! Calling pTg-equilibrium for finishing phase-change module
                    ! if (relax_model == 6 &
                    !     .and. (pS < pCr) .and. (TS < TCr) &
                    !     .and. q_cons_vf(lp + advxb - 1)%sf(j, k, l) >= vof_threshold &
                    !     .and. q_cons_vf(lp + advxb - 1)%sf(j, k, l) <= 1.0d0 - vof_threshold) then
                        
                    !     call s_TSat(pS, TSatSL, TS, 1)
                    ! endif
                    ! if (relax_model == 6 .and. num_fluids == 2 &
                    !         .and. TS > TSatSL &
                    !         ! .and. abs(TS - TSatSL)/TSatSL > 1.0d-2 &
                    !         .and. (pS < pCr) .and. (TS < TCr) &
                    !         .and. q_cons_vf(lp + advxb - 1)%sf(j, k, l) >= vof_threshold &
                    !         .and. q_cons_vf(lp + advxb - 1)%sf(j, k, l) <= 1.0d0 - vof_threshold &
                    !         ) then
                    !     ! print *, 'Phase change is happening'
                    !     ! calling the pTg-equilibrium solver
                    !     ! print *, 'p_pT: ','pS:', pS, 'TS:', TS
                    !     call s_infinite_ptg_Wakimura(j, k, l, pS, p_infpT, rM, rhoe, q_cons_vf, TS, p_before_thermal)
                    !     ! call s_infinite_ptg_anyrate(j, k, l, pS, p_infpT, rM, rhoe, q_cons_vf, TS)
                    !     ! print *, 'p_pT_pTg: ','pS:', pS, 'TS:', TS

                    !     ! Calculations AFTER equilibrium

                    !     ! density
                    !     rhok(1:num_fluids) = (pS + ps_inf(1:num_fluids)) &
                    !         /((gs_min(1:num_fluids) - 1)*cvs(1:num_fluids)*TS)

                    !     ! volume fractions
                    !     alphak(1) = (rM - rhok(vp))/(rhok(lp) - rhok(vp))
                    !     alphak(2) = 1.0d0 - alphak(1)

                    !     ! internal energy
                    !     ek(1:num_fluids) = (pS + gs_min(1:num_fluids) &
                    !         *ps_inf(1:num_fluids))/(gs_min(1:num_fluids) - 1d0) &
                    !         /rhok(1:num_fluids) + qvs(1:num_fluids)

                    !     ! entropy
                    !     sk(1:num_fluids) = cvs(1:num_fluids)*DLOG((TS**gs_min(1:num_fluids)) &
                    !         /((pS + ps_inf(1:num_fluids))**(gs_min(1:num_fluids) - 1.0d0))) + qvps(1:num_fluids)

                    !     ! enthalpy
                    !     hk(1:num_fluids) = gs_min(1:num_fluids)*cvs(1:num_fluids)*TS &
                    !         + qvs(1:num_fluids)

                    !     ! Gibbs-free energy
                    !     gk(1:num_fluids) = hk(1:num_fluids) - TS*sk(1:num_fluids)

                    !     ! calculating volume fractions, internal energies, and total entropy
                    !     rhos = 0.0d0
                    !     !$acc loop seq
                    !     do i = 1, num_fluids

                    !         ! volume fractions
                    !         q_cons_vf(i + advxb - 1)%sf(j, k, l) = alphak(i)

                    !         ! alpha*rho
                    !         q_cons_vf(i + contxb - 1)%sf(j, k, l) = alphak(i)*rhok(i)

                    !         ! alpha*E
                    !         q_cons_vf(i + intxb - 1)%sf(j, k, l) = alphak(i)*rhok(i)*ek(i) &
                    !                                             + dynE/rho*alphak(i)*rhok(i)

                    !         ! Total entropy
                    !         rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)

                    !     end do

                    !     ! q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(intxb)%sf(j, k, l) + q_cons_vf(intxe)%sf(j, k, l)

                    !     do i = 1, sys_size
                    !         if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l)) .or. &
                    !             ieee_is_nan(pS) .or. ieee_is_nan(TS)) then
                    !             print *, 'troubled cell caused by NAD in chemical relaxation step', i
                    !             print *, 'pS before chemical relaxation', p_before_chemical
                    !             print *, 'TS before chemical relaxation', T_before_chemical
                    !             do q = 1, sys_size
                    !                 print *, 'q_before_chemical', q_before_chemical(q)
                    !             end do
                    !             print *, 'pS after chemical relaxation', pS
                    !             print *, 'TS after chemical relaxation', TS
                    !             do q = 1, sys_size
                    !                 print *, 'q_after_chemical', q_cons_vf(q)%sf(j, k, l)
                    !             end do
                    !             call s_mpi_abort('NAD criteria violated')
                    !         end if
                    !     end do


                    ! end if

                    ! check if pTg-equilibrium is required
                    ! NOTE that NOTHING else needs to be updated OTHER than the individual partial densities
                    ! given the outputs from the pT- and pTg-equilibrium solvers are just p and one of the partial masses
                    ! (pTg- case)
                    if (relax_model == 6 &
                        .and. (pS < pCr) .and. (TS < TCr) &
                        .and. q_cons_vf(lp + advxb - 1)%sf(j, k, l) >= vof_threshold &
                        .and. q_cons_vf(lp + advxb - 1)%sf(j, k, l) <= 1.0d0 - vof_threshold) then
                        
                        call s_TSat(pS, TSatSL, TS, 1)
                    endif
                    if ((relax_model == 6) &
                        ! .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) > mixM*rM) &
                        ! .and. (q_cons_vf(vp + contxb - 1)%sf(j, k, l) > mixM*rM)) &
                        .and. (pS < pCr) .and. (TS < TCr) &
                        .and. q_cons_vf(lp + advxb - 1)%sf(j, k, l) >= vof_threshold &
                        .and. q_cons_vf(lp + advxb - 1)%sf(j, k, l) <= 1.0d0 - vof_threshold &
                        .and. (TS > TSatSL) &
                        ! .and. abs(TS - TSatSL)/TSatSL > 1.0d-2 &
                        ) then

                        p_before_chemical = pS
                        T_before_chemical = TS
                        do q = 1, sys_size
                            q_before_chemical(q) = q_cons_vf(q)%sf(j, k, l)
                        end do

                        ! Checking if phase change is needed, by checking whether the final solution is either subcoooled
                        ! liquid or overheated vapor.

                        ! overheated vapor case
                        ! depleting the mass of liquid
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! tranferring the total mass to vapor
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        ! calling pT-equilibrium for overheated vapor, which is MFL = 0
                        call s_infinite_pt_relaxation_k(j, k, l, 0, pSOV, p_infOV, rM, q_cons_vf, rhoe, TSOV)
                        ! call s_infinite_pt_Wakimura2024(j, k, l, pSOV, p_infOV, mixM*rM, (1.0d0 - mixM)*rM, &
                        !                             rM, q_cons_vf, rhoe, TSOV, alphak)
                        
                        ! calculating Saturation temperature
                        call s_TSat(pSOV, TSatOV, TSOV, 2)

                        ! subcooled liquid case
                        ! tranferring the total mass to liquid
                        q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        ! depleting the mass of vapor
                        q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        ! calling pT-equilibrium for subcooled liquid, which is MFL = 1
                        call s_infinite_pt_relaxation_k(j, k, l, 1, pSSL, p_infSL, rM, q_cons_vf, rhoe, TSSL)
                        ! call s_infinite_pt_Wakimura2024(j, k, l, pSSL, p_infSL, (1.0d0 - mixM)*rM, mixM*rM, &
                        !                             rM, q_cons_vf, rhoe, TSSL, alphak)
                        
                        ! calculating Saturation temperature
                        call s_TSat(pSSL, TSatSL, TSSL, 3)

                        ! checking the conditions for overheated vapor and subcooled liquide
                        if (TSOV > TSatOV) then

                            ! Assigning pressure
                            pS = pSOV

                            ! Assigning Temperature
                            TS = TSOV

                            ! correcting the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                        elseif (TSSL < TSatSL) then

                            ! Assigning pressure
                            pS = pSSL

                            ! Assigning Temperature
                            TS = TSSL

                            ! correcting the liquid partial density
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM

                            ! correcting the vapor partial density
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                        else

                            ! returning partial pressures to what they were from the homogeneous solver
                            ! liquid
                            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = m1

                            ! vapor
                            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = m2

                            ! calling the pTg-equilibrium solver
                            call s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)
                            ! call s_infinite_ptg_anyrate(j, k, l, pS, p_infpT, rM, rhoe, q_cons_vf, TS)
                            ! call s_infinite_ptg_Wakimura(j, k, l, pS, p_infpT, rM, rhoe, q_cons_vf, TS, p_before_thermal)
                            

                        end if
                        
                        ! Calculations AFTER equilibrium

                        ! density
                        rhok(1:num_fluids) = (pS + ps_inf(1:num_fluids)) &
                            /((gs_min(1:num_fluids) - 1)*cvs(1:num_fluids)*TS)

                        ! volume fractions
                        alphak(1) = (rM - rhok(vp))/(rhok(lp) - rhok(vp))
                        alphak(2) = 1.0d0 - alphak(1)

                        ! internal energy
                        ek(1:num_fluids) = (pS + gs_min(1:num_fluids) &
                            *ps_inf(1:num_fluids))/(gs_min(1:num_fluids) - 1d0) &
                            /rhok(1:num_fluids) + qvs(1:num_fluids)

                        ! entropy
                        sk(1:num_fluids) = cvs(1:num_fluids)*DLOG((TS**gs_min(1:num_fluids)) &
                            /((pS + ps_inf(1:num_fluids))**(gs_min(1:num_fluids) - 1.0d0))) + qvps(1:num_fluids)

                        ! enthalpy
                        hk(1:num_fluids) = gs_min(1:num_fluids)*cvs(1:num_fluids)*TS &
                            + qvs(1:num_fluids)

                        ! Gibbs-free energy
                        gk(1:num_fluids) = hk(1:num_fluids) - TS*sk(1:num_fluids)

                        ! calculating volume fractions, internal energies, and total entropy
                        rhos = 0.0d0
                        !$acc loop seq
                        do i = 1, num_fluids

                            ! volume fractions
                            q_cons_vf(i + advxb - 1)%sf(j, k, l) = alphak(i)

                            ! alpha*rho
                            q_cons_vf(i + contxb - 1)%sf(j, k, l) = alphak(i)*rhok(i)

                            ! alpha*E
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = alphak(i)*rhok(i)*ek(i) &
                                                                + dynE/rho*alphak(i)*rhok(i)

                            ! Total entropy
                            rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)

                        end do

                        q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(intxb)%sf(j, k, l) + q_cons_vf(intxe)%sf(j, k, l)

                        do i = 1, sys_size
                            if (ieee_is_nan(q_cons_vf(i)%sf(j, k, l)) .or. &
                                ieee_is_nan(pS) .or. ieee_is_nan(TS)) then
                                print *, 'troubled cell caused by NAD in chemical relaxation step', i
                                print *, 'pS before chemical relaxation', p_before_chemical
                                print *, 'TS before chemical relaxation', T_before_chemical
                                do q = 1, sys_size
                                    print *, 'q_before_chemical', q_before_chemical(q)
                                end do
                                print *, 'pS after chemical relaxation', pS
                                print *, 'TS after chemical relaxation', TS
                                do q = 1, sys_size
                                    print *, 'q_after_chemical', q_cons_vf(q)%sf(j, k, l)
                                end do
                                call s_mpi_abort('NAD criteria violated')
                            end if
                        end do

                    end if

                end do
            end do
        end do

    end subroutine s_infinite_relaxation_k ! ----------------

    
    !>  This auxiliary subroutine is created to activate the pT-equilibrium for liquid-vapor fluids
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param MFL flag that tells whether the fluid is pure gas (0), pure liquid (1), or a mixture (2)
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rM sum of the reacting masses
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
        !!  @param TS equilibrium temperature at the interface
        !!  @reference Wakimura, 2024, JCP, this paper use a simple quadratic equation to solve the pT-equilibrium
    subroutine s_infinite_pt_Wakimura2024(j, k, l, pS, p_infpT, m1, m2, rM, q_cons_vf, rhoe, TS, alphak)

        !$acc routine seq

        ! initializing variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(OUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(OUT) :: p_infpT, alphak
        real(kind(0.0d0)), intent(IN) :: m1, m2, rM, rhoe
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: pk !< individual initial pressures
        integer, dimension(num_fluids) :: ig !< flags to toggle the inclusion of fluids for the pT-equilibrium

        real(kind(0.0d0)) :: mCP, mQ !< variables for the Newton Solver

        real(kind(0.0d0)) :: a0, b0, d0, Delta0, Y01, Y02, rhoe_min

        integer :: i, ns !< generic loop iterators

        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0d0; mQ = 0.0d0; p_infpT = ps_inf; 
        ! Performing tests before initializing the pT-equilibrium
        !$acc loop seq
        do i = 1, num_fluids

            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i)*gs_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

        end do

        Y01 = m1/rM
        Y02 = m2/rM
        rhoe_min = rhoe - mQ
        a0 = Y01*cvs(lp) + Y02*cvs(vp)
        b0 = Y01*cvs(lp)*(gs_min(lp)*ps_inf(lp) + ps_inf(vp) - (gs_min(lp) - 1.0d0)*rhoe_min) &
           + Y02*cvs(vp)*(gs_min(vp)*ps_inf(vp) + ps_inf(lp) - (gs_min(vp) - 1.0d0)*rhoe_min)
        d0 = Y01*cvs(lp)*ps_inf(vp)*(gs_min(lp)*ps_inf(lp) - (gs_min(lp) - 1.0d0)*rhoe_min) &
           + Y02*cvs(vp)*ps_inf(lp)*(gs_min(vp)*ps_inf(vp) - (gs_min(vp) - 1.0d0)*rhoe_min)

        Delta0 = b0**2 - 4d0*a0*d0
        if (b0 > 0) then

            pS = 2d0*d0/(-b0 - DSQRT(Delta0))

        else

            pS = (-b0 + DSQRT(Delta0))/(2d0*a0)

        end if

        ! calculate the equilibrium volume fraction first to make sure the total volume fraction is 1

        alphak(1) = cvs(1)*(gs_min(1) - 1.0d0)*(pS + ps_inf(2))*m1 / &
                    (cvs(1)*(gs_min(1) - 1.0d0)*(pS + ps_inf(2))*m1 + &
                     cvs(2)*(gs_min(2) - 1.0d0)*(pS + ps_inf(1))*m2)
        alphak(2) = 1.0d0 - alphak(1)

        ! common temperature
        TS = (pS + ps_inf(1))*alphak(1)/(cvs(1)*(gs_min(1) - 1.0d0)*m1)

        
        ! print *, 'Delta0', Delta0
        ! print *, 'a0, b0, d0', a0, b0, d0
        ! print *, 'two roots:', (-b0 + DSQRT(Delta0))/(2d0*a0), (-b0 - DSQRT(Delta0))/(2d0*a0)
        ! print *, 'pS from Wakimura2024', pS
        ! print *, 'TS from Wakimura2024', TS
        ! print *, 'alphak from Wakimura2024', alphak
        ! ! print 
        ! call s_mpi_abort('test pT relaxation Wakimura2024')

        if (alphak(1) >1d0 .or. alphak(1) < 0d0) then 
            print *, 'alpha out of range in pT-equilibrium', alphak
        end if
        if (ieee_is_nan(pS) .or. pS < 0d0) then
            print *, 'troubled cell caused by NAD in thermal relaxation step'
            print *, 'delta0', Delta0
            print *, 'equilibrium pressure after thermal relaxation', pS
            print *, 'equilibrium temperature after thermal relaxation', TS
            print *, 'alphak', alphak
            call s_mpi_abort('NAD criteria violated in thermal relaxation model')
        end if

    end subroutine s_infinite_pt_Wakimura2024 ! -----------------------

    
    !>  This auxiliary subroutine is created to activate the pT-equilibrium for liquid-vapor fluids with
        !!       There is a final common p and T after relaxation.
        !!       Method corresponding to paper `Pelanti 2022 IJMF 6eq phase change anyrate`
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param MFL flag that tells whether the fluid is pure gas (0), pure liquid (1), or a mixture (2)
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rM sum of the reacting masses
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_pt_anyrate(j, k, l, pS, p_infpT, m1, m2, rM, q_cons_vf, rhoe, TS, Tk)

        !$acc routine seq

        ! initializing variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(OUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(OUT) :: p_infpT, Tk
        real(kind(0.0d0)), intent(IN) :: m1, m2, rM, rhoe
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: pk,alphak,rhok !< individual initial pressures
        integer, dimension(num_fluids) :: ig !< flags to toggle the inclusion of fluids for the pT-equilibrium

        real(kind(0.0d0)) :: mCP, mQ !< variables for the Newton Solver

        integer :: i, ns !< generic loop iterators

        ! anyrate relaxation variables
        real(kind(0.0d0)), dimension(num_fluids) :: xi_T0, xi_T
        real(kind(0.0d0)), dimension(num_fluids) :: phi_T0, phi_T
        real(kind(0.0d0)), dimension(num_fluids) :: zeta_T0, zeta_T
        real(kind(0.0d0)), dimension(num_fluids) :: T_0, Gamma_pel
        real(kind(0.0d0)), dimension(num_fluids) :: rhoc2
        real(kind(0.0d0)) :: delta_T0, delta_T
        real(kind(0.0d0)) :: mix_gamma, mix_pi
        real(kind(0.0d0)) :: p0, D_thermal, Z, K_T, theta_T

        ! auxiliary variables for the pT-equilibrium solver
        mQ = 0.0d0; p_infpT = ps_inf; 
        mix_gamma = 0.0d0; mix_pi = 0.0d0;
        ! Performing tests before initializing the pT-equilibrium
        !$acc loop seq
        do i = 1, num_fluids
            alphak(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)
            rhok(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/alphak(i)
            Gamma_pel(i) = gs_min(i) - 1d0

            mix_gamma = mix_gamma + alphak(i)*gammas(i)
            mix_pi = mix_pi + alphak(i)*pi_infs(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

        end do

        ! calculate the initial pressure after pressure relaxation step
        p0 = (rhoe - mQ - mix_pi)/mix_gamma

        ! calculate the initial density after pressure relaxation step
        !$acc loop seq
        do i = 1, num_fluids
            ! calculate the initial temperature after pressure relaxation step
            T_0(i) = (p0 + p_infpT(i))/((gs_min(i) - 1)*cvs(i)*rhok(i))
            ! calculate the initial phi, and zeta after pressure relaxation step
            phi_T0(i) = -rhok(i)/T_0(i)
            zeta_T0(i) = 1d0/((gs_min(i) - 1)*cvs(i)*T_0(i))
            rhoc2(i) = gs_min(i)*(p0 + p_infpT(i))
        end do

        delta_T0 = T_0(2) - T_0(1)

        ! calculate the initial xi after pressure relaxation step
        Z = alphak(2)*Gamma_pel(1) + alphak(1)*Gamma_pel(2)
        D_thermal = alphak(1)*rhoc2(2) + alphak(2)*rhoc2(1)

        !$acc loop seq
        do i = 1, num_fluids
            xi_T0(i) = -(phi_T0(i)*D_thermal)/&
                        (rhok(i)/alphak(i)*Z + (-1)**(i-1)*zeta_T0(i)*&
                        (Gamma_pel(1)*rhoc2(2) - Gamma_pel(2)*rhoc2(1)))
            xi_T(i) = xi_T0(i)   
        end do
        theta_T = 2000d0
        K_T = theta_T*(1d0/xi_T(1) + 1d0/xi_T(2))
        ! print *, 'K_T: ', exp(-K_T*dt)

        ! new alpha after thermal relaxation
        alphak(1) = min(max(alphak(1) + Z/D_thermal*delta_T0/(1d0/xi_T(1) + 1d0/xi_T(2))*&   !!! finite-rate
                    (1d0 - exp(-K_T*dt)),1d-8),1d0-1d-8)
        ! alphak(1) = alphak(1) + Z/D_thermal*delta_T0/(1d0/xi_T(1) + 1d0/xi_T(2))   !!! infinite-rate
        alphak(2) = 1d0 - alphak(1)
        ! mixture variables after thermal relaxation
        mQ = 0.0d0; mix_gamma = 0.0d0; mix_pi = 0.0d0;
        ! Performing tests before initializing the pT-equilibrium
        !$acc loop seq
        do i = 1, num_fluids
            mix_gamma = mix_gamma + alphak(i)*gammas(i)
            mix_pi = mix_pi + alphak(i)*pi_infs(i)
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)
        end do
        ! calculate the pressure and temperature after thermal relaxation step
        TS = (xi_T0(1)*T_0(1) + xi_T0(2)*T_0(2))/(xi_T0(1) + xi_T0(2))
        pS = (rhoe - mQ - mix_pi)/mix_gamma

        !$acc loop seq
        do i = 1, num_fluids
            rhok(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/alphak(i)
            Tk(i) = (pS + p_infpT(i))/((gs_min(i) - 1)*cvs(i)*rhok(i))
            ! Tk(i) = TS
        end do
        ! print *, 'T1-TS: ', Tk(1)-TS, 'T2-TS: ', Tk(2)-TS

    end subroutine s_infinite_pt_anyrate ! -----------------------


    !>  This auxiliary subroutine is created to activate the pTg-equilibrium for liquid-vapor fluids.
        !!       There is a final common p and T during relaxation, and g1=g2 finally.
        !!       Method corresponding to paper `Hiro Wakimura 2024 JCP 6eq phase change BVD`
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rhoe mixture energy
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_ptg_Wakimura(j, k, l, pS, p_infpT, rM, rhoe, q_cons_vf, TS,p_before_thermal)

        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infpT
        real(kind(0.0d0)), intent(INOUT) :: pS, TS
        real(kind(0.0d0)), intent(IN) :: rhoe, rM, p_before_thermal
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: p_infpTg !< stiffness for the participating fluids for pTg-equilibrium
        real(kind(0.0d0)) :: Om ! underrelaxation factor
        real(kind(0.0d0)) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD ! auxiliary variables for the pTg-solver
        real(kind(0.0d0)) :: ap, bp, dp, Deltap, Tp
        real(kind(0.0d0)) :: gp, gpp
        real(kind(0.0d0)) :: dadp, dbdp, dddp, dsqrtDeltadp, dTdp
        real(kind(0.0d0)) :: pO !< old pressure
        real(kind(0.0d0)) :: TO !< old temperature
        real(kind(0.0d0)) :: p_before_chem, T_before_chem
        real(kind(0.0d0)) :: Teq_at_fault_p
        real(kind(0.0d0)) :: T_root_1, T_root_2

        !< Generic loop iterators
        integer :: i, ns
        ! pTg-equilibrium solution procedure
        ! Newton Solver parameters
        ! counter
        ns = 0

        ! Relaxation factor
        Om = 1.0d-3
        ! Om = 1.0d0

        p_infpTg = p_infpT

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        pO = pS
        TO = TS
        Tp = TS

        p_before_chem = pS
        T_before_chem = TS

        do while ((DABS(pS - pO) > palpha_eps) .and. (DABS((pS - pO)/pO) > palpha_eps/1e4) &
                    .and. DABS(TO - Tp) >1d-2 .or. (ns == 0))

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS
            TO = Tp

            ! coefficients for the equation of pressure
            ap = rM*cvs(lp)*cvs(vp)*((gs_min(vp) - 1d0)*(pS+gs_min(lp)*ps_inf(lp)) - &
                                     (gs_min(lp) - 1d0)*(pS+gs_min(vp)*ps_inf(vp)) )

            bp = cvs(vp)*(pS+ps_inf(lp))*(pS+gs_min(vp)*ps_inf(vp) - (gs_min(vp) - 1d0) &
                            *(rhoe - rM*qvs(lp))) - &
                 cvs(lp)*(pS+ps_inf(vp))*(pS+gs_min(lp)*ps_inf(lp) - (gs_min(lp) - 1d0) &
                            *(rhoe - rM*qvs(vp)))

            dp = (qvs(vp)-qvs(lp))*(pS+ps_inf(vp))*(pS+ps_inf(lp))

            Deltap = bp**2 - 4d0*ap*dp
            ! print *, 'ap: ', ap, 'bp: ', bp, 'dp: ', dp
            ! if (Deltap < 0) then
            !     print *, 'Deltap: ', Deltap
            !     print *, 'ap: ', ap
            !     print *, 'bp: ', bp
            !     print *, 'dp: ', dp
            ! end if
            ! print *, 'Deltap: ', Deltap

            if (Deltap < 0d0) then
                if (dsqrt(dabs(Deltap))/dabs(bp) < 1d-3) then
                    Deltap = dabs(Deltap)
                    print *, 'Deltap: ', Deltap, 'is lower than 0 and has been set to abs(Deltap)'
                else
                    print *, 'quadratic equation failed in pTg solver at pressure:', pO
                    call s_TSat(pO, Teq_at_fault_p, TS, 4)
                    print *, 'Teq_at_fault_p: ', Teq_at_fault_p
                    print *, 'error:', dsqrt(dabs(Deltap))/dabs(bp)
                    print *, 'Deltap: ', Deltap
                    print *, 'ap: ', ap
                    print *, 'bp: ', bp
                    print *, 'dp: ', dp
                    call s_mpi_abort('NAD criteria violated')
                end if

            end if

            ! derivatives of these coefficients with respect to pressure
            dadp = rM*cvs(lp)*cvs(vp)*(gs_min(vp) - gs_min(lp))
            dbdp = cvs(vp)*(pS+ps_inf(lp)+pS+gs_min(vp)*ps_inf(vp) - (gs_min(vp)-1d0) &
                            *(rhoe - rM*qvs(lp))) - &
                   cvs(lp)*(pS+ps_inf(vp)+pS+gs_min(lp)*ps_inf(lp) - (gs_min(lp)-1d0) &
                            *(rhoe - rM*qvs(vp)))
            dddp = (qvs(vp)-qvs(lp))*(pS+ps_inf(vp)+pS+ps_inf(lp))
            dsqrtDeltadp = (bp*dbdp - 2d0*ap*dddp - 2d0*dadp*dp)/DSQRT(Deltap)
            
            ! temperature and its derivative with respect to pressure
            T_root_1 = (-bp - DSQRT(Deltap))/(2d0*ap)
            T_root_2 = (-bp + DSQRT(Deltap))/(2d0*ap)

            if (abs(T_root_1 - T_before_chem) < abs(T_root_2 - T_before_chem) .and. T_root_1 > 0d0) then
                Tp = T_root_1
                dTdp = (-dbdp - dsqrtDeltadp)/(2d0*ap) &
                        - (-bp - DSQRT(Deltap))*dadp/(2d0*ap**2)
            else
                Tp = T_root_2
                dTdp = (-dbdp + dsqrtDeltadp)/(2d0*ap) &
                        - (-bp + DSQRT(Deltap))*dadp/(2d0*ap**2)
            end if

            ! if (bp > 0) then

            !     Tp = 2d0*dp/(-bp - DSQRT(Deltap))
            !     dTdp = 2d0*dddp/(-bp - DSQRT(Deltap)) &
            !             - 2d0*dp*(-dbdp - dsqrtDeltadp)/(-bp - DSQRT(Deltap))**2

            ! else

            !     Tp = (-bp + DSQRT(Deltap))/(2d0*ap)
            !     dTdp = (-dbdp + dsqrtDeltadp)/(2d0*ap) &
            !             - (-bp + DSQRT(Deltap))*dadp/(2d0*ap**2)

            ! end if
            
            ! function and its derivative with respect to pressure
            ! Note that g(pS) = 0 is the equilibrium condition
            gp = A + B/Tp + C*DLOG(Tp) + D*DLOG(pS + ps_inf(lp)) - DLOG(pS + ps_inf(vp))
            gpp = (-B/Tp**2 + C/Tp)*dTdp + D/(pS + ps_inf(lp)) - 1d0/(pS + ps_inf(vp))

            ! gp = A*Tp + B + C*Tp*DLOG(Tp) + D*Tp*DLOG(pS + ps_inf(lp)) - Tp*DLOG(pS + ps_inf(vp))
            ! gpp = (A + C + C*DLOG(Tp) + D*DLOG(pS + ps_inf(lp)) - DLOG(pS + ps_inf(vp)))*dTdp + &
            !        Tp*(D/(pS + ps_inf(lp)) - 1d0/(pS + ps_inf(vp)))

            ! updating common pressure for the newton solver
            pS = pO - Om*gp/gpp

            if (ns > 10000) then
                print *, 'old temperature:', Tp
                print *, 'old pressure:', pO
                print *, 'Updated pressure:', pS
                print *, 'coeffs under updated pressure:', 'ap', ap, 'bp', bp, 'dp', dp, '-b/2a', -bp/(2d0*ap),'b2-4ad', bp**2-4d0*ap*dp
                if (bp**2-4d0*ap*dp >= 0d0) then
                    print *, 'the solution of temperature under updated pressure is:', &
                            (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap), &
                            (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap)
                end if
                print *, 'gp:', gp
                print *, 'gpp:', gpp
                print *, 'pS-pO:', pS-pO
                print *, 'p_before_thermal', p_before_thermal
                print *, 'p_before_chem', p_before_chem
                print *, 'T_before_chem', T_before_chem
                print *, 'rho0', rM
                print *, 'rhoe0', rhoe
                print *, 'ns:', ns
                call s_mpi_abort('too many iteration times!')
            end if

        end do

        ! common temperature
        ! coefficients for the equation of pressure
        ap = rM*cvs(lp)*cvs(vp)*((gs_min(vp) - 1d0)*(pS+gs_min(lp)*ps_inf(lp)) - &
                                 (gs_min(lp) - 1d0)*(pS+gs_min(vp)*ps_inf(vp)) )

        bp = cvs(vp)*(pS+ps_inf(lp))*(pS+gs_min(vp)*ps_inf(vp) - (gs_min(vp)-1d0) &
                *(rhoe - rM*qvs(lp))) - &
             cvs(lp)*(pS+ps_inf(vp))*(pS+gs_min(lp)*ps_inf(lp) - (gs_min(lp)-1d0) &
                *(rhoe - rM*qvs(vp)))

        dp = (qvs(vp)-qvs(lp))*(pS+ps_inf(vp))*(pS+ps_inf(lp))

        Deltap = bp**2 - 4d0*ap*dp

        if (Deltap < 0d0) then
            ! if (dsqrt(dabs(Deltap))/dabs(bp) < 1d-5) then
                Deltap = dabs(Deltap)
                print *, 'Deltap: ', Deltap, 'is lower than 0 and has been set to abs(Deltap)'
            ! else
            !     print *, 'quadratic equation failed in pTg solver at pressure:', pS
            !     call s_TSat(pS, Teq_at_fault_p, TS, 4)
            !     print *, 'Teq_at_fault_p: ', Teq_at_fault_p
            !     print *, 'error:', dsqrt(dabs(Deltap))/dabs(bp)
            !     print *, 'f(T,p): ', ap*Teq_at_fault_p**2 + bp*Teq_at_fault_p + dp
            !     print *, 'Deltap: ', Deltap
            !     print *, 'ap: ', ap
            !     print *, 'bp: ', bp
            !     print *, 'dp: ', dp
            !     print *, 'pS before chemical relaxation', p_before_chem
            !     print *, 'TS before chemical relaxation', T_before_chem
            !     call s_mpi_abort('NAD criteria violated')
            ! end if

        end if

        ! temperature

        T_root_1 = (-bp - DSQRT(Deltap))/(2d0*ap)
        T_root_2 = (-bp + DSQRT(Deltap))/(2d0*ap)

        if (abs(T_root_1 - T_before_chem) < abs(T_root_2 - T_before_chem) .and. T_root_1 > 0d0) then
            TS = T_root_1
        else
            TS = T_root_2
        end if

        ! if (bp > 0) then

        !     TS = 2d0*dp/(-bp - DSQRT(Deltap))

        ! else

        !     TS = (-bp + DSQRT(Deltap))/(2d0*ap)

        ! end if

        if (ieee_is_nan(pS) .or. ieee_is_nan(TS)) then
            print *, 'troubled cell caused by NAD in chemical relaxation step'
            print *, 'pS before chemical relaxation', p_before_chem
            print *, 'TS before chemical relaxation', T_before_chem
            do i = 1, sys_size
                print *, 'q_before_chemical', q_cons_vf(i)%sf(j, k, l)
            end do
            print *, 'Deltap: ', Deltap
            print *, 'ap: ', ap
            print *, 'bp: ', bp
            print *, 'dp: ', dp
            print *, 'rho0', rM
            print *, 'rhoe0', rhoe
            print *, 'pS after chemical relaxation', pS
            print *, 'TS after chemical relaxation', TS
            print *, 'Tp within the iteration', Tp
            call s_mpi_abort('NAD criteria violated')
        end if

    end subroutine s_infinite_ptg_Wakimura !-------------------------------


    !>  This auxiliary subroutine is created to activate the pTg-equilibrium for liquid-vapor fluids.
        !!       There is a final common p and T during relaxation, and g1=g2 finally.
        !!       Method corresponding to paper `Hiro Wakimura 2024 JCP 6eq phase change BVD` with some modifications to enhance robustness
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rhoe mixture energy
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_ptg_Wakimura_new(j, k, l, pS, p_infpT, rM, rhoe, q_cons_vf, TS,p_before_thermal)

        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infpT
        real(kind(0.0d0)), intent(INOUT) :: pS, TS
        real(kind(0.0d0)), intent(IN) :: rhoe, rM, p_before_thermal
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: p_infpTg !< stiffness for the participating fluids for pTg-equilibrium
        real(kind(0.0d0)) :: Om ! underrelaxation factor
        real(kind(0.0d0)) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD ! auxiliary variables for the pTg-solver
        real(kind(0.0d0)) :: ap, bp, dp, Tp
        real(kind(0.0d0)) :: gT, gTT
        real(kind(0.0d0)) :: TO !< old temperature
        real(kind(0.0d0)) :: pO !< old pressure
        real(kind(0.0d0)) :: p_before_chem, T_before_chem
        real(kind(0.0d0)) :: Teq_at_fault_p
        real(kind(0.0d0)) :: S_index

        !< Generic loop iterators
        integer :: i, ns
        ! pTg-equilibrium solution procedure
        ! Newton Solver parameters
        ! counter
        ns = 0

        ! Relaxation factor
        Om = 1.0d-2
        ! Om = 1.0d-2

        p_infpTg = p_infpT

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        TO = TS

        p_before_chem = pS
        T_before_chem = TS

        ! initialize old temperature and pressure
        TO = TS
        pO = pS

        ! if (TO < 300d0 .or. TO > 600d0) then
        !     print *, 'initial temperature in pTg solver is out of range:', TO
            TS = 400d0
            TO = TS
        !     print *, 'we therefore reset the temperature to ', TO
        ! end if
        ! calculate the saturation pressure under given temperature
        call s_pSat(pS, TO, 2d1)
        ! print *, 'Initial pressure:', pO, 'Initial temperature:', TO, 'Saturation pressure under initial temperature:', pS
        ! print *, 'p_before_thermal', p_before_thermal
        ! call s_mpi_abort('test')

        ! initial scaling factor for the Newton solver
        S_index = 1.0d0
        do
            ! increasing counter
            ns = ns + 1

            TO = TS
            pO = pS

            ! coefficients for the equation of pressure
            ap = rM*cvs(lp)*cvs(vp)*((gs_min(vp) - 1d0)*(pS+gs_min(lp)*ps_inf(lp)) - &
                                     (gs_min(lp) - 1d0)*(pS+gs_min(vp)*ps_inf(vp)) )

            bp = cvs(vp)*(pS+ps_inf(lp))*(pS+gs_min(vp)*ps_inf(vp) - (gs_min(vp) - 1d0) &
                            *(rhoe - rM*qvs(lp))) - &
                 cvs(lp)*(pS+ps_inf(vp))*(pS+gs_min(lp)*ps_inf(lp) - (gs_min(lp) - 1d0) &
                            *(rhoe - rM*qvs(vp)))

            dp = (qvs(vp)-qvs(lp))*(pS+ps_inf(vp))*(pS+ps_inf(lp))
            
            ! function and its derivative with respect to temperature
            ! Note that g(TS) = 0 is the equilibrium condition
            if (ns == 1) then
                S_index = max( 1.0d0, max( abs(ap)*TS*TS, max(abs(bp)*TS, abs(dp)) ) )
            end if
            gT = (ap*TS**2d0 + bp*TS + dp)/S_index
            gTT = (2d0*ap*TS + bp)/S_index

            ! updating common temperature for the newton solver
            TS = TO - Om*gT/gTT
            if (TS > 681d0 .or. TS < 273d0) then
                print *, 'old temperature:', TO
                print *, 'Updated temperature:', TS
                print *, 'old pressure:', pO
                print *, 'Updated pressure:', pS
                print *, 'coeffs under updated pressure:', 'ap', ap, 'bp', bp, 'dp', dp, '-b/2a', -bp/(2d0*ap),'b2-4ad', bp**2-4d0*ap*dp
                if (bp**2-4d0*ap*dp >= 0d0) then
                    print *, 'the solution of temperature under updated pressure is:', &
                            (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap), &
                            (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap)

                    print *, 'the difference between the result and the old temperature is:', &
                            (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap) - TO, &
                            (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap) - TO
                end if
                print *, 'gT:', gT
                print *, 'gTT:', gTT
                print *, 'TS-TO:', TS-TO
                print *, 'pS-pO:', pS-pO
                print *, 'p_before_thermal', p_before_thermal
                print *, 'p_before_chem', p_before_chem
                print *, 'T_before_chem', T_before_chem
                print *, 'rho0', rM
                print *, 'rhoe0', rhoe
                print *, 'ns:', ns
                call s_mpi_abort('calculated temperature out of range!')
            end if


            call s_pSat(pS, TS, 2d1)
            if ( &
                    ((DABS(TS - TO) < palpha_eps) .or. (DABS((TS - TO)/TO) < 1d-5)) &
                    .or. (((dabs(pS-pO) < 1d-1) .or. (dabs(pS-pO)/pO < 1d-1/1d5))) &
                    .or. abs(gT) < 1d-4 &
                    .and. ns > 1) then
                
                ! print *, 'old temperature:', TO
                ! print *, 'Updated temperature:', TS
                ! print *, 'old pressure:', pO
                ! print *, 'Updated pressure:', pS
                ! print *, 'coeffs under updated pressure:', 'ap', ap, 'bp', bp, 'dp', dp, '-b/2a', -bp/(2d0*ap),'b2-4ad', bp**2-4d0*ap*dp
                ! if (bp**2-4d0*ap*dp >= 0d0) then
                !     print *, 'the solution of temperature under updated pressure is:', &
                !             (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap), &
                !             (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap)

                !     print *, 'the difference between the result and the old temperature is:', &
                !             (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap) - TO, &
                !             (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap) - TO
                ! end if
                ! print *, 'gT:', gT
                ! print *, 'gTT:', gTT
                ! print *, 'TS-TO:', TS-TO
                ! print *, 'pS-pO:', pS-pO
                ! print *, 'ns:', ns
                ! call s_mpi_abort('TEST iteration checker!')
                    
                exit
            end if
            
            if (ns > 2000) then
                print *, 'old temperature:', TO
                print *, 'Updated temperature:', TS
                print *, 'old pressure:', pO
                print *, 'Updated pressure:', pS
                print *, 'coeffs under updated pressure:', 'ap', ap, 'bp', bp, 'dp', dp, '-b/2a', -bp/(2d0*ap),'b2-4ad', bp**2-4d0*ap*dp
                if (bp**2-4d0*ap*dp >= 0d0) then
                    print *, 'the solution of temperature under updated pressure is:', &
                            (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap), &
                            (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap)

                    print *, 'the difference between the result and the old temperature is:', &
                            (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap) - TO, &
                            (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap) - TO
                end if
                print *, 'gT:', gT
                print *, 'gTT:', gTT
                print *, 'TS-TO:', TS-TO
                print *, 'pS-pO:', pS-pO
                print *, 'p_before_thermal', p_before_thermal
                print *, 'p_before_chem', p_before_chem
                print *, 'T_before_chem', T_before_chem
                print *, 'rho0', rM
                print *, 'rhoe0', rhoe
                print *, 'ns:', ns
                call s_mpi_abort('too many iteration times!')
            end if
        end do

        if (TS < 200d0 .or. TS > 1000d0) then
            print *, 'old temperature:', TO
            print *, 'Updated temperature:', TS
            print *, 'old pressure:', pO
            print *, 'Updated pressure:', pS
            print *, 'coeffs under updated pressure:', 'ap', ap, 'bp', bp, 'dp', dp, '-b/2a', -bp/(2d0*ap),'b2-4ad', bp**2-4d0*ap*dp
            if (bp**2-4d0*ap*dp >= 0d0) then
                print *, 'the solution of temperature under updated pressure is:', &
                        (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap), &
                        (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap)

                print *, 'the difference between the result and the old temperature is:', &
                        (-bp + DSQRT(bp**2-4d0*ap*dp))/(2d0*ap) - TO, &
                        (-bp - DSQRT(bp**2-4d0*ap*dp))/(2d0*ap) - TO
            end if
            print *, 'gT:', gT
            print *, 'gTT:', gTT
            print *, 'TS-TO:', TS-TO
            print *, 'pS-pO:', pS-pO
            print *, 'ns:', ns
            print *, 'p_before_thermal', p_before_thermal
            print *, 'p_before_chem', p_before_chem
            print *, 'T_before_chem', T_before_chem
            call s_mpi_abort('totally wrong temperature!')
        end if

        if (ieee_is_nan(pS) .or. ieee_is_nan(TS)) then
            print *, 'troubled cell caused by NAD in chemical relaxation step'
            print *, 'pS before chemical relaxation', p_before_chem
            print *, 'TS before chemical relaxation', T_before_chem
            do i = 1, sys_size
                print *, i, 'q_before_chemical', q_cons_vf(i)%sf(j, k, l)
            end do
            print *, 'ap: ', ap
            print *, 'bp: ', bp
            print *, 'dp: ', dp
            print *, 'rho0', rM
            print *, 'rhoe0', rhoe
            print *, 'pS after chemical relaxation', pS
            print *, 'TS after chemical relaxation', TS
            call s_mpi_abort('NAD criteria violated')
        end if

    end subroutine s_infinite_ptg_Wakimura_new !-------------------------------

    
    !>  This auxiliary subroutine is created to activate the pTg-equilibrium for liquid-vapor fluids.
        !!       There is a final common p and T during relaxation, and g1=g2 finally.
        !!       Method corresponding to paper `Hiro Wakimura 2024 JCP 6eq phase change BVD`
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rhoe mixture energy
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_ptg_anyrate(j, k, l, pS, p_infpT, rM, rhoe, q_cons_vf, TS)

        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infpT
        real(kind(0.0d0)), intent(INOUT) :: pS, TS
        real(kind(0.0d0)), intent(IN) :: rhoe, rM
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: p_infpTg !< stiffness for the participating fluids for pTg-equilibrium
        real(kind(0.0d0)) :: mCP, mQ ! auxiliary variables for the pTg-solver

        !< Generic loop iterators
        integer :: i, ns
        ! pTg-equilibrium solution procedure

        ! Variables for the chemical relaxation step
        real(kind(0.0d0)), dimension(num_fluids) :: xi_g0, xi_g
        real(kind(0.0d0)), dimension(num_fluids) :: s_k0, s_k
        real(kind(0.0d0)), dimension(num_fluids) :: rhok, alphak, alphak0
        real(kind(0.0d0)), dimension(num_fluids) :: phi_k0, phi_k
        real(kind(0.0d0)), dimension(num_fluids) :: zeta_k0, zeta_k
        real(kind(0.0d0)), dimension(num_fluids) :: Gamma_pel, chi_pel
        real(kind(0.0d0)), dimension(num_fluids) :: rhoc2
        real(kind(0.0d0)), dimension(num_fluids) :: gk0, gk
        real(kind(0.0d0)), dimension(num_fluids) :: alpharhok0, alpharhok
        real(kind(0.0d0)) :: S_alpha, S_p, S_T, D_T
        real(kind(0.0d0)) :: delta_g0, delta_g
        real(kind(0.0d0)) :: mix_gamma, mix_pi
        real(kind(0.0d0)) :: K_g, hat_K_g, nu_g
        ! Newton Solver parameters

        ! auxiliary variables for the pT-equilibrium solver
        mQ = 0.0d0; p_infpTg = p_infpT; 
        mix_gamma = 0.0d0; mix_pi = 0.0d0;
        ! calculate initial variables after thermal relaxation step
        !$acc loop seq
        do i = 1, num_fluids
            ! alpha, rho
            alphak(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)
            alphak0(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)
            rhok(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/alphak(i)
            alpharhok0(i) = alphak(i)*rhok(i)

            ! Gamma, chi
            Gamma_pel(i) = (gs_min(i) - 1)
            chi_pel(i) = -(gs_min(i) - 1)*qvs(i)

            ! phi, zeta, rho*c^2
            phi_k0(i) = -rhok(i)/TS
            zeta_k0(i) = 1d0/((gs_min(i) - 1)*cvs(i)*TS)
            rhoc2(i) = gs_min(i)*(pS + p_infpTg(i))

            mix_gamma = mix_gamma + alphak0(i)*gammas(i)
            mix_pi = mix_pi + alphak0(i)*pi_infs(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ +alpharhok0(i)*qvs(i)

        end do

        D_T = alphak(1)*alphak(2)*&
                    (rhoc2(1)/Gamma_pel(1) - rhoc2(2)/Gamma_pel(2))*&
                    (phi_k0(1)*zeta_k0(2) - phi_k0(2)*zeta_k0(1)) + &
                    (alphak(1)/Gamma_pel(1) + alphak(2)/Gamma_pel(2))*&
                    (alphak(1)*phi_k0(1)*rhok(2) + alphak(2)*phi_k0(2)*rhok(1))

        S_alpha = 1d0/D_T*(&
                    (alphak(1)/Gamma_pel(1) + alphak(2)/Gamma_pel(2))*&
                    (alphak(1)*phi_k0(1) + alphak(2)*phi_k0(2)) + alphak(1)*alphak(2)*&
                    (phi_k0(1)*zeta_k0(2) - phi_k0(2)*zeta_k0(1))*&
                    (chi_pel(1)/Gamma_pel(1) - chi_pel(2)/Gamma_pel(2)))

        S_T = 1d0/D_T*(&
                    (chi_pel(2)/Gamma_pel(2) - chi_pel(1)/Gamma_pel(1))*&
                    (alphak(1)*zeta_k0(1)*rhok(2) + alphak(2)*zeta_k0(2)*rhok(1)) +&
                    (rhoc2(1)/Gamma_pel(1) - rhoc2(2)/Gamma_pel(2))*&
                    (alphak(1)*zeta_k0(1) + alphak(2)*zeta_k0(2)) +&
                    (rhok(2) - rhok(1))*(alphak(1)/Gamma_pel(1) + alphak(2)/Gamma_pel(2)))

        S_p = 1d0/D_T*(&
                    (chi_pel(1)/Gamma_pel(1) - chi_pel(2)/Gamma_pel(2))*&
                    (alphak(1)*phi_k0(1)*rhok(2) + alphak(2)*phi_k0(2)*rhok(1)) +&
                    (rhoc2(2)/Gamma_pel(2) - rhoc2(1)/Gamma_pel(1))*&
                    (alphak(1)*phi_k0(1) + alphak(2)*phi_k0(2)))

        ! constant variables for the ODE
        !$acc loop seq
        do i = 1, num_fluids

            s_k0(i) = cvs(i)*log(TS**gs_min(i)/(pS + p_infpTg(i))**(gs_min(i) - 1)) + qvps(i)

            gk0(i) = (gs_min(i)*cvs(i) - qvps(i))*TS - &
                    cvs(i)*TS*log(TS**gs_min(i)/(pS + p_infpTg(i))**(gs_min(i) - 1)) + qvs(i)

            xi_g0(i) = 1d0/((-1)**(i-1)/rhok(i)*S_p - s_k0(i)*S_T)

            xi_g(i) = xi_g0(i)

        enddo

        nu_g = 1d-4
        K_g = nu_g*(1d0/xi_g(1) + 1d0/xi_g(2))
        delta_g0 = gk0(2) - gk0(1)

        !$acc loop seq
        do i = 1, num_fluids

            alpharhok(i) = alpharhok0(i) + (-1)**(i-1)*delta_g0/(1d0/xi_g(1) + 1d0/xi_g(2))*&
                            (1d0 - exp(-K_g*dt))

            alphak(i) = alphak0(i) + (-1)**(i-1)*S_alpha*delta_g0/(1d0/xi_g(1) + 1d0/xi_g(2))*&
                            (1d0 - exp(-K_g*dt))

        enddo

        ! modify the volume fraction if it is not reasonable
        if ( alphak(1) < 1d-8 .or. alphak(1) > 1d0 - 1d-8 ) then
            if (s_k0(2) > s_k0(1)) then
                alphak(1) = 1d-8
                alphak(2) = 1d0 - 1d-8
            else
                alphak(1) = 1d0 - 1d-8
                alphak(2) = 1d-8
            end if

            hat_K_g = -1d0/dt*&
                log(1d0 - (alphak(1) - alphak0(1))/&
                (S_alpha*delta_g0/(1d0/xi_g(1) + 1d0/xi_g(2))))

            alpharhok(1) = alpharhok0(1) + delta_g0/(1d0/xi_g(1) + 1d0/xi_g(2))*&
                (1d0 - exp(-hat_K_g*dt))
            alpharhok(2) = alpharhok0(2) - delta_g0/(1d0/xi_g(1) + 1d0/xi_g(2))*&
                (1d0 - exp(-hat_K_g*dt))

        endif

        ! calculate the mixture parameters after chemical relaxation step
        mCP = 0.0d0; mQ = 0.0d0; mix_gamma = 0.0d0; mix_pi = 0.0d0;
        !$acc loop seq
        do i = 1, num_fluids

            mCP = mCP + alpharhok(i)*cvs(i)*gs_min(i)

            mQ = mQ + alpharhok(i)*qvs(i)

            mix_gamma = mix_gamma + alphak(i)*gammas(i)

            mix_pi = mix_pi + alphak(i)*pi_infs(i)

        enddo

        ! update the partial densities after chemical relaxation step
        q_cons_vf(1 + contxb - 1)%sf(j, k, l) = alpharhok(1)
        q_cons_vf(2 + contxb - 1)%sf(j, k, l) = alpharhok(2)

        ! calculate the pressure after chemical relaxation step
        pS = (rhoe - mQ - mix_pi)/mix_gamma

        ! calculate the temperature after chemical relaxation step
        TS = (rhoe + pS - mQ)/mCP


    end subroutine s_infinite_ptg_anyrate !-------------------------------

    !>  This auxiliary subroutine is created to activate the pT-equilibrium for N fluids
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param MFL flag that tells whether the fluid is pure gas (0), pure liquid (1), or a mixture (2)
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rM sum of the reacting masses
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_pt_relaxation_k(j, k, l, MFL, pS, p_infpT, rM, q_cons_vf, rhoe, TS)

        !$acc routine seq

        ! initializing variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(OUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(OUT) :: p_infpT
        real(kind(0.0d0)), intent(IN) :: rM, rhoe
        integer, intent(IN) :: j, k, l, MFL
        real(kind(0.0d0)), dimension(num_fluids) :: pk !< individual initial pressures
        integer, dimension(num_fluids) :: ig !< flags to toggle the inclusion of fluids for the pT-equilibrium

        real(kind(0.0d0)) :: gp, gpp, hp, pO, mCP, mQ !< variables for the Newton Solver

        integer :: i, ns !< generic loop iterators

        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0d0; mQ = 0.0d0; p_infpT = ps_inf; 
        ! Performing tests before initializing the pT-equilibrium
        !$acc loop seq
        do i = 1, num_fluids

            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i)*gs_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

        end do

        ! Checking energy constraint
        if ((rhoe - mQ - minval(p_infpT)) < 0.0d0) then

            if ((MFL == 0) .or. (MFL == 1)) then

                ! Assigning zero values for mass depletion cases
                ! pressure
                pS = 0.0d0

                ! temperature
                TS = 0.0d0

                return
            end if

        end if

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        pO = 0.0d0

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(ps_inf)
        ! and infinity, a solution should be able to be found.
        pS = 1.0d4

        ! Newton Solver for the pT-equilibrium
        ns = 0
        ! change this relative error metric. 1E4 is just arbitrary
        do while ((DABS(pS - pO) > palpha_eps) .and. (DABS((pS - pO)/pO) > palpha_eps/1e4) .or. (ns == 0))

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            gpp = 0.0d0; gp = 0.0d0; hp = 0.0d0
            !$acc loop seq
            do i = 1, num_fluids

                gp = gp + (gs_min(i) - 1.0d0)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                     *(rhoe + pS - mQ)/(mCP*(pS + p_infpT(i)))

                gpp = gpp + (gs_min(i) - 1.0d0)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                      *(p_infpT(i) - rhoe + mQ)/(mCP*(pS + p_infpT(i))**2)

            end do

            hp = 1.0d0/(rhoe + pS - mQ) + 1.0d0/(pS + minval(p_infpT))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0d0 - gp)/gpp)/(1.0d0 - (1.0d0 - gp + DABS(1.0d0 - gp)) &
                                          /(2.0d0*gpp)*hp)
        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP

    end subroutine s_infinite_pt_relaxation_k ! -----------------------

    !>  This auxiliary subroutine is created to activate the pTg-equilibrium for N fluids under pT
        !!      and 2 fluids under pTg-equilibrium. There is a final common p and T during relaxation
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rhoe mixture energy
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rhoe, q_cons_vf, TS)

        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infpT
        real(kind(0.0d0)), intent(INOUT) :: pS, TS
        real(kind(0.0d0)), intent(IN) :: rhoe
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(num_fluids) :: p_infpTg !< stiffness for the participating fluids for pTg-equilibrium
        real(kind(0.0d0)), dimension(2, 2) :: Jac, InvJac, TJac !< matrices for the Newton Solver
        real(kind(0.0d0)), dimension(2) :: R2D, DeltamP !< residual and correction array
        real(kind(0.0d0)) :: Om ! underrelaxation factor
        real(kind(0.0d0)) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD ! auxiliary variables for the pTg-solver

        !< Generic loop iterators
        integer :: i, ns
        ! pTg-equilibrium solution procedure
        ! Newton Solver parameters
        ! counter
        ns = 0

        ! Relaxation factor
        Om = 1.0d-1

        p_infpTg = p_infpT

        if (((pS < 0.0d0) .and. ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
                                  + q_cons_vf(vp + contxb - 1)%sf(j, k, l)) > ((rhoe &
                                                                                - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1))/qvs(lp)))) .or. &
            ((pS >= 0.0d0) .and. (pS < 1.0d-1))) then

            ! improve this initial condition
            pS = 1.0d4

        end if

        ! Loop until the solution for F(X) is satisfied
        ! Check whether I need to use both absolute and relative values
        ! for the residual, and how to do it adequately.
        ! Dummy guess to start the pTg-equilibrium problem.
        ! improve this initial condition
        R2D(1) = 0.0d0; R2D(2) = 0.0d0
        DeltamP(1) = 0.0d0; DeltamP(2) = 0.0d0
        do while (((DSQRT(R2D(1)**2 + R2D(2)**2) > ptgalpha_eps) &
                   .and. ((DSQRT(R2D(1)**2 + R2D(2)**2)/rhoe) > (ptgalpha_eps/1d6))) &
                  .or. (ns == 0))

            ! Updating counter for the iterative procedure
            ns = ns + 1

            ! Auxiliary variables to help in the calculation of the residue
            mCP = 0.0d0; mCPD = 0.0d0; mCVGP = 0.0d0; mCVGP2 = 0.0d0; mQ = 0.0d0; mQD = 0.0d0
            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure
            !$acc loop seq
            do i = 1, num_fluids

                ! sum of the total alpha*rho*cp of the system
                mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                      *cvs(i)*gs_min(i)

                ! sum of the total alpha*rho*q of the system
                mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                ! These auxiliary variables now need to be updated, as the partial densities now
                ! vary at every iteration
                if ((i /= lp) .and. (i /= vp)) then

                    mCVGP = mCVGP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                            *cvs(i)*(gs_min(i) - 1)/(pS + ps_inf(i))

                    mCVGP2 = mCVGP2 + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                             *cvs(i)*(gs_min(i) - 1)/((pS + ps_inf(i))**2)

                    mQD = mQD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                    ! sum of the total alpha*rho*cp of the system
                    mCPD = mCPD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                           *gs_min(i)

                end if

            end do

            ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
            call s_compute_jacobian_matrix(InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2, pS, q_cons_vf, TJac)

            ! calculating correction array for Newton's method
            DeltamP = -1.0d0*matmul(InvJac, R2D)

            ! updating two reacting 'masses'. Recall that inert 'masses' do not change during the phase change
            ! liquid
            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + Om*DeltamP(1)

            ! gas
            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = q_cons_vf(vp + contxb - 1)%sf(j, k, l) - Om*DeltamP(1)

            ! updating pressure
            pS = pS + Om*DeltamP(2)

            ! calculating residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
            ! and (ii) the energy before and after the phase-change process.
            call s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)

        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP
    end subroutine s_infinite_ptg_relaxation_k ! -----------------------

    !>  This auxiliary subroutine compute the current pressure and temperature of the system
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
    subroutine s_compute_p_T(pS, TS, q_cons_vf, rhoe, j, k, l)
        !$acc routine seq

        !> @name variables for the computation of the pressure and temperature
        !> @{
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: rhoe
        real(kind(0.0d0)), intent(OUT) :: pS, TS
        real(kind(0.0d0)) :: rho, gamma, pi_inf, mQ, mCP
        integer, intent(IN) :: j, k, l

        !< Generic loop iterators
        integer :: i
        !> @}

        rho = 0.0d0; gamma = 0.0d0; pi_inf = 0.0d0; mQ = 0.0d0; mCP = 0.0d0

        !$acc loop seq
        do i = 1, num_fluids

            rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)
            gamma = gamma + q_cons_vf(i + advxb - 1)%sf(j, k, l)*gammas(i)
            pi_inf = pi_inf + q_cons_vf(i + advxb - 1)%sf(j, k, l)*pi_infs(i)
            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i)*gs_min(i)
            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

        end do

        pS = (rhoe - mQ - pi_inf)/gamma

        TS = (rhoe + pS - mQ)/mCP

    end subroutine s_compute_p_T

    !>  This auxiliary subroutine corrects the partial densities of the REACTING fluids in case one of them is negative
        !!      but their sum is positive. Inert phases are not corrected at this moment
        !!  @param MCT partial density correction parameter
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rM sum of the reacting masses
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
    subroutine s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)
        !$acc routine seq

        !> @name variables for the correction of the reacting partial densities
        !> @{
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), intent(INOUT) :: rM
        real(kind(0.0d0)), intent(OUT) :: MCT
        integer, intent(IN) :: j, k, l
        !> @}
        if (rM < 0.0d0) then

            if ((q_cons_vf(lp + contxb - 1)%sf(j, k, l) >= -1.0d0*mixM) .and. &
                (q_cons_vf(vp + contxb - 1)%sf(j, k, l) >= -1.0d0*mixM)) then

                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = 0.0d0

                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = 0.0d0

                rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

            end if

        end if

        ! Defining the correction in terms of an absolute value might not be the best practice.
        ! Maybe a good way to do this is to partition the partial densities, giving a small percentage of the total reacting density
        MCT = 2*mixM

        ! correcting the partial densities of the reacting fluids. What to do for the nonreacting ones?
        if (q_cons_vf(lp + contxb - 1)%sf(j, k, l) < 0.0d0) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = MCT*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

        elseif (q_cons_vf(vp + contxb - 1)%sf(j, k, l) < 0.0d0) then

            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - MCT)*rM

            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = MCT*rM

        end if
    end subroutine s_correct_partial_densities

    !>  This auxiliary subroutine calculates the 2 x 2 Jacobian and, its inverse and transpose
        !!      to be used in the pTg-equilibirium procedure
        !!  @param InvJac Inverse of the Jacobian Matrix
        !!  @param j generic loop iterator for x direction
        !!  @param Jac Jacobian Matrix
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param mCPD  sum of the total alpha*rho*cp
        !!  @param mCVGP auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press
        !!  @param mCVGP2 auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press^2
        !!  @param pS equilibrium pressure at the interface
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TJac Transpose of the Jacobian Matrix
    subroutine s_compute_jacobian_matrix(InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2, pS, q_cons_vf, TJac)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, mCPD, mCVGP, mCVGP2
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(2, 2), intent(OUT) :: Jac, InvJac, TJac
        real(kind(0.0d0)) :: ml, mT, TS, dFdT, dTdm, dTdp ! mass of the reacting fluid, total reacting mass, and auxiliary variables

        ! mass of the reacting liquid
        ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
             + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

        TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        dFdT = &
            -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*DLOG(TS) &
            - (qvps(lp) - qvps(vp)) &
            + cvs(lp)*(gs_min(lp) - 1)*DLOG(pS + ps_inf(lp)) &
            - cvs(vp)*(gs_min(vp) - 1)*DLOG(pS + ps_inf(vp))

        dTdm = -(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                 - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))*TS**2

        dTdp = (mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                + mCVGP2)*TS**2

        ! F = (F1,F2) is the function whose roots we are looking for
        ! x = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
        ! F1 = 0 is the Gibbs free energy quality
        ! F2 = 0 is the enforcement of the thermodynamic (total - kinectic) energy
        ! dF1dm
        Jac(1, 1) = dFdT*dTdm

        ! dF1dp
        Jac(1, 2) = dFdT*dTdp + TS &
                    *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))

        ! dF2dm
        Jac(2, 1) = (qvs(vp) - qvs(lp) &
                     + (cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                     /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP) &
                     - (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                        - mT*cvs(vp)*gs_min(vp) - mCPD) &
                     *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                       - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                     /((ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                            - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                        + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2))/1
        ! dF2dp
        Jac(2, 2) = (1 + (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                          - mT*cvs(vp)*gs_min(vp) - mCPD) &
                     *(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 + mCVGP2) &
                     /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2)/1

        ! intermediate elements of J^{-1}
        InvJac(1, 1) = Jac(2, 2)
        InvJac(1, 2) = -1.0d0*Jac(1, 2)
        InvJac(2, 1) = -1.0d0*Jac(2, 1)
        InvJac(2, 2) = Jac(1, 1)

        ! elements of J^{T}
        TJac(1, 1) = Jac(1, 1)
        TJac(1, 2) = Jac(2, 1)
        TJac(2, 1) = Jac(1, 2)
        TJac(2, 2) = Jac(2, 2)

        ! dividing by det(J)
        InvJac = InvJac/(Jac(1, 1)*Jac(2, 2) - Jac(1, 2)*Jac(2, 1))

    end subroutine s_compute_jacobian_matrix

    !>  This auxiliary subroutine computes the residue of the pTg-equilibrium procedure
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param mCPD  sum of the total alpha*rho*cp
        !!  @param mCVGP auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press
        !!  @param mQD sum of the total alpha*rho*qv
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param pS equilibrium pressure at the interface
        !!  @param rhoe mixture energy
        !!  @param R2D (2D) residue array
    subroutine s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, rhoe, mCPD, mCVGP, mQD
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(2), intent(OUT) :: R2D
        real(kind(0.0d0)) :: ml, mT, TS !< mass of the reacting liquid, total reacting mass, equilibrium temperature

        ! mass of the reacting liquid
        ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
             + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

        TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        ! Gibbs Free Energy Equality condition (DG)
        R2D(1) = TS*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                     *(1 - DLOG(TS)) - (qvps(lp) - qvps(vp)) &
                     + cvs(lp)*(gs_min(lp) - 1)*DLOG(pS + ps_inf(lp)) &
                     - cvs(vp)*(gs_min(vp) - 1)*DLOG(pS + ps_inf(vp))) &
                 + qvs(lp) - qvs(vp)

        ! Constant Energy Process condition (DE)
        R2D(2) = (rhoe + pS &
                  + ml*(qvs(vp) - qvs(lp)) - mT*qvs(vp) - mQD &
                  + (ml*(gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
                     - mT*gs_min(vp)*cvs(vp) - mCPD) &
                  /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                        - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                    + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP))/1

    end subroutine s_compute_pTg_residue

    !>  This auxiliary subroutine finds the Saturation temperature for a given
        !!      saturation pressure through a newton solver
        !!  @param pSat Saturation Pressure
        !!  @param TSat Saturation Temperature
        !!  @param TSIn equilibrium Temperature
    subroutine s_TSat(pSat, TSat, TSIn, location_of_calculated_TSat)
        !$acc routine seq

        real(kind(0.0d0)), intent(OUT) :: TSat
        real(kind(0.0d0)), intent(IN) :: pSat, TSIn
        integer, intent(IN) :: location_of_calculated_TSat
        real(kind(0.0d0)) :: dFdT, FT, Om !< auxiliary variables

        ! Generic loop iterators
        integer :: ns

        if ((pSat == 0.0d0) .and. (TSIn == 0.0d0)) then

            ! assigning Saturation temperature
            TSat = 0.0d0

        else

            ! calculating initial estimate for temperature in the TSat procedure. I will also use this variable to
            ! iterate over the Newton's solver
            TSat = TSIn

            ! iteration counter
            ns = 0

            ! underrelaxation factor
            ! Om = 1.0d-3
            Om = 1.0d0
            FT = 0.0d0
            do while ((DABS(FT) > ptgalpha_eps) .or. (ns == 0))
                ! increasing counter
                ns = ns + 1

                ! calculating residual
                FT = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                           *(1 - DLOG(TSat)) - (qvps(lp) - qvps(vp)) &
                           + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

                ! calculating the jacobian
                dFdT = &
                    -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*DLOG(TSat) &
                    - (qvps(lp) - qvps(vp)) &
                    + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                    - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))

                ! ! calculating residual
                ! FT = A + B/TSat + C*DLOG(TSat) + D*DLOG(pSat + ps_inf(lp)) - &
                !         DLOG(pSat + ps_inf(vp))

                ! ! calculating the jacobian
                ! dFdT = -B/TSat**2 + C/TSat

                ! updating saturation temperature
                TSat = TSat - Om*FT/dFdT

                if (ns > 1000) then
                    print *, 'Iteration counter: ', ns,'TSIn:',TSIn,'pSat:', pSat,'TSat:', TSat
                    print *, 'FT:', FT, 'dFdT',dFdT
                    print *, 'A', A, 'B', B, 'C', C, 'D', D, 'B/C', B/C
                    FT = A + B/100d0 + C*DLOG(100d0) + D*DLOG(pSat + ps_inf(lp)) - &
                            DLOG(pSat + ps_inf(vp))
                    print *, 'T = 100:, F(T)=', FT
                    FT = A + B/(B/C) + C*DLOG(B/C) + D*DLOG(pSat + ps_inf(lp)) - &
                            DLOG(pSat + ps_inf(vp))
                    print *, 'T =',B/C, ':, F(T)=', FT
                    print *, 'location of difficult step to calculate TSat', location_of_calculated_TSat
                    call s_mpi_abort('too large iteration number. Exiting ...')
                endif

            end do

        end if

    end subroutine s_TSat


    !>  This auxiliary subroutine finds the Saturation pressure for a given
        !!      saturation temperature through a newton solver
        !!  @param pSat Saturation Pressure
        !!  @param TSat Saturation Temperature
        !!  @param pSIn equilibrium Pressure
    subroutine s_pSat(pSat, TSat, pSIn)
        !$acc routine seq

        real(kind(0.0d0)), intent(OUT) :: pSat
        real(kind(0.0d0)), intent(IN) :: TSat, pSIn
        real(kind(0.0d0)) :: dFdp, Fp, Om !< auxiliary variables

        ! Generic loop iterators
        integer :: ns

        if ((TSat == 0.0d0) .and. (pSIn == 0.0d0)) then

            ! assigning Saturation pressure
            pSat = 0.0d0

        else

            ! calculating initial estimate for pressure in the pSat procedure. I will also use this variable to
            ! iterate over the Newton's solver
            pSat = pSIn

            ! iteration counter
            ns = 0

            ! underrelaxation factor
            ! Om = 1.0d-3
            Om = 1.0d0
            Fp = 0.0d0
            do while ((DABS(Fp) > ptgalpha_eps) .or. (ns == 0))
                ! increasing counter
                ns = ns + 1

                ! resetting pSat if it goes negative
                if (pSat < 0.0d0) pSat = 1e5

                ! calculating residual
                Fp = (A + B/TSat + C*DLOG(TSat) + D*DLOG(pSat + ps_inf(lp)) - &
                        DLOG(pSat + ps_inf(vp)))*TSat*(cvs(vp)*(gs_min(vp) - 1))

                ! calculating the jacobian
                dFdp = &
                    (D/(pSat + ps_inf(lp)) - 1d0/(pSat + ps_inf(vp)))*TSat*(cvs(vp)*(gs_min(vp) - 1))

                ! updating saturation temperature
                pSat = pSat - Om*Fp/dFdp

                if (ns > 1000) then
                    print *, 'Iteration counter: ', ns,'pSIn:',pSIn,'pSat:', pSat,'TSat:', TSat
                    print *, 'Fp:', Fp, 'dFdp',dFdp
                    print *, 'A', A, 'B', B, 'C', C, 'D', D, '(pi1 - D*pi2)/(D-1)', (ps_inf(lp) - D*ps_inf(vp))/(D - 1d0)
                    Fp = A + B/TSat + C*DLOG(TSat) + D*DLOG((ps_inf(lp) - D*ps_inf(vp))/(D - 1d0) + ps_inf(lp)) - &
                            DLOG((ps_inf(lp) - D*ps_inf(vp))/(D - 1d0) + ps_inf(vp))
                    print *, 'p =',(ps_inf(lp) - D*ps_inf(vp))/(D - 1d0), ':, F(p)=', Fp
                    call s_mpi_abort('too large iteration number when calculating pSat with given TSat. Exiting ...')
                endif

            end do

            
        end if

    end subroutine s_pSat

    !>  This subroutine finalizes the phase change module
    subroutine s_finalize_relaxation_solver_module()
    end subroutine

#endif

end module m_phase_change
