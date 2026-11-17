# m_start_up.fpp

## s_read_input_file

Reads the configuration file pre_process.inp

Variables:

- file_loc
- file_check
- iostatus
- line

## s_check_input_file

Check whether the user inputs make sense.

1. check whether the case file exist or not
2. s_check_inputs
3. s_check_patches

## s_read_serial_grid_data_files

Read any preexisting grid data

Variables:

- $x_{cb}$: cell-boundary location
- $x_{cc}$: cell-center location
- x_domain%beg = $x_{cb}(-1)$
- x_domain%end = $x_{cb}(m)$

1. check existence of files `x_cb.dat` in `p_all/p0`
2. read $x_{cb}$
3. compute $x_{cc}$

$$
x_{cc}(i) = \frac{x_{cb}(i) + x_{cb}(i-1)}{2} \\
dx = x_{cb}(0) - x_{cb}(-1)
$$

4. $y_{cc},z_{cc}$ similar to $x_{cc}$

## s_check_grid_data_files

check positivity of cell-widths

## s_read_serial_ic_data_files

1. check existence of files `q_cons_vf.dat` in `p_all/p0`
2. read `q_cons_vf(i)%sf`

## s_read_parallel_grid_data_files

Read any preexisting grid data, parallel using `MPI`

Variables:

- $x_{cb}$: cell-boundary location
- $x_{cc}$: cell-center location
- x_domain%beg = $x_{cb}(-1)$
- x_domain%end = $x_{cb}(m)$

1. check existence of files `x_cb.dat` in `reastart_data/0`
2. read $x_{cb}$
3. compute $x_{cc}$

$$
x_{cc}(i) = \frac{x_{cb}(i) + x_{cb}(i-1)}{2} \\
dx = x_{cb}(0) - x_{cb}(-1)
$$

4. $y_{cc},z_{cc}$ similar to $x_{cc}$

## s_read_parallel_ic_data_files

`MPI` data I/O

1. check existence of files `q_cons_vf.dat` in `restart_data/0`
2. read `q_cons_vf(i)%sf`

## s_initialize_modules

initialize many modules

create output directory for serial I/O 'D'

## s_apply_initial_condition

1. s_generate_initial_condition (achieve conservative variables within each cell)
2. s_infinite_relaxation_k (consider `phase-change` if required)
3. s_write_data_files (restore the conservative variables in 'D')

## s_save_data (save time data)

## s_initialize_mpi_domain

## s_finalize_modules

1. disassociate pointers for serial and parallel I/O
2. deallocation procedures for the modules
3. finalize MPI environment