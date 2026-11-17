# m_data_output.fpp

write the grid and initial condition data files into the "0" time step

## s_write_serial_data_files

write data files into D/ and p_all/pk/0/

1. write x_cb.dat into p_all/pk/0/, and y_cb.dat, z_cb.dat similarly
    - use variable x_cb, y_cb, z_cb

2. write q_cons_vfi.dat into p_all/pk/0/
    - use q_cons_vf(i)%sf

3. create D/ if needed, and write prim., cons. into D/

    - 1D: both prim and cons (prim only for 5-equation model)
    - 2D: only cons
    - 3D: only cons

## s_write_parallel_data_files

- `MPI` parallel for I/O
- write data into file restart_data/lustre_0

## s_initialize_data_output_module

1. create directory p_all/pk/0/ for serial data and restart_data/
2. point s_write_data_files to s_write_serial/parallel_data_files

## s_finalize_data_output_module (release s_write_data_files)
