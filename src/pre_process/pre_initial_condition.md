# m_initial_condition.fpp

## s_initialize_initial_condition_module

Allocate variables and set default values for initialization.

Variables:

- q_prim_vf:    (:) of scalar_field
- q_cons_vf:    (:) of scalar_field
- patch_id_fp   (:,:,:) of integer, the identity of patch for initial regions, each patch has its own geometry and physical conditions.

## s_generate_initial_condition

Setup initial conditions on the grid according to patch parameters. Finally initialize the initial condtions within each grid cell (conservative variables).

1. set primitive variables within each cell

    $3D:p>0; 2D:n>0; 1D:m>0$

    | #    | Function               | Dim.  |
    | ---: | :----:             | :---: |
    | 1    | s_line_segment 	    | 1     |
    | 2    | s_circle 		    | 2     |
    | 3    | s_rectangle 	        | 2     |
    | 4    | s_sweep_line 		| 2     |
    | 5    | s_ellipse 		    | 2     |
    | 6    | N/A (no defination)		    | 2     |
    | 7    | s_2D_analytical 	    | 2     |
    | 8    | s_sphere 		    | 3     |
    | 9    | s_cuboid 		    | 3     |
    | 10   | s_cylinder 		    | 3     |
    | 11   | s_sweep_plane 	    | 3     |
    | 12   | s_ellipsoid 		    | 3     |
    | 13   | s_3D_analytical 	    | 3     |
    | 14   | s_spherical_harmonic | 3     |
    | 15   | s_1d_analytical      | 1     |
    | 16   | s_1d_bubble_pulse    | 1     |
    | 17   | s_spiral             | 2     |
    | 18   | s_varcircle       | 2     |
    | 19   | s_3dvarcircle       | 3     |
    | 20   | s_2D_TaylorGreen_vortex  | 2  |
    | 21   | s_model              | 2 & 3 |

2. adding perturb to the field if any

2. convert primitive variables to conservative ones

## s_perturb_sphere, ...

...

## s_finalize_initial_condition_module

deallocate variables

Variables:

- q_prim_vf
- q_cons_vf
- patch_id_fp