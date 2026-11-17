# m_patches.fpp

## s_line_segment (1D)

Variables:

- patch_id
- q_prim_vf

User-defined geometric parameters:

- x_centroid
- length_x

1. load geometric parameters

- x_centroid
- length_x
- x_boundary%beg: x_centroid - length/2
- x_boundary%end: x_centroid + length/2

2. assign primitive variables within the region of corresponding patch

    - consider smoothing around the patch boundary

    eta = $\eta = 1$ (no smoothing)

    - s_assign_patch_primitive_variables

## s_spiral (2D)

User-defined geometric parameters:

- x_centroid, y_centroid
- radius
- length_x, length_y

...

## s_2D_TaylorGreen_Vortex (2D)

User-defined geometric parameters:

- x_centroid, y_centroid
- length_x, length_y

1. load geometric parameters

- x_centroid
- length_x
- x_boundary%beg: x_centroid - length/2
- x_boundary%end: x_centroid + length/2
- $U_0 = u$
- $L_0 = v$

2. assign primitive variables within the region of corresponding patch

    - consider smoothing around the patch boundary

    eta = $\eta = 1$ (no smoothing)

    - s_assign_patch_primitive_variables
    - @:analytical()
    - assign other specific parameters

        $$
        u = U_0\sin(\frac{x_i}{L_0})\cos(\frac{y_j}{L_0}) \\
        v = -U_0\cos(\frac{x_i}{L_0})\sin(\frac{y_j}{L_0}) \\
        p = p + (\frac{\cos(2x_j)}{L_0} + \frac{\cos(2y_j)}{L_0})\frac{\alpha_1\rho_1 U_0^2}{16}
        $$

## s_1D_analytical (1D)

self-defined analytical functions for initialization

1. load initial parameters

   - @:Hardcoded1DVariables()
   - x_centroid
   - length_x
   - x_boundary%beg: x_centroid - length/2
   - x_boundary%end: x_centroid + length/2
   - $\eta = 1$ (no smoothing)

2. assign primitive variables within region of each patch

    - s_assign_patch_primitive_variables
    - @:hardcoded1D()

## s_2D_analytical (2D)

self-defined 2D analytical functions for initialization

1. load initial parameters

   - @:Hardcoded2DVariables()
   - x_centroid, y_centroid
   - length_x, length_y
   - x_boundary%beg: x_centroid - length_x/2
   - x_boundary%end: x_centroid + length_x/2
   - y_boundary%beg: y_centroid - length_y/2
   - y_boundary%end: y_centroid + length_y/2
   - $\eta = 1$ (no smoothing)

2. assign primitive variables within region of each patch

    - s_assign_patch_primitive_variables
    - @:hardcoded2D()

## s_convert_cylindrical_to_cartesian_coord

$$
cart_x = cyl_x \\
cart_y = cyl_y\sin(cyl_z) \\
cart_z = cyl_y\cos(cyl_z)
$$
