A description of the experiments in the Distance Estimation paper

`flow_dist_moon`: Finding the closest distance for a non-convex moon-shaped set on the Flow system. Figure 8, section IV.A.

`twist_distance_half_sphere`: Distance betweens points starting in a sphere following the Twist system and unsafe points in a half-sphere. Figure 10a, Section IV.B.

`twist_distance_half_sphere`: Same scenario as `twist_distance_half_sphere`, but now using the L4 distance. Figure 10b, Section IV.B.

`twist_distance_half_sphere_csp`: Applying correlative sparsity to the `twist_distance_half_sphere` example.

`twist_distance_compare_csp`: Compare runtimes and objectives for twist L2 systems. Tables III and IV, section V.D.

`flow_dist_shape_test`: Distance between points on a translating shape (square) and the half-circle unsafe set. Figure 15, section VI.D.

`flow_dist_shape_rotate_test`: The square is now rotating at a constant angular velocity of 1 radian/second. This optimization problem is significantly harder than in `flow_dist_shape_test` due to the shape measure. Figure 16, section VI.D.


`flow_distance_test_w`: Flow-half-circle distance test where dynamics have a time-varying uncertainty process. Figure 17, section VII.A.

`flow_distance_test_l1.m`: Flow-half-circle distance test with an L1 distance, enforced by a polytopic lift. Figure 18, section VII.B.