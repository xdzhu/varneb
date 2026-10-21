"""Author-approved finite design; no calculation or imports with side effects."""
POLICY = dict(conditions=['zero_q1', 'plus1_q0'], sign='minus',
    held=['Qx', 'Qz', 'Qy', 'polar_proxy_z'], maximum_calls_per_condition=40,
    maximum_new_calls=80, gradient_tolerance_eV_A=.03, initial_trust_radius_A=.05,
    maximum_coordinate_norm_A=1., minimum_distance_A=1.6, maximum_deformation=.25,
    optimizer='safeguarded_bfgs', maximum_iterations=200, maximum_allocations=2,
    MPI_ranks=32, partition='hfacnormal01', walltime='02:00:00',
    pressure_GPa=0., physical_inputs_changed=False, automatic_extension=False,
    full_stability_or_TS_or_global_barrier_certified=False)
