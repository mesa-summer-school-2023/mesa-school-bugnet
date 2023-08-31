Minilab 2: Rotation and its effect on stellar pulsations
===================================

Now that we have MESA and GYRE running, we want to investigate the impact of rotation. All stars rotate - however, describing how rotation impacts the structure and evolution of a star is a complicated business. Rotation induces two major changes to stellar evolution calculations. First, it introduces a new mechanism for transporting angular momentum throughout the star, and second it introduces a new mechanism for transporting chemicals throughout a star. In this Minilab, we want to investigate how we can use MESA to incorporate rotation in our stellar models and how different implementations modify the structure, and hence, the asteroseismic signature of a typical red giant. We will follow two cases: 1) Where we impose a constant viscosity to approximate rigid rotation, and 2) where we use a physical approach to compute the efficiency of angular momentum (AM) transport, both for a rotation rate at the zero-age main-sequence (ZAMS) of 20\% the critical rotation rate. Both of these will require modifications to the standard inlist that we will follow below.


Exercise 0: Setup
--------

We start this minilab from a work directory that you can download `here
<https://github.com/mesa-summer-school-2023/mesa-school-bugnet/blob/main/work_directories/work_mini2.zip>`__.

Exercise 1: Uniform rotation
--------

Asteroseismology of red giants (RG) has taught us that AM takes place in stars, but the theory of AM is not able to explain all observations. Generally, AM transport is treated by a combination of advective and diffusive processes (and possibly mechanisms like internal gravity waves generated at the convective core boundary).
In MESA, the transport of AM is done in a fully diffusive way by solving the following equation with rotation profile :math:`\Omega(m)` (with :math:`m` the mass coordinate),

.. math::

    \begin{split}
    \left(\frac{{\partial \Omega}}{\partial t}\right)_m &= \frac{1}{i}\left( \frac{\partial }{\partial m} \right)_t \left[ (4 \pi r^2 \rho)^2 i \nu_{\rm AM} \left( \frac{\partial \Omega}{\partial m} \right)_t \right] \\
    &- \frac{ \Omega}{r} \left( \frac{\partial r }{\partial t} \right)_m \left(\frac{{\rm d} \ln i }{{\rm d} \ln r} \right)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(1)
    \end{split}

Here, :math:`i` is the specific moment of inertia of a shell at mass coordinate :math:`m`, and :math:`\rho` is the local density. The second term on the right-hand side accounts for the change in the radius of the star as it evolves, while the first term describes the actual transport of AM.
The efficiency is encompassed in the viscosity :math:`\nu_{\rm AM}(m)`. First, we will set this viscosity to a very high value, such that AM is redistributed almost instanteously, and thus the star rotates as a solid body.

We need to make some small additions to the inlist. Specifically, we need to tell MESA how to deal with rotation.

First, we need to set some flags in the ``&star_job`` section of the inlist so MESA knows that we want to activate rotation.

.. code-block:: console

  ! Rotation

  new_rotation_flag = .true.
  change_rotation_flag = .false.
  change_initial_rotation_flag = .true.


We tell MESA that we want to set a new rotation rate, but we only want to set a new rotation rate at the beginning of the run. Therefore, we set the ``change_initial_rotation_flag`` to  ``.true.``.

Instead of doing complicated maths to figure out what the appropriate rotation rate is in radians per second, we will utilise the functionality of MESA to set the rotation rate as a fraction of the critical rotation rate.
The (Keplerian) critical rotation rate that MESA uses is defined as

.. math::

  \Omega_{\rm crit} = \sqrt{\frac{\Gamma G M_\star}{R_\star^3}}

where :math:`M_\star` and :math:`R_\star` are the mass and radius, respectively. The factor :math:`\Gamma` takes into account the radiation pressure and is 1 when the star is well below the Eddington luminosity. Set an initial rotation rate by adding the following lines to the ``&star_job`` section.

.. code-block:: console

  near_zams_relax_omega_div_omega_crit = .true.

  set_omega_div_omega_crit = .false.
  set_initial_omega_div_omega_crit = .true.
  new_omega_div_omega_crit = 0.2d0

  num_steps_to_relax_rotation = 50

  change_D_omega_flag = .true.
  new_D_omega_flag = .true.


Now that we've set the options in the ``&star_job`` section, we need to set the options in the ``&controls`` section. In this section, we want to set an arbitrary viscosity that is constant throughout the whole star, and in time. This can be achieved by setting the following options:

.. code-block:: console

       set_uniform_am_nu_non_rot = .true.
       uniform_am_nu_non_rot = 1d20 ! In cm^2/s


.. note::

    The ``run_star_extras.f90`` file has already been modified in ``extras_finish_step`` to terminate when the model reaches :math:`\nu_{\rm max}=180\,\mu{\rm Hz}`, and to start writing profiles only on the RGB.

    .. code-block:: console

        if (s% nu_max < 250.) s% write_profiles_flag = .true.
        if (s% nu_max < 180.) extras_finish_step = terminate

Now, do ``./clean && ./mk`` and run ``./rn``.
Look at the rotation profile ``log_omega`` in the PGplot. Is the rotation indeed uniform?

In the next step, we will be passing the stellar profiles to GYRE. The following lines in the ``&controls`` tell MESA to output also a separate input file for GYRE along with the profiles (you cannot output pulse data without also outputting the profiles).

.. code-block:: console

    write_pulse_data_with_profile = .true.
    pulse_data_format = 'GYRE'


Now, we will make the changes to the GYRE inlist ``gyre_mix_minilab2.in``.
In the GYRE inlist, we set

.. code-block:: console

     Omega_rot_source = 'MODEL'

Using this option, GYRE will use the rotation profile of the MESA model to account for the effect of rotation on the stellar pulsations. Next, tell GYRE to use the last MESA model as input

.. code-block:: console

     file = 'xxx.data.GYRE'

and give a name for the output (summary) file

.. code-block:: console

     summary_file = 'xxx_const_visc.summary'

In this minilab, we will run GYRE stand-alone like you have done during Tuesday's lab. To run GYRE, use

.. code-block:: console

    $GYRE_DIR/bin/gyre gyre_mix_minilab2.in

You can comment out the :math:`\ell = 2` mode in the inlist. If you feel bold, you can try increasing ``freq_min``. 
To have a quick inspection of the GYRE summary file, we will use the online `MESA explorer <https://billwolf.space/mesa-explorer/>`__ designed by Bill Wolf. Upload your summary file, and plot ``n_pg`` (the radial order) vs. ``Re(freq)`` (the real part of the mode frequency). Plot the data in a scatter plot.


Exercise 2: Physical approach
--------

Now, we want to take a more physical approach and compute the viscosity from the six (magneto)hydrodynamical processes implemented in MESA that can induce turbulence (and thus transport angular momentum).
The physics regarding these (magneto)hydrodynamical processes is described in `Heger et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...528..368H/abstract>`__. For a short summary, have a look at `Mombarg (2023) <https://ui.adsabs.harvard.edu/abs/2023arXiv230617211M/abstract>`__.

First, we now disable using a uniform viscosity in the inlist

.. code-block:: console

       set_uniform_am_nu_non_rot = .false.

In MESA, each process can be turned on and off separately. To enable all of them without any additional scaling, set all diffusion coefficients equal to 1,

.. code-block:: console

        D_DSI_factor = 1 ! Dynamical shear instability
        D_SH_factor  = 1 ! Solberg-Høiland instability
        D_SSI_factor = 1 ! Secular shear instability
        D_ES_factor  = 1 ! Eddington-Sweet circulation
        D_GSF_factor = 1 ! Goldreich-Schubert-Fricke instability
        D_ST_factor  = 1 ! Spuit-Tayler dynamo

Run MESA again with this other way of AM transport. Do you see any changes in the rotation profile?

.. warning::

    Do not forgot to change the name of your output directory through ``log_directory`` in the ``&controls`` section!

In the mixing panel of ``PGstar``, you should also be able to see the predicted viscosity (or diffusion coefficient) for each of the six processes.
However, because we set ``am_D_mix_factor = 0`` in ``&controls``, we only study the effect of AM transport and not on the transport of chemical elements.

Run GYRE again at the same age (again, remember to provide a different name for the summary file!), and compare the pulsations. Upload also this summary file to `MESA explorer <https://billwolf.space/mesa-explorer/>`__ and toggle between the two. You can fix the ranges of the x- and y-axis to make it easier to see the differences. Assuming we start with a certain initial rotation frequency, the final rotation profile at the point where we stop the models will depend on our choice of the treatment for AM transport. 
Are the final rotation profiles for the two cases different enough to observe with asteroseismology?

.. admonition:: Bonus exercise 1

    Try computing also the modes for :math:`m=-1` and :math:`m=1` (:math:`\ell=1`).

.. admonition:: Bonus exercise 2

    Try enabling only a subset of the hydrodynamical processes and run again. Which ones are the most important?
