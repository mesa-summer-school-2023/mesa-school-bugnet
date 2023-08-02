Maxilab: Measuring the magnetic field
===================================

Contents
--------

Under the influence of rotation, pulsation modes with the same spherical degree, :math:`\ell`, (number of surface nodal lines), split into multiplets :math:`m = -\ell, -(\ell-1), ..., \ell-1, \ell` where the azimuthal order :math:`m` indicates the number of surface nodal lines intersecting the rotational axis. On top of that, the presence of a magnetic field in the core can introduce asymmetries in the splittings (difference in frequency) of these pulsation modes. 
In this Maxilab, we are going to infer the internal magnetic field of the red giant (RG) KIC11515377, observed with the NASA Kepler mission. We follow the methodology of Li et al. (2022, Nature).  The squared radial magnetic field averaged in the horizontal direction is inferred by,

.. math::

    \left< B_r^2\right> = \frac{\mu_0 \delta \omega_g (2 \pi \nu_{\rm max})^3}{\mathcal{I}},~~~~~~~~~~~(1)

where :math:`\mu_0 = 4\pi \cdot 10^{-6} \,{\rm kG\,cm\,A^{-1} }` is the magnetic permeability in vacuum, :math:`\delta \omega_g` is the observed frequency shift of the g modes, and :math:`\nu_{\rm max}` is the frequency of maximum power. The factor :math:`\mathcal{I}` in the denominator is defined as,

.. math::

    \mathcal{I} = \frac{\int \left(\frac{N}{r}\right)^3 \frac{dr}{\rho}}{\int \frac{N}{r}dr},~~~~~~~~~~~(2)

where :math:`N` is the Brunt-V\"ais\"al\"a frequency, :math:`\rho` is the density, and :math:`r` is the radial coordinate.

The denominator in the integral relates to the asymptotic period spacing of modes with spherical degree :math:`\ell = 1` as follows,

.. math::

    \Delta \Pi_{\ell = 1} = \frac{2 \pi^2}{\sqrt{2}}\left( \int \frac{N}{r}dr \right)^{-1}.~~~~~~~~~~~(3)

We are going to compute the quantity :math:`\mathcal{I}` of a MESA model in the ``run_star_extras.f90``. 

Exercise 0 
--------
As a first step, copy the work directory over and verify it runs. 

At each step of the evolution we want to compute :math:`\left< B_r^2\right>^{1/2}` and store it in the output of the history file. 

Exercise 1 
--------
Prepare in the ``run_star_extras.f90`` three additional history columns named ``I``, ``Br_mean``, and ``Delta_Pi1``. Set the correct number of additional columns in the ``how_many_extra_history_columns`` function, and add the following to ``data_for_extra_history_columns`` for each additional column,

.. code-block:: console

    names(1) = '...'
    vals(1) = ...

You can set the values to 0 for now.
Do a ``./clean`` and ``./mk`` and check this works.

Exercise 2
--------
The first step is to compute the two integrals in Eq (2). For the Brunt-Väisälä frequency, we need to first ensure it is zero in convective regions and so we compute a new array with all elements >0. A new array of a variable length is defined as follows,

.. code:: fortran

    double precision, allocatable :: brunt_N(:)

    allocate(brunt_N(s% nz))

This defines an array with the same length as the number of cells at each time step.
Here, the declaration of the (double precision) variable goes right below the ``subroutine`` statement, and the allocate statement after all other variable declarations and the call to the ``star_info`` structure has been made. These are the lines

.. code:: fortran

    call star_ptr(id,s,ierr)
    if(ierr/=0) return

We can then access variables part of the ``star_info`` structure such as the radius, density, and the squared Brunt-Väisälä frequency (:math:`N^2`)

.. code:: console

    s% r
    s% rho
    s% brunt_N2

You can check out ``MESA_DIR/star_data/public/star_data_work.inc`` to see what variables are accessible this way.
Moreover, ``s\% r(k)`` will give you the k-th element of the array.

Compute :math:`N` from the values of :math:`N^2` defined in MESA, but set negative values to zero.

.. code:: console

    sqrt(max(0._dp, s% brunt_N2))

In Fortran, the function ``max()`` will element-wise return the larger element of the two arguments. The ``_dp`` indicates we are dealing with double precision here.
At the end of the subroutine, you can deallocate the array to free up memory.

.. code:: console

    deallocate(brunt_N)

If your model has a high enough spatial resolution, you can assume,

.. math::

    \int x\,{\rm d}x \approx \sum_i x_i\,\Delta x_i,

where the index :math:`i` runs over the cells.
First, define two quantities in which you store the values of the two integrals. For the summation (integral), you will have to something like

.. code:: fortran

    sum = 0._dp
    do k = 1, s% nz-1
      sum = sum + delta(k)
    end do

where ``delta(k)`` is the function we want to integrate (:math:`x_i \Delta x_i`). Remember :math:`k=1` is the outermost cell.
In MESA, there are quantities that are defined at the mass centre of the cell, and there are quantities that are defined at the edge of the cell. Think about this when you compute the integrals.

.. tip::

   In ``star_info``, ``s% r`` is defined at the cell edge, while ``s% rmid`` is defined at the centre.


Once you have computed :math:`\mathcal{I}`, write this value out to the first extra column in history.

Exercise 3
--------
Next, we want to pass on the value of :math:`\delta \omega_g` to the ``run_star_extras.f90``. In your inlist, you can set

.. code:: console

    x_ctrl(1) = ...

to a value that you can then access in the ``run_star_extras.f90`` through,

.. code:: console

    s% x_ctrl(1)

Add a control in your inlist to do this. The observed value for KIC11515377 is :math:`\delta \omega_g / (2 \pi) = 126` nHz. The value of math:`\nu_{\rm max}` you can get from ``star_info``. Pay attention to the correct units. In ``MESA_DIR/star_data/public/star_data_work.inc`` you can also find the units of each quantity in ``star_info``. Unless specified, MESA works in cgs units.

Finally, write :math:`\left< B_r^2\right>^{1/2}` and :math:`\Delta \Pi_1` also to your history file. Recompile and verify that on the RGB you find an average magnetic field of the order of 100 kG.

Exercise 4
--------
Finally, we want to stop the evolution when the model has roughly reached the observed values of :math:`\nu_{\rm max, obs} = 191.6 \pm 1\,\mu{\rm Hz}` and :math:`\Delta \Pi_{\rm 1, obs} = 83.16 \pm 1\,{\rm s}`. Add two additional controls to your inlist to pass these two values on to the ``run_star_extras.f90`` and define

.. math::

   \chi^2 = (\nu_{\rm max} - \nu_{\rm max, obs})^2 + (\Delta \Pi_1 - \Delta \Pi_{1, \rm obs})^2.

Change the inlist to start the evolution from the zero-age main sequence instead of loading in a precomputed RGB model. Be sure to properly set the initial composition by setting

.. code:: console

    set_uniform_initial_composition = .true.

Once on the RGB, after each time step, check whether the :math:`\chi^2` is smaller or bigger than the previous value. If it is bigger, terminate. First, define a global variable in which you store the value of :math:`\chi^2`. A global variable means this variable can be accessed by all subroutines in the ``run_star_extras.f90``, and is declared at the start of the ``run_star_extras.f90``, right below ``implicit none``. Now, in ``data_for_extra_history_columns`` you can set the value of :math:`\chi^2`.
In addition, also define a global variable which stores the previous value of :math:`\chi^2`. For the first time step, we need to initialise this variable to a large value (e.g. 1e99).

.. code:: console

    chi2_old = 1d99

Have a look at the flowchart in Fig.~\ref{fig:flowchart} and see which subroutine is only called once at the start of a run.
Lastly, check in the flowchart where MESA decides to keep going or terminate. Here, add a condition that will terminate the run if the new :math:`\chi^2` is larger than the previous value. Else, update the previous value to the new one. To make sure we are on the RG branch, add the following second condition

.. code:: console

    safe_log10(s% Teff) < 3.7

Add to your PGstar inlist the target values, so that you can see how close your models gets to the observations. To do this, have a look at the controls in ``inlist_pgstar`` that are currently commented out.
Pick a value for the initial mass from the spreadsheet and note down the lowest found :math:`\chi^2` value and the corresponding value of the internal magnetic field (in kG).

.. image:: flowchart_run_star_extras.png
   :alt: Flowchart
   :scale: 100%
   :align: left

   Flowchart of the ``run_star_extras``, taken from the MESA docs.

.. admonition:: Solution

    The part where you compute and add the additional history columns should look something like this.

    .. code:: fortran

        subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
           integer, intent(in) :: id, n
           character (len=maxlen_history_column_name) :: names(n)
           real(dp) :: vals(n), integral_N, integral_N3, I, mu_0, Br_mean
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           double precision, allocatable :: brunt_N(:)
           integer :: k
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           mu_0 = 4d-6*pi

           ! note: do NOT add the extras names to history_columns.list
           ! the history_columns.list is only for the built-in history column options.
           ! it must not include the new column names you are adding here.

           allocate(brunt_N(s% nz))
           names(1) = 'I'
           names(2) = 'Br_mean'
           names(3) = 'Delta_Pi1'
           brunt_N = sqrt(max(0._dp,s% brunt_N2))
           integral_N3 = 0.0_dp
           integral_N = 0.0_dp
           do k = 1, s%nz-1
             integral_N3 = integral_N3 + (brunt_N(k)**3/(s% rho(k)))*abs(s% rmid(k+1) - s% rmid(k)) / (s% r(k))**3
             integral_N  = integral_N + brunt_N(k)*abs(s% rmid(k+1) - s% rmid(k)) / s% r(k)
           end do
           I = integral_N3 / integral_N
           vals(1) = I
           omega_max = 2 * pi * s% nu_max * 1d-6
           Br_mean = sqrt(mu_0 * (2*pi*delta_omega_g*1d-9) * omega_max**3 / I) ! In kG.
           vals(2) = Br_mean
           Delta_Pi1 = (2._dp*pi**2)/integral_N / (sqrt(2._dp))
           vals(3) = Delta_Pi1
           write(*,*) 'Br_mean [kG] = ', Br_mean, 'Delta_Pi1 [s] = ', Delta_Pi1, 'nu_max [uHz] = ', s% nu_max, 'delta_nu [uHz]', s% delta_nu,   'I = ', I
           chi2 = (Delta_Pi1 - s% x_ctrl(2))**2 + (s% nu_max - s% x_ctrl(3))**2
           write(*,*) 'chi2', chi2
           deallocate(brunt_N)

        end subroutine data_for_extra_history_columns
