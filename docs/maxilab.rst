Maxilab
===================================

Contents
--------

Under the influence of rotation, pulsation modes with the same spherical degree, $\ell$, (number of surface nodal lines), split into multiplets :math:`m = -\ell, -(\ell-1), ..., \ell-1, \ell` where the azimuthal order :math:`m` indicates the number of surface nodal lines intersecting the rotational axis. On top of that, the presence of a magnetic field in the core can introduce asymmetries in the splittings (difference in frequency) of these pulsation modes. 
In this Maxilab, we are going to infer the internal magnetic field of the red giant (RG) KIC11515377, observed with the NASA Kepler mission. We follow the methodology of Li et al. (2022, Nature).  The squared radial magnetic field averaged in the horizontal direction is inferred by,
.. math::
    \left< B_r^2\right> = \frac{\mu_0 \delta \omega_g (2 \pi \nu_{\rm max})^3}{\mathcal{I}},

where :math:`\mu_0 = 4\pi \cdot 10^{-6} \,{\rm kG\,cm\,A^{-1} }` is the magnetic permeability in vacuum, :math:`\delta \omega_g` is the observed frequency shift of the g modes, and :math:`\nu_{\rm max}` is the frequency of maximum power. The factor :math:`\mathcal{I}` in the denominator is defined as,
.. math::
    \mathcal{I} = \frac{\int \left(\frac{N}{r}\right)^3 \frac{dr}{\rho}}{\int \frac{N}{r}dr},

where :math:`N` is the Brunt-V\"ais\"al\"a frequency, :math:`\rho` is the density, and :math:`r` is the radial coordinate.


