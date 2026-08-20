.. _incompressibleneohooke:

Incompressible Neo-Hooke model
===============================

Theory
------

The classical incompressible (purely isochoric) Neo-Hookean model has energy density

.. math::

   \Psi_{\rm iso}(\bar\lambda_1,\bar\lambda_2,\bar\lambda_3)
   =
   \frac{\mu}{2} \left( \bar\lambda_1^2 + \bar\lambda_2^2 + \bar\lambda_3^2 - 3 \right)
   =
   \frac{\mu}{2} \left( \bar I_1 - 3 \right),

where :math:`\bar\lambda_i = J^{-1/3}\lambda_i` are the isochoric principal stretches
(:math:`J=\det\boldsymbol F`) and :math:`\bar I_1` is the first invariant of the isochoric right
Cauchy-Green tensor. Unlike :ref:`compressibleneohooke`, there is no volumetric energy term
:math:`U(J)`: the material resists no volumetric deformation on its own and must be paired with a
mixed displacement-pressure element, e.g. :ref:`displacementpressurefinitestrainelement`, that
enforces incompressibility as an independent constraint.

This is exactly the one-term :ref:`ogden` model with exponent :math:`\alpha=2` and modulus
:math:`\mu_1=\mu`, and is implemented as such: rather than re-deriving the same principal-stretch
stress and tangent, this material reuses Marmot::Materials::Ogden's validated spectral machinery
internally, exposing only the familiar single-parameter Neo-Hookean interface.

Implementation
--------------

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`\mu`
     - Shear modulus
   * - 1 (optional)
     - :math:`\rho`
     - Density

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**

   * - ---
     - ---
     - No state variables required.

.. doxygenclass:: Marmot::Materials::IncompressibleNeoHooke
   :allow-dot-graphs:
