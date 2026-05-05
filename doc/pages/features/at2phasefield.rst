AT2 Phase Field Fracture
========================

A linear elastic phase-field model for brittle fracture based on the AT2 formulation.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`E`
     - Young's modulus
   * - 1
     - :math:`\nu`
     - Poisson's ratio
   * - 2
     - :math:`G_c`
     - Critical fracture energy (energy release rate)
   * - 3
     - :math:`l`
     - Internal length scale

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**
   * - :math:`\mathcal{H}`
     - ``maxCrackDrivingForce``
     - Crack driving force history variable (maximum elastic strain energy density)
   * - :math:`\boldsymbol{\varepsilon}`
     - ``strain``
     - Total strain tensor (Voigt notation)

Theory
------

The AT2 model couples a linear elastic material with a phase-field variable
:math:`\varphi \in [0, 1]`, where :math:`\varphi = 0` represents intact material and
:math:`\varphi = 1` represents fully fractured material.

Degraded Stress
^^^^^^^^^^^^^^^

The stress is degraded by a quadratic degradation function :math:`g(\varphi)`:

.. math::
   \boldsymbol{\sigma} = g(\varphi)\, \mathbb{C} : \boldsymbol{\varepsilon}

with the isotropic elastic stiffness tensor :math:`\mathbb{C}` and the quadratic degradation function

.. math::
   g(\varphi) = (1 - \varphi)^2.

Phase-Field Evolution Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The phase-field variable :math:`\varphi` is governed by the balance equation

.. math::
   \varphi - l^2 \Delta\varphi = \frac{2l}{G_c}(1 - \varphi)\,\mathcal{H}(\boldsymbol{\varepsilon}),

where :math:`\mathcal{H}` is the crack driving force history variable, defined as the
maximum positive elastic strain energy density encountered up to the current time:

.. math::
   \mathcal{H}^{n+1} = \max\!\left( \mathcal{H}^{n},\, \psi(\boldsymbol{\varepsilon}^{n+1}) \right),
   \quad \psi(\boldsymbol{\varepsilon}) = \tfrac{1}{2}\boldsymbol{\varepsilon} : \mathbb{C} : \boldsymbol{\varepsilon}.

The irreversibility of crack growth is enforced through the history variable :math:`\mathcal{H}`,
which ensures the phase field can only grow monotonically.

Implementation
--------------

The model is implemented within the general gradient-enhanced hypoelastic framework
(:cpp:class:`MarmotMaterialGeneralGradientEnhancedHypoElastic`), treating the phase-field
variable :math:`\varphi` as the single nonlocal variable (:math:`\kappa_1 = \varphi`).
The gradient enhancement coefficient is constant: :math:`c = l^2`.

.. doxygenclass:: Marmot::Materials::AT2PhaseField
   :allow-dot-graphs:
