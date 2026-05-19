Compressible Neo Hooke model (Automatic Differentiation)
=========================================================

An automatic differentiation implementation of the compressible Neo-Hookean hyperelastic
material model using the Pence–Gou potential (variant B).

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`K`
     - Bulk modulus
   * - 1
     - :math:`G`
     - Shear modulus
   * - 2 *(optional)*
     - :math:`\rho`
     - Density

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Description**
   * - —
     - No state variables required.

Theory
------

See :ref:`compressibleneohooke`.

Implementation
--------------

Implementation follows the :ref:`compressibleneohooke`, however, automatic differentiation
is used to compute the consistent algorithmic tangent operator.

.. doxygenclass:: Marmot::Materials::ADCompressibleNeoHooke
   :allow-dot-graphs:
