#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_materials( nb::module_& m );
void bind_finite_strain_solver( nb::module_& m );
void bind_hypo_elastic_solver( nb::module_& m );

NB_MODULE( _marmot, m )
{
  m.doc() = "Marmot Python bindings";

  nb::module_ core = m.def_submodule( "core", "Core types and functionality" );

  nb::module_ solvers = m.def_submodule( "solvers", "Material point solvers" );
  bind_finite_strain_solver( solvers );
  bind_hypo_elastic_solver( solvers );

  nb::module_ materials = m.def_submodule( "materials", "Material models and factories" );
  bind_materials( materials );
}
