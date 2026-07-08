#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_materials( nb::module_& m )
{
  // Basic material bindings could go here
  // Currently, solvers instantiate materials internally via factories
  // using the material name and properties array.
}
