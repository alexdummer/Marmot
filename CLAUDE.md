# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Marmot (MAteRialMOdellingToolbox) is a C++20 shared library providing finite elements and constitutive
material models for quasi-brittle materials (concrete, rock, soil). It's built as a single shared
library (`libMarmot`) assembled at CMake-configure time from independently-pluggable modules, and is
consumed by external FE codes (EdelweissFE, Abaqus, MOOSE, OpenSees) through small, code-agnostic
interfaces (`MarmotMaterial*`, `MarmotElement`).

Dependencies: Eigen (>3.3.8), autodiff (>0.6.0), Fastor (>6.4.0). All discovered via
`find_package`/`find_path` in the top-level `CMakeLists.txt`; there is no vendoring.

## Build & test

```bash
mkdir build && cd build
cmake ..
make install                       # builds libMarmot and installs it (CMAKE_INSTALL_PREFIX)
ctest --output-on-failure          # run all registered tests
ctest -R TestLinearElastic         # run a single test by name
ctest --output-on-failure --no-tests=error   # CI form: fail if zero tests were registered
```

In-source builds (`cmake .` from repo root) are explicitly rejected — always build from `build/`.
Default `CMAKE_BUILD_TYPE` is `Release`.

Formatting/static checks run via pre-commit (clang-format on C++, black/isort/flake8/autoflake on
Python):

```bash
pre-commit install                 # one-time
pre-commit run --all-files
```

`scripts/projectmanager.py <cmd...>` runs an arbitrary command inside the repo root and inside every
directory that contains a `module.cmake` (3 levels under `modules/`) — useful for repo-wide,
per-module operations (e.g. `git`, linting).

## Module architecture

`modules/` is organized into categories, each scanned independently by the top-level `CMakeLists.txt`:
`core`, `materials`, `elements`, `particles`, `materialpoints`, `cells`, `cellelements`. Every
immediate subdirectory containing a `module.cmake` is auto-discovered and pulled into the build —
**no top-level CMakeLists.txt edits are needed to add a new material or element.** Filter variables
(`CORE_MODULES`, `MATERIAL_MODULES`, `ELEMENT_MODULES`, ...) can restrict which modules get built;
they default to `"all"`.

A module's `module.cmake` conventionally does three things:
```cmake
list(APPEND INSTALLED_MODULE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
file(GLOB module_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${module_sources})
```
and an optional `test.cmake` (auto-included if present) registers its tests via the
`add_marmot_test(<name> <source>)` helper defined in the top-level `CMakeLists.txt`. Categories are
scanned in a fixed order (core → materials → elements → particles → materialpoints → cells →
cellelements), so a module can only depend on modules from an earlier category. Modules that need a
specific other module present (elements depending on specific material/core modules) guard their
`module.cmake` body with a `MODULE_DEPENDENCIES` check against the already-populated `INSTALLED_MODULES`
list — see `modules/elements/DisplacementFiniteElement/module.cmake` for the pattern. A module may also
append to `publicheaders` to control which of its headers get installed (see
`modules/core/MarmotMechanicsCore/module.cmake`).

Standard module layout:
```
modules/<category>/<ModuleName>/
  module.cmake
  test.cmake                # optional
  include/Marmot/*.h        # public headers, installed flat under Marmot/
  src/*.cpp
  test/test.cpp             # single main(), uses MarmotTesting.h helpers
  README.md                 # optional
```

### Self-registering factories

Materials and elements are looked up by string name at runtime through factory singletons in the
`MarmotLibrary` namespace — there is no central registry file to edit. Each module registers itself
via a static-initialization trick in a dedicated `<ModuleName>Registration.cpp`:

```cpp
#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {
  namespace Registration {
    using namespace MarmotLibrary;
    const static bool LinearElasticIsRegistered =
      MarmotMaterialHypoElasticFactory::registerMaterial< LinearElastic >( "LINEARELASTIC" );
  } // namespace Registration
}
```

There are several parallel factories depending on the base class/kinematics a module implements —
pick the one matching the base class you derive from:
- `MarmotMaterialHypoElasticFactory` — small-strain materials deriving `MarmotMaterialHypoElastic`
- `MarmotMaterialFiniteStrainFactory` — finite-strain materials
- `MarmotMaterialGeneralGradientEnhancedHypoElasticFactory` — gradient-enhanced (nonlocal) materials
- `MarmotElementFactory` — finite elements deriving `MarmotElement`

When adding a new material/element, follow the existing sibling module exactly (header in `include/`,
implementation in `src/*.cpp`, registration in `src/*Registration.cpp`, test in `test/test.cpp` plus
`test.cmake`) rather than inventing a new structure.

### Core module responsibilities (`modules/core/`)

- `MarmotUtilitiesCore` — logging (`MarmotJournal`), typedefs, `MarmotTesting.h` test helpers
- `MarmotMathCore` — general math utilities
- `MarmotMechanicsCore` — small-strain (hypoelastic) material base classes, elasticity/viscoelasticity
  helpers, Voigt notation, Haigh-Westergaard invariants, stress return-mapping substeppers
  (`PerezFouget*`, `AdaptiveSubstepper*`), yield surface utilities
- `MarmotFiniteStrainMechanicsCore` — finite-strain material base classes and kinematics
- `MarmotGradientMechanicsCore` — gradient-enhanced (nonlocal damage) material base classes
- `MarmotFiniteElementCore` — `MarmotElement` base class, shape functions, geometry, DOF layout,
  enhanced-assumed-strain utilities

`particles`, `materialpoints`, `cellelements` exist as empty category directories under `modules/`
(registered in the build but with no modules yet). `cells` is registered as a category in the
top-level `CMakeLists.txt` but `modules/cells/` doesn't exist at all — `register_module_category`
silently no-ops when the path is missing.

## Testing conventions

- Each test is a standalone executable with one `main()` in `modules/<category>/<Module>/test/test.cpp`
  that runs all checks for that module and exits nonzero on failure.
- Use helpers from `Marmot/MarmotTesting.h` rather than hand-rolling comparisons.
- Register new tests in the module's `test.cmake` via `add_marmot_test("TestName" "${CURR_TEST_SOURCE_DIR}/....cpp")`.
- Tests must be deterministic and fast (they run under plain `ctest`, no fixtures/harness beyond that).

## Commit & PR conventions (from CONTRIBUTING.md)

- Commits **must** follow [Conventional Commits](https://www.conventionalcommits.org)
  (`feat`, `fix`, `docs`, `test`, `refactor`, `perf`, `style`, `build`, `ci`, `chore`, `revert`), with
  an optional scope, e.g. `fix(mechanicscore): guard haighWestergaard() rho derivative singularity`.
- Feature PRs must include documentation updates and ctest-registered tests for new behavior.
- Doxygen comments (`@brief`, `@param`, `@return`, `@tparam`, `@details`) are expected on all public
  classes/functions/members — this codebase documents heavily, including LaTeX math in comments; match
  that style when touching public headers.
- New material models / elements require Sphinx documentation (`module.rst` under
  `doc/pages/features/...`) referenced from `materials.rst`/`elements.rst` — check `doc/` for the
  registration pattern when adding one.
