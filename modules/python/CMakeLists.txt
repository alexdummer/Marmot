# Fetch nanobind
include(FetchContent)
FetchContent_Declare(
  nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind.git
  GIT_TAG        v1.9.2
)
FetchContent_MakeAvailable(nanobind)

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

# Define the python extension module
nanobind_add_module(marmot
    src/bind_main.cpp
    src/bind_materials.cpp
    src/bind_solver_finite_strain.cpp
    src/bind_solver_hypo_elastic.cpp
)

# Link against the Marmot C++ library
target_link_libraries(marmot PRIVATE Marmot)

# Set RPATH so it finds libMarmot.so in the build tree instead of conda env
set_target_properties(marmot PROPERTIES
    BUILD_RPATH "${CMAKE_BINARY_DIR}/lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
)

# Add include directory for fastor_nanobind_caster.h
target_include_directories(marmot PRIVATE src)
