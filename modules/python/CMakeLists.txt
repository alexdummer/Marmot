# Find Python first, so that the very same interpreter is used to patch nanobind
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

# Fetch nanobind
include(FetchContent)
FetchContent_Declare(
  nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind.git
  GIT_TAG        v1.9.2
  PATCH_COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/patch_nanobind.py" "include/nanobind/nb_attr.h"
)
FetchContent_MakeAvailable(nanobind)

# Define the python extension module
nanobind_add_module(_marmot
    src/bind_main.cpp
    src/bind_materials.cpp
    src/bind_solver_finite_strain.cpp
    src/bind_solver_hypo_elastic.cpp
)

# Link against the Marmot C++ library
target_link_libraries(_marmot PRIVATE Marmot)

# Set RPATH so it finds libMarmot.so in the build tree instead of conda env
set_target_properties(_marmot PROPERTIES
    BUILD_RPATH "${CMAKE_BINARY_DIR}/lib"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
)

# Add include directory for fastor_nanobind_caster.h
target_include_directories(_marmot PRIVATE src)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import sysconfig, os; print(os.path.relpath(sysconfig.get_path('platlib'), sysconfig.get_path('data')))"
    OUTPUT_VARIABLE PYTHON_RELATIVE_INSTALL_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

install(TARGETS _marmot DESTINATION ${PYTHON_RELATIVE_INSTALL_DIR}/marmot)
install(DIRECTORY marmot/ DESTINATION ${PYTHON_RELATIVE_INSTALL_DIR}/marmot)



