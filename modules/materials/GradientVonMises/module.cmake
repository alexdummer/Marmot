set(MODULE_NAME
    "GradientVonMises")

set(MODULES_DEPENDENCIES
    MarmotMechanicsCore
    MarmotGradientMechanicsCore
    )

set(DEPENDENCIES_FULFILLED TRUE)
foreach( DEPENDENCY ${MODULES_DEPENDENCIES} )
    if (NOT DEPENDENCY IN_LIST INSTALLED_MODULES)
        message("----> " "module ${MODULE_NAME} dependency not fulfilled: ${DEPENDENCY}")
        set(DEPENDENCIES_FULFILLED FALSE)
    endif()
endforeach()

if ( DEPENDENCIES_FULFILLED )
    include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
    file(GLOB sources_material "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
    list(APPEND sources ${sources_material})
endif()
