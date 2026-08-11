list(APPEND INSTALLED_MODULE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
file(GLOB module_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
list(APPEND sources ${module_sources})
list(APPEND publicheaders
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialGradientPlasticityFiniteStrain.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialGradientPlasticityFiniteStrainAD.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialGradientPlasticityFiniteStrainFactory.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialGradientPlasticityHypoElastic.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/Marmot/MarmotMaterialGradientPlasticityHypoElasticFactory.h")
