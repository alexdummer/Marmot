#include "Marmot/CompressibleFiniteStrainLinearViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainSubstepped.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                        // scope

  const static bool CompressibleFiniteStrainLinearViscoelasticityRegistered = MarmotMaterialFiniteStrainFactory::
    registerMaterial< CompressibleFiniteStrainLinearViscoelasticity >(
      "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY" );

  const static bool
    CompressibleFiniteStrainLinearViscoelasticitySubsteppedRegistered = MarmotMaterialFiniteStrainFactory::
      registerMaterial< MarmotMaterialFiniteStrainSubstepped< CompressibleFiniteStrainLinearViscoelasticity > >(
        "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY_SUBSTEPPED" );
} // namespace Marmot::Materials::Registration
