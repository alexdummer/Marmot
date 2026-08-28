#include "Marmot/FiniteStrainIsotropicBiotViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                        // scope

  const static bool FiniteStrainIsotropicBiotViscoelasticityRegistered = MarmotMaterialFiniteStrainFactory::
    registerMaterial< FiniteStrainIsotropicBiotViscoelasticity >( "FINITESTRAINISOTROPICBIOTVISCOELASTICITY" );

} // namespace Marmot::Materials::Registration
