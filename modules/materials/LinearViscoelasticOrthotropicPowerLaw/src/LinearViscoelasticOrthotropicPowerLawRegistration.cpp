#include "Marmot/LinearViscoelasticOrthotropicPowerLaw.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                        // scope

  const static bool LinearViscoelasticOrthotropicPowerLawisRegistered = MarmotMaterialHypoElasticFactory::
    registerMaterial< LinearViscoelasticOrthotropicPowerLaw >( "LINEARVISCOELASTICORTHOTROPICPOWERLAW" );

} // namespace Marmot::Materials::Registration
