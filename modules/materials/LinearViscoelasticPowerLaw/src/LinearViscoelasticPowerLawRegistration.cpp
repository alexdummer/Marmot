#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                        // scope

  const static bool LinearViscoelasticPowerLawisRegistered = MarmotMaterialHypoElasticFactory::registerMaterial<
    LinearViscoelasticPowerLaw >( "LINEARVISCOELASTICPOWERLAW" );

} // namespace Marmot::Materials::Registration
