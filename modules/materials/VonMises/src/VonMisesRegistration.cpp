#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/VonMises.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool VonMisesIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< VonMisesModel >(
    "VONMISES" );

} // namespace Marmot::Materials::Registration
