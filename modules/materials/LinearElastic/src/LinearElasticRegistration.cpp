#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Factory;

  const static bool LinearElasticIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< LinearElastic >(
    "LINEARELASTIC" );

} // namespace Marmot::Materials::Registration
