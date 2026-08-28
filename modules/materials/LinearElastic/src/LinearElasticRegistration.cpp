#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                          // scope

    const static bool LinearElasticIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< LinearElastic >(
      "LINEARELASTIC" );

  } // namespace Registration
} // namespace Marmot::Materials
