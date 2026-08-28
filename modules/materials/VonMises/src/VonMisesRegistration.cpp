#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/VonMises.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                          // scope

    const static bool VonMisesIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< VonMisesModel >(
      "VONMISES" );

  } // namespace Registration
} // namespace Marmot::Materials
