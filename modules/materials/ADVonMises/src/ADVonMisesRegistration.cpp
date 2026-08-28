#include "Marmot/ADVonMises.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                          // scope

    const static bool ADVonMisesIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< ADVonMises >(
      "ADVONMISES" );

  } // namespace Registration
} // namespace Marmot::Materials
