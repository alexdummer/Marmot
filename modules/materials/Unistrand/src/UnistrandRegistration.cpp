#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/Unistrand.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool UnistrandIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< Unistrand >(
      "UNISTRAND" );

  } // namespace Registration
} // namespace Marmot::Materials
