#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/TsaiWuPlasticity.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      TsaiWuPlasticityIsRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< TsaiWuPlasticityModel >(
        "TSAIWUPLASTICITY" );

  } // namespace Registration

} // namespace Marmot::Materials
