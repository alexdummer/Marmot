#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/WoodCreepDamagePlasticity.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int WoodCreepDamagePlasticityCode = 11930000 + 20;

    using namespace MarmotLibrary;

    const static bool WoodCreepDamagePlasticityIsRegistered = MarmotMaterialFactory::
      registerMaterial( WoodCreepDamagePlasticityCode,
                        "WOODCREEPDAMAGEPLASTICITY",
                        makeDefaultMarmotMaterialFactoryFunction< class WoodCreepDamagePlasticity >() );

  } // namespace Registration
} // namespace Marmot::Materials
