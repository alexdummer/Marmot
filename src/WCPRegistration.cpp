#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/WCP.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int WoodCreepPlasticityCode = 11930000 + 20;

    using namespace MarmotLibrary;

    const static bool WoodCreepPlasticityIsRegistered = MarmotMaterialFactory::
      registerMaterial( WoodCreepPlasticityCode,
                        "WOODCREEPPLASTICITY",
                        makeDefaultMarmotMaterialFactoryFunction< class WoodCreepPlasticity >() );

  } // namespace Registration
} // namespace Marmot::Materials
