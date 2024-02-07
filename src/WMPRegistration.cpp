#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/WMP.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int WoodMultisurfacePlasticityCode = 11930000 + 20;

    using namespace MarmotLibrary;

    const static bool WoodMultisurfacePlasticityIsRegistered = MarmotMaterialFactory::
      registerMaterial( WoodMultisurfacePlasticityCode,
                        "WOODMULTISURFACEPLASTICITY",
                        makeDefaultMarmotMaterialFactoryFunction< class WoodMultisurfacePlasticity >() );

  } // namespace Registration
} // namespace Marmot::Materials
