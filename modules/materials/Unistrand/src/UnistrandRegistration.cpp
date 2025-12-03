#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/Unistrand.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int UnistrandCode = 11930000 + 27;

    using namespace MarmotLibrary;

    const static bool UnistrandIsRegistered = MarmotMaterialFactory::
      registerMaterial( UnistrandCode, "UNISTRAND", makeDefaultMarmotMaterialFactoryFunction< class Unistrand >() );

  } // namespace Registration
} // namespace Marmot::Materials
