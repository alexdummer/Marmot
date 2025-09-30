#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/ViscoelasticNeoHooke.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int ViscoelasticNeoHookeCode = 11930000 + 13;

    using namespace MarmotLibrary;

    const static bool ViscoelasticNeoHookeRegistered = MarmotMaterialFactory::
      registerMaterial( ViscoelasticNeoHookeCode,
                        "VISCOELASTICNEOHOOKE",
                        makeDefaultMarmotMaterialFactoryFunction< class ViscoelasticNeoHooke >() );

  } // namespace Registration
} // namespace Marmot::Materials
