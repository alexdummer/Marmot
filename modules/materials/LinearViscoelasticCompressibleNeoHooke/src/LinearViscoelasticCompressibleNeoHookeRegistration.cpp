#include "Marmot/LinearViscoelasticCompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearViscoelasticCompressibleNeoHookeRegistered = MarmotMaterialFiniteStrainFactory::
      registerMaterial< LinearViscoelasticCompressibleNeoHooke >( "LINEARVISCOELASTICCOMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
