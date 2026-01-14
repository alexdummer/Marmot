#include "Marmot/LinearViscoelasticCompressibleMooneyRivlin.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool LinearViscoelasticCompressibleMooneyRivlinRegistered = MarmotMaterialFiniteStrainFactory::
      registerMaterial< LinearViscoelasticCompressibleMooneyRivlin >( "LINEARVISCOELASTICCOMPRESSIBLEMOONEYRIVLIN" );

  } // namespace Registration
} // namespace Marmot::Materials
