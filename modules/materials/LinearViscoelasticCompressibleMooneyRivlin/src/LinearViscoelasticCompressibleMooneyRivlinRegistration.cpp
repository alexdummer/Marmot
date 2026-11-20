#include "Marmot/LinearViscoelasticCompressibleMooneyRivlin.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace {
  using namespace Marmot::Materials;

  const static bool LinearViscoelasticCompressibleMooneyRivlinRegistered = Marmot::Registry::
    MaterialFiniteStrainFactory::registerMaterial< LinearViscoelasticCompressibleMooneyRivlin >(
      "LINEARVISCOELASTICCOMPRESSIBLEMOONEYRIVLIN" );
} // namespace
