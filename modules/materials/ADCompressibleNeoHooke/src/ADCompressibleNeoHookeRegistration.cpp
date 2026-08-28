#include "Marmot/ADCompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials::Registration {

  using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                        // scope

  const static bool
    ADCompressibleNeoHookeRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial< ADCompressibleNeoHooke >(
      "ADCOMPRESSIBLENEOHOOKE" );

} // namespace Marmot::Materials::Registration
