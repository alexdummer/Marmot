#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace Marmot::Registration; // factory/registry namespace, distinct from this file's own ::Registration
                                          // scope

    const static bool
      CompressibleNeoHookeRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial< CompressibleNeoHooke >(
        "COMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
