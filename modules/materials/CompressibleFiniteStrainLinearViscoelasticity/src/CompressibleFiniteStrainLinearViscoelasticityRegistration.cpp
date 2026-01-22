#include "Marmot/CompressibleFiniteStrainLinearViscoelasticity.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool CompressibleFiniteStrainLinearViscoelasticityRegistered = MarmotMaterialFiniteStrainFactory::
      registerMaterial< CompressibleFiniteStrainLinearViscoelasticity >( "LINEARVISCOELASTICCOMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
