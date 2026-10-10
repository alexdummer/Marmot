#include "Marmot/FiniteStrainGradientVonMises.h"
#include "Marmot/MarmotMaterialGradientPlasticityFiniteStrainFactory.h"

namespace Marmot::Materials {
  const bool isFiniteStrainGradientVonMisesRegistered = MarmotLibrary::
    MarmotMaterialGradientPlasticityFiniteStrainFactory< 1 >::registerMaterial< FiniteStrainGradientVonMises >(
      "FINITESTRAINGRADIENTVONMISES" );
  // the former spelling, kept as an alias
  const bool isFiniteStrainGradientVonMisesAliasRegistered = MarmotLibrary::
    MarmotMaterialGradientPlasticityFiniteStrainFactory< 1 >::registerMaterial< FiniteStrainGradientVonMises >(
      "FiniteStrainGradientVonMises" );
} // namespace Marmot::Materials