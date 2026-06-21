#include "Marmot/FiniteStrainGradientVonMises.h"
#include "Marmot/MarmotMaterialGradientPlasticityFiniteStrainFactory.h"

namespace Marmot::Materials {
  const bool isFiniteStrainGradientVonMisesRegistered = MarmotLibrary::
    MarmotMaterialGradientPlasticityFiniteStrainFactory< 1 >::registerMaterial< FiniteStrainGradientVonMises >(
      "FiniteStrainGradientVonMises" );
}