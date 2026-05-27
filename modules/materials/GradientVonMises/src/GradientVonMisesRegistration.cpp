#include "Marmot/GradientVonMises.h"
#include "Marmot/MarmotMaterialGradientPlasticityHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool GradientVonMisesIsRegistered = MarmotMaterialGradientPlasticityHypoElasticFactory::
      registerMaterial< GradientVonMises >( "GRADIENTVONMISES" );

  } // namespace Registration

} // namespace Marmot::Materials
