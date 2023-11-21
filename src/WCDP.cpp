#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/WCP.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  WoodCreepPlasticity::WoodCreepPlasticity( const double* materialProperties,
                                            int           nMaterialProperties,
                                            int           materialNumber )
    : MarmotMaterialHypoElastic::MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      ER( materialProperties[0] ),
      ET( materialProperties[1] ),
      EL( materialProperties[2] ),
      nuTR( materialProperties[3] ),
      nuLR( materialProperties[4] ),
      nuLT( materialProperties[5] ),
      GRT( materialProperties[6] ),
      GRL( materialProperties[7] ),
      GTL( materialProperties[8] ),
      nR( { materialProperties[9], materialProperties[10], materialProperties[11] } ),
      nT( { materialProperties[12], materialProperties[13], materialProperties[14] } )
  {
    localElasticStiffnessTensor = Orthotropic::stiffnessTensor( ER, ET, EL, nuTR, nuLR, nuLT, GRT, GTL, GRL );
    localCoordinateSystem       = Marmot::Math::orthonormalCoordinateSystem( nR, nT );

    using namespace ContinuumMechanics::VoigtNotation::Transformations;
    transformationMatrixStrain = transformationMatrixStrainVoigt( localCoordinateSystem );
    transformationMatrixStress = transformationMatrixStressVoigt( localCoordinateSystem );

    std::cout << "Cel = " << localElasticStiffnessTensor << std::endl;
    std::cout << "Cel ^ -1 = " << localElasticStiffnessTensor.inverse() << std::endl;
  }

  void WoodCreepPlasticity::computeStress( double*       stress,
                                           double*       dStressDDStrain,
                                           const double* dStrain,
                                           const double* timeOld,
                                           const double  dT,
                                           double&       pNewDT )
  {
    mVector6d             S( stress );
    Map< const Vector6d > dE( dStrain );
    mMatrix6d             mC( dStressDDStrain );

    mC = transformationMatrixStress.inverse() * localElasticStiffnessTensor * transformationMatrixStrain;

    // Zero strain  increment check
    if ( ( dE.array() == 0 ).all() )
      return;

    const Vector6d dStrainLocal   = transformationMatrixStrain * dE;
    const Vector6d stressLocalOld = transformationMatrixStress * S;

    Vector6d stressLocalTrial = stressLocalOld + localElasticStiffnessTensor * dStrainLocal;

    // transform stress into global coordinate system
    S = transformationMatrixStress.inverse() * stressLocalTrial;
  }
} // namespace Marmot::Materials
