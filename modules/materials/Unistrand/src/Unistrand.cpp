#include "Marmot/Unistrand.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/UnistrandPlasticity.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  Unistrand::Unistrand( const double* materialProperties, int nMaterialProperties, int materialNumber )
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
      nT( { materialProperties[12], materialProperties[13], materialProperties[14] } ),
      plasticity( { { materialProperties[15],      // r11t
                      materialProperties[16],      // r11c
                      materialProperties[17],      // r22t
                      materialProperties[18],      // r22c
                      materialProperties[19],      // r33t
                      materialProperties[20],      // r33c
                      materialProperties[21],      // r12
                      materialProperties[22],      // r13
                      materialProperties[23] } } ) // r23
  {

    localElasticStiffnessTensor  = Orthotropic::stiffnessTensor( ER, ET, EL, nuTR, nuLR, nuLT, GRT, GTL, GRL );
    localElasticComplianceTensor = Orthotropic::complianceTensor( ER, ET, EL, nuTR, nuLR, nuLT, GRT, GTL, GRL );
    localCoordinateSystem        = Marmot::Math::orthonormalCoordinateSystem( nR, nT );

    using namespace ContinuumMechanics::VoigtNotation::Transformations;
    transformationMatrixStrain    = transformationMatrixStrainVoigt( localCoordinateSystem );
    transformationMatrixStrainInv = transformationMatrixStrain.colPivHouseholderQr().inverse();
    transformationMatrixStress    = transformationMatrixStressVoigt( localCoordinateSystem );
    transformationMatrixStressInv = transformationMatrixStress.colPivHouseholderQr().inverse();

    globalElasticStiffnessTensor = transformationMatrixStrainInv * localElasticStiffnessTensor *
                                   transformationMatrixStrain;
  }

  void Unistrand::computeStress( double*       stress,
                                 double*       dStressDDStrain,
                                 const double* dStrain,
                                 const double* timeOld,
                                 const double  dT,
                                 double&       pNewDT )
  {
    mVector6d             S( stress );
    Map< const Vector6d > dE( dStrain );
    mMatrix6d             mC( dStressDDStrain );
    auto&                 alpha = managedStateVars->alpha;

    // Zero strain  increment check
    if ( ( dE.array() == 0 ).all() ) {
      mC = globalElasticStiffnessTensor;
      return;
    }

    // transform old stress and strain increment to local coordinate system
    const Vector6d dStrainLocal   = transformationMatrixStrain * dE;
    const Vector6d stressLocalOld = transformationMatrixStress * S;

    // compute trial stress
    Vector6d stressLocalTrial = stressLocalOld + localElasticStiffnessTensor * dStrainLocal;

    // initialize trial state for plasticity with trial stress and old alpha
    UnistrandPlasticity::MaterialState trialState = { stressLocalTrial, alpha };

    if ( plasticity.checkIfYielding( trialState ) ) {
      // plastic step

      // initialize old state with old stress and old alpha for the possible use in the smart return mapping
      UnistrandPlasticity::MaterialState oldState = { stressLocalOld, alpha };
      try {

        // perform (smart) return mapping
        UnistrandPlasticity::ReturnMapResult result = plasticity
                                                        .performSmartReturnMapping( trialState,
                                                                                    oldState,
                                                                                    localElasticStiffnessTensor,
                                                                                    localElasticComplianceTensor );
        // update stress and alpha
        S     = transformationMatrixStressInv * result.materialState.stress;
        alpha = result.materialState.alpha;

        // update algorithmic tangent
        mC = transformationMatrixStressInv * result.algorithmicModuli.dStress_dStrain * transformationMatrixStrain;
      }
      catch ( UnistrandPlasticity::ReturnMappingFailedException& e ) {
        pNewDT = 0.5;
        return;
      }
    }
    else {
      // elastic step
      S  = transformationMatrixStressInv * stressLocalTrial;
      mC = globalElasticStiffnessTensor;
    }
  }

  StateView Unistrand::getStateView( const std::string& stateName )
  {
    return managedStateVars->getStateView( stateName );
  }

  void Unistrand::assignStateVars( double* stateVars, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->managedStateVars = std::make_unique< UnistrandStateVarManager >( stateVars );

    return MarmotMaterialHypoElastic::assignStateVars( stateVars, nStateVars );
  }
} // namespace Marmot::Materials
