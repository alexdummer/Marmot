#include "Marmot/WCP.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/WCPPlasticity.h"

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
      nT( { materialProperties[12], materialProperties[13], materialProperties[14] } ),
      plasticityParams( {
        { materialProperties[15], materialProperties[16] }, // strength parameters radial direction (compression)
        { materialProperties[17],
          materialProperties[18],
          materialProperties[19] },                         // hardening parameters radial direction
        { materialProperties[20], materialProperties[21] }, // strength parameters tangential direction (compression)
        { materialProperties[22],
          materialProperties[23],
          materialProperties[24] },                         // hardening parameters tangential direction
        { materialProperties[25], materialProperties[26] }, // strength parameters longitudinal direction (compression)
        { materialProperties[27],
          materialProperties[28],
          materialProperties[29] },                         // hardening parameters longitudinal direction
        { materialProperties[30], materialProperties[31] }, // strength parameters longitudinal direction (tension)
        { materialProperties[32], materialProperties[33] }, // strength parameters shear in RT plane
        { materialProperties[34], materialProperties[35] }, // strength parameters shear in RL plane
        { materialProperties[36], materialProperties[37] }, // strength parameters shear in TL plane
        materialProperties[38],                             // actual moisture content
      } )
  {
    localElasticStiffnessTensor  = Orthotropic::stiffnessTensor( ER, ET, EL, nuTR, nuLR, nuLT, GRT, GTL, GRL );
    localElasticComplianceTensor = Orthotropic::complianceTensor( ER, ET, EL, nuTR, nuLR, nuLT, GRT, GTL, GRL );
    localCoordinateSystem        = Marmot::Math::orthonormalCoordinateSystem( nR, nT );

    using namespace ContinuumMechanics::VoigtNotation::Transformations;
    transformationMatrixStrain    = transformationMatrixStrainVoigt( localCoordinateSystem );
    transformationMatrixStrainInv = transformationMatrixStrain.colPivHouseholderQr().inverse();
    transformationMatrixStress    = transformationMatrixStressVoigt( localCoordinateSystem );
    transformationMatrixStressInv = transformationMatrixStress.colPivHouseholderQr().inverse();
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
    double&               alphaR = managedStateVars->alphaR;
    double&               alphaT = managedStateVars->alphaT;
    double&               alphaL = managedStateVars->alphaL;

    mC = transformationMatrixStress.inverse() * localElasticStiffnessTensor * transformationMatrixStrain;

    // Zero strain  increment check
    if ( ( dE.array() == 0 ).all() )
      return;

    const Vector6d dStrainLocal   = transformationMatrixStrain * dE;
    const Vector6d stressLocalOld = transformationMatrixStress * S;

    Vector6d stressLocalTrial = stressLocalOld + localElasticStiffnessTensor * dStrainLocal;

    WCPPlasticity::MaterialState trialState = { stressLocalTrial, { alphaR, alphaT, alphaL } };

    WCPPlasticity plasticity = WCPPlasticity( plasticityParams );

    if ( plasticity.checkIfYielding( trialState ) ) {
      WCPPlasticity::MaterialState oldState = { stressLocalOld, { alphaR, alphaT, alphaL } };
      try {
        WCPPlasticity::ReturnMapResult result = plasticity.performSmartReturnMapping( trialState,
                                                                                      oldState,
                                                                                      localElasticStiffnessTensor,
                                                                                      localElasticComplianceTensor );
        S                                     = transformationMatrixStressInv * result.materialState.stress;
        alphaR                                = result.materialState.alpha( 0 );
        alphaT                                = result.materialState.alpha( 1 );
        alphaL                                = result.materialState.alpha( 2 );
        mC = transformationMatrixStressInv * result.algorithmicModuli.dStress_dStrain * transformationMatrixStrain;
      }
      catch ( WCPPlasticity::ReturnMappingFailedException& e ) {
        pNewDT = 0.5;
        return;
      }
    }
    else {
      S  = transformationMatrixStressInv * stressLocalTrial;
      mC = transformationMatrixStressInv * localElasticStiffnessTensor * transformationMatrixStrain;
    }
  }

  StateView WoodCreepPlasticity::getStateView( const std::string& stateName )
  {
    return managedStateVars->getStateView( stateName );
  }

  void WoodCreepPlasticity::assignStateVars( double* stateVars, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->managedStateVars = std::make_unique< WCPStateVarManager >( stateVars );

    return MarmotMaterialHypoElastic::assignStateVars( stateVars, nStateVars );
  }
} // namespace Marmot::Materials
