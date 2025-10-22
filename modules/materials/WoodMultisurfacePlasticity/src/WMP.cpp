#include "Marmot/WMP.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/WMPPlasticity.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  WoodMultisurfacePlasticity::WoodMultisurfacePlasticity( const double* materialProperties,
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
      plasticity( { { materialProperties[15],   // ftR
                      materialProperties[16],   // fcR
                      materialProperties[17],   // ftT
                      materialProperties[18],   // fcT
                      materialProperties[19],   // ftL
                      materialProperties[20],   // fcL
                      materialProperties[21],   // fsRT
                      materialProperties[22],   // fsRL
                      materialProperties[23] }, // fsTL
                    {
                      { materialProperties[24],
                        materialProperties[25],
                        materialProperties[26],
                        materialProperties[27],
                        materialProperties[28] }, // hardening radial tension
                      { materialProperties[29],
                        materialProperties[30],
                        materialProperties[31],
                        materialProperties[32],
                        materialProperties[33] }, // hardening radial computpression
                      { materialProperties[34],
                        materialProperties[35],
                        materialProperties[36],
                        materialProperties[37],
                        materialProperties[38] }, // hardening tangential tension
                      { materialProperties[39],
                        materialProperties[40],
                        materialProperties[41],
                        materialProperties[42],
                        materialProperties[43] }, // hardening tangential compression
                      { materialProperties[44],
                        materialProperties[45],
                        materialProperties[46],
                        materialProperties[47],
                        materialProperties[48] }, // hardening longitudinal tension
                      { materialProperties[49],
                        materialProperties[50],
                        materialProperties[51],
                        materialProperties[52],
                        materialProperties[53] }, // hardening longitudinal compression
                      { materialProperties[54],
                        materialProperties[55],
                        materialProperties[56],
                        materialProperties[57],
                        materialProperties[58] }, // hardening shear
                    } } )
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

  void WoodMultisurfacePlasticity::computeStress( double*       stress,
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
    WMPPlasticity::MaterialState trialState = { stressLocalTrial, alpha };

    if ( plasticity.checkIfYielding( trialState ) ) {
      // plastic step

      // initialize old state with old stress and old alpha for the possible use in the smart return mapping
      WMPPlasticity::MaterialState oldState = { stressLocalOld, alpha };
      try {

        // perform (smart) return mapping
        WMPPlasticity::ReturnMapResult result = plasticity.performSmartReturnMapping( trialState,
                                                                                      oldState,
                                                                                      localElasticStiffnessTensor,
                                                                                      localElasticComplianceTensor );
        // update stress and alpha
        S     = transformationMatrixStressInv * result.materialState.stress;
        alpha = result.materialState.alpha;

        // update algorithmic tangent
        mC = transformationMatrixStressInv * result.algorithmicModuli.dStress_dStrain * transformationMatrixStrain;
      }
      catch ( WMPPlasticity::ReturnMappingFailedException& e ) {
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

  StateView WoodMultisurfacePlasticity::getStateView( const std::string& stateName )
  {
    return managedStateVars->getStateView( stateName );
  }

  void WoodMultisurfacePlasticity::assignStateVars( double* stateVars, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->managedStateVars = std::make_unique< WMPStateVarManager >( stateVars );

    return MarmotMaterialHypoElastic::assignStateVars( stateVars, nStateVars );
  }
} // namespace Marmot::Materials
