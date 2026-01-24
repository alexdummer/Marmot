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
      nT( { materialProperties[12], materialProperties[13], materialProperties[14] } )
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

    initializeStateLayout();
  }

  double Unistrand::getDensity()
  {
    throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": Density not implemented yet!" );
  }

  void Unistrand::computeStress( state3D&        state,
                                 double*         dStress_dStrain,
                                 const double*   dStrain,
                                 const timeInfo& timeInfo ) const
  {
    mVector6d             S( state.stress.data() );
    Map< const Vector6d > dE( dStrain );
    mMatrix6d             mC( dStress_dStrain );
    Map< Vector9d >       alpha( stateLayout.getPtr( state.stateVars, "kappa" ) );

    // Zero strain  increment check
    if ( ( dE.array() == 0 ).all() ) {
      mC = globalElasticStiffnessTensor;
      return;
    }

    const double& viscosity = materialProperties[33];

    // transform old stress and strain increment to local coordinate system
    const Vector6d dStrainLocal   = transformationMatrixStrain * dE;
    const Vector6d stressLocalOld = transformationMatrixStress * S;

    // compute trial stress
    Vector6d stressLocalTrial = stressLocalOld + localElasticStiffnessTensor * dStrainLocal;

    // initialize trial state for plasticity with trial stress and old alpha
    UnistrandPlasticity::MaterialState trialState = { stressLocalTrial, alpha };

    // plasticity
    UnistrandPlasticity plasticity = UnistrandPlasticity(
      // 9 strengths
      { materialProperties[15],
        materialProperties[16],
        materialProperties[17],
        materialProperties[18],
        materialProperties[19],
        materialProperties[20],
        materialProperties[21],
        materialProperties[22],
        materialProperties[23] },
      // 9 hardening moduli
      { materialProperties[24],
        materialProperties[25],
        materialProperties[26],
        materialProperties[27],
        materialProperties[28],
        materialProperties[29],
        materialProperties[30],
        materialProperties[31],
        materialProperties[32] },
      this->characteristicElementLength );

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

        auto& dT = timeInfo.dT;
        // update stress and alpha with Duvaut-Lions viscoplasticity
        S = transformationMatrixStressInv * ( trialState.stress + dT / viscosity * result.materialState.stress ) /
            ( 1.0 + dT / viscosity );
        alpha = ( trialState.alpha + dT / viscosity * result.materialState.alpha ) / ( 1.0 + dT / viscosity );

        // update algorithmic tangent
        mC = transformationMatrixStressInv *
             ( localElasticStiffnessTensor + dT / viscosity * result.algorithmicModuli.dStress_dStrain ) /
             ( 1. + dT / viscosity ) * transformationMatrixStrain;
      }
      catch ( UnistrandPlasticity::ReturnMappingFailedException& e ) {
        throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": Return mapping failed for material number "
                                               << materialNumber );
      }
    }
    else {
      // elastic step
      S  = transformationMatrixStressInv * stressLocalTrial;
      mC = globalElasticStiffnessTensor;
    }
  }
} // namespace Marmot::Materials
