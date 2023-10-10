#include "Marmot/WoodCreepDamagePlasticity.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;
  using namespace ContinuumMechanics::Elasticity;

  WoodCreepDamagePlasticity::WoodCreepDamagePlasticity( const double* materialProperties,
                                                        int           nMaterialProperties,
                                                        int           materialNumber )
    : MarmotMaterialHypoElastic::MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      anisotropicType( static_cast< Type >( nMaterialProperties ) ),
      E1( materialProperties[0] ),
      E2( materialProperties[1] ),
      E3( materialProperties[2] ),
      nu12( materialProperties[3] ),
      nu23( materialProperties[4] ),
      nu13( materialProperties[5] ),
      G12( materialProperties[6] ),
      G23( materialProperties[7] ),
      G13( materialProperties[8] ),
      n1( { materialProperties[9], materialProperties[10], materialProperties[11] } ),
      n2( { materialProperties[12], materialProperties[13], materialProperties[14] } )
  {
  }

  void WoodCreepDamagePlasticity::computeStress( double*       stress,
                                                 double*       dStressDDStrain,
                                                 const double* dStrain,
                                                 const double* timeOld,
                                                 const double  dT,
                                                 double&       pNewDT )
  {
    Matrix6d localStiffnessTensor = Orthotropic::stiffnessTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G13 );

    Matrix3d localCoordinateSystem = Marmot::Math::orthonormalCoordinateSystem( n1, n2 );

    // strain and stress transformation matrices
    using namespace ContinuumMechanics::VoigtNotation::Transformations;
    Matrix6d transformationStrainInv = transformationMatrixStrainVoigt( localCoordinateSystem ).inverse();
    Matrix6d transformationStress    = transformationMatrixStressVoigt( localCoordinateSystem );

    // transformation Cel into global coordinate system
    C = transformationStrainInv * localStiffnessTensor * transformationStress;

    mVector6d             S( stress );
    Map< const Vector6d > dE( dStrain );
    mMatrix6d             mC( dStressDDStrain );
    mC = C;

    // Zero strain  increment check
    if ( ( dE.array() == 0 ).all() )
      return;

    S += mC * dE;
  }
} // namespace Marmot::Materials
