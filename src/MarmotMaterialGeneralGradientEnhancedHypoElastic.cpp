#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"

using namespace Eigen;

void MarmotMaterialGeneralGradientEnhancedHypoElastic::computeStress( double*       stress_,
                                                               double*       K_local,
                                                               double&       nonLocalRadius,
                                                               double*       dStressDDDeformationGradient_,
                                                               double*       dK_localDDeformationGradient_,
                                                               double*       dStressDK,
                                                               double*       dKlocal_dK,
                                                               const double* FOld_,
                                                               const double* FNew_,
                                                               const double*  KOld,
                                                               const double*  dK,
                                                               const double time,
                                                               const double  dT)
{
  // Standard implemenation of the Abaqus like Hughes-Winget algorithm
  // Approximation of the algorithmic tangent in order to
  // facilitate the dCauchy_dStrain tangent provided by
  // small strain material models

  using namespace Marmot;
  const Map< const Matrix3d > FNew( FNew_ );
  const Map< const Matrix3d > FOld( FOld_ );
  Marmot::mVector6d           stress( stress_ );

  Matrix6d CJaumann                = Matrix6d::Zero();
  Vector6d dK_LocalDStretchingRate = Vector6d::Zero();

  Marmot::NumericalAlgorithms::HughesWinget
    hughesWingetIntegrator( FOld, FNew, Marmot::NumericalAlgorithms::HughesWinget::Formulation::AbaqusLike );

  auto dEps = hughesWingetIntegrator.getStrainIncrement();
  stress    = hughesWingetIntegrator.rotateTensor( stress );

  computeStress( stress.data(),
                 K_local,
                 nonLocalRadius,
                 CJaumann.data(),
                 dK_LocalDStretchingRate.data(),
                 dStressDK,
                 dKlocal_dK,
                 dEps.data(),
                 KOld,
                 dK,
                 time,
                 dT);

  TensorMap< Eigen::Tensor< double, 3 > > dS_dF( dStressDDDeformationGradient_, 6, 3, 3 );
  Map< Matrix3d >                         dKLocal_dF( dK_localDDeformationGradient_ );

  Matrix3d FInv = FNew.inverse();
  dS_dF         = hughesWingetIntegrator.compute_dS_dF( stress, FInv, CJaumann );
  dKLocal_dF    = hughesWingetIntegrator.compute_dScalar_dF( FInv, dK_LocalDStretchingRate );
}
