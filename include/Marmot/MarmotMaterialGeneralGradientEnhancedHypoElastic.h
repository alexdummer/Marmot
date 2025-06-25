/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once
#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedMechanical.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Eigen;
using namespace Marmot;

template < int nNonlocalVariables >
class MarmotMaterialGeneralGradientEnhancedHypoElastic
  : public MarmotMaterialGeneralGradientEnhancedMechanical< nNonlocalVariables > {

public:
  using MarmotMaterialGeneralGradientEnhancedMechanical<
    nNonlocalVariables >::MarmotMaterialGeneralGradientEnhancedMechanical;

  using response  = typename MarmotMaterialGeneralGradientEnhancedMechanical< nNonlocalVariables >::response;
  using tangents  = typename MarmotMaterialGeneralGradientEnhancedMechanical< nNonlocalVariables >::tangents;
  using increment = typename MarmotMaterialGeneralGradientEnhancedMechanical< nNonlocalVariables >::increment;

  virtual void computeStress( double*       stress_,
                              double*       KLocal_,
                              double*       nonLocalRadius,
                              double*       dStress_dDeformationGradient,
                              double*       dKLocal_dDeformationGradient,
                              double*       dStress_dK,
                              double*       dKlocal_dK,
                              const double* FOld_,
                              const double* FNew_,
                              const double* KOld,
                              const double* dK,
                              const double  time,
                              const double  dT )
  {
    // Standard implemenation of the Abaqus like Hughes-Winget algorithm
    // Approximation of the algorithmic tangent in order to
    // facilitate the dCauchy_dStrain tangent provided by
    // small strain material models

    using namespace Marmot;
    const Map< const Matrix3d > FNew( FNew_ );
    const Map< const Matrix3d > FOld( FOld_ );
    Marmot::mVector6d           stress( stress_ );

    using mVectorNd      = Map< Vector< double, nNonlocalVariables > >;
    using constmVectorNd = const Map< const Vector< double, nNonlocalVariables > >;
    using mMatrixNd      = Map< Matrix< double, nNonlocalVariables, nNonlocalVariables > >;
    using mMatrix6Nd     = Map< Matrix< double, 6, nNonlocalVariables > >;

    Matrix6d CJaumann                = Matrix6d::Zero();
    Vector6d dK_LocalDStretchingRate = Vector6d::Zero();

    Marmot::NumericalAlgorithms::HughesWinget
      hughesWingetIntegrator( FOld, FNew, Marmot::NumericalAlgorithms::HughesWinget::Formulation::AbaqusLike );

    auto dEps = hughesWingetIntegrator.getStrainIncrement();
    stress    = hughesWingetIntegrator.rotateTensor( stress );

    response        res = { stress, mVectorNd( KLocal_ ), mVectorNd( nonLocalRadius ) };
    tangents        tan = { mMatrix6d( dStress_dDeformationGradient ),
                            mMatrix6Nd( dKLocal_dDeformationGradient ),
                            mMatrix6Nd( dStress_dK ),
                            mMatrixNd( dKlocal_dK ) };
    const increment inc = { dEps, constmVectorNd( KOld ), constmVectorNd( dK ), dT, time };

    computeStress( res, tan, inc );

    TensorMap< Eigen::Tensor< double, 3 > > dS_dF( tan.dStressddStrain.data(), 6, 3, 3 );
    Map< Matrix3d >                         dKLocal_dF( tan.dKLocalddStrain.data() );

    Matrix3d FInv = FNew.inverse();
    dS_dF         = hughesWingetIntegrator.compute_dS_dF( stress, FInv, CJaumann );
    dKLocal_dF    = hughesWingetIntegrator.compute_dScalar_dF( FInv, dK_LocalDStretchingRate );
  }

  virtual void computeStress( response& res, tangents& tan, const increment& inc ) = 0;

  using MarmotMaterialGeneralGradientEnhancedMechanical< nNonlocalVariables >::computePlaneStress;
  virtual void computePlaneStress( response& res, tangents& tan, const increment& inc )
  {
    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    Map< VectorXd > stateVars( this->stateVars, this->nStateVars );

    VectorXd  stateVarsOld = stateVars;
    response  resTemp      = { res.stress, res.KLocal, res.c };
    increment incTemp      = { inc.dStrain, inc.K, inc.dK, inc.time, inc.dT };
    // assumption of isochoric deformation for initial guess

    double residual          = 1;
    double tangentCompliance = 1.;
    // assumption of isochoric deformation for initial guess
    double strainIncrement = ( -incTemp.dStrain( 0 ) - incTemp.dStrain( 1 ) );
    incTemp.dStrain( 2 )   = strainIncrement;

    int planeStressCount = 1;
    while ( true ) {

      // set old response
      resTemp = { res.stress, res.KLocal, res.c };
      // set old state variables
      stateVars = stateVarsOld;
      // compute stress
      computeStress( resTemp, tan, incTemp );

      // evauate residual
      residual = std::abs( resTemp.stress.array().abs()[2] / std::max( resTemp.stress.array().abs().maxCoeff(), 1. ) );
      /* std::cout << "residual: " << residual << std::endl; */
      if ( ( residual < 1e-10 && std::abs( strainIncrement ) < 1e-8 ) || ( planeStressCount > 7 && residual < 1e-5 ) ) {
        break;
      }

      // correct strain increment
      tangentCompliance = 1. / tan.dStressddStrain( 2, 2 );
      if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 ) {
        tangentCompliance = 1e10;
      }

      strainIncrement = -resTemp.stress( 2 ) * tangentCompliance;
      incTemp.dStrain( 2 ) += strainIncrement;

      planeStressCount += 1;
      if ( planeStressCount > 13 ) {
        throw std::runtime_error( "PlaneStressWrapper requires cutback" );
      }
    }

    res = { resTemp.stress, resTemp.KLocal, resTemp.c };
  }
};
