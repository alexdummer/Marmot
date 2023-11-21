/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 *
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/YieldSurfaceCombinationManager.h"

namespace Marmot::Materials {

  class WCPPlasticity {
  public:
    struct ReturnMappingFailedException : std::exception {};

    struct MaterialParameters {
      double fcL;
      double fcT;
      double fcR;
      double ftL;
      double fsRT;
      double fsRL;
      double fsTL;
    };

    struct MaterialState {
      Marmot::Vector6d stress;
      double           alphaL;
      double           alphaR;
      double           alphaT;
    };

    struct AlgorithmicModuli {
      Marmot::Matrix6d dSdE;
      Marmot::Matrix6d dEpdE;
    };

    struct ReturnMapResult {
      MaterialState    materialState;
      Marmot::Vector6d dStrainPlastic;
      Eigen::MatrixXd  Jacobian;
    };

    WCPPlasticity( const MaterialParameters& );

    bool checkIfYielding( const MaterialState& trialState );

    ReturnMapResult performReturnMapping( const MaterialState&    trialState,
                                          const Marmot::Matrix6d& elasticStiffness,
                                          const Marmot::Matrix6d& elasticCompliance );

    AlgorithmicModuli computeAlgorithmicModuli( const ReturnMapResult&,
                                                const Marmot::Matrix6d&,
                                                const Marmot::Matrix6d& );

  private:
    constexpr static int nStress            = 6;
    constexpr static int nInternalVariables = 3;
    constexpr static int nYieldSurfaces     = 3;

    constexpr static int idxS        = 0;
    constexpr static int idxInternal = idxS + nStress;
    constexpr static int idxSurface  = idxInternal + nInternalVariables;

    const MaterialParameters materialParameters;

    const Vector6d aL, aT, aR;
    const Matrix6d bL, bT, bR;

    Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces > yieldSurfCombiManager;

    typedef Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces >::YieldSurfFlagArr
      YieldSurfFlagArr;
    typedef Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces >::YieldSurfResArr
      YieldSurfResArr;

    YieldSurfResArr evaluateYieldFunctions( const MaterialState& materialState );

    template < typename T >
    Marmot::VectorXt< T > F( const Marmot::VectorXt< T >& X,
                             const Marmot::Matrix6d&      elasticStiffness,
                             const YieldSurfFlagArr&      activeSurfaces )
    {

      const Marmot::VectorXt< T >& stress = X.segment< nStress >( idxS );
      const T                      alphaL = X( idxInternal );
      const T                      alphaR = X( idxInternal + 1 );
      const T                      alphaT = X( idxInternal + 2 );

      const Matrix6d& Cel = elasticStiffness;

      Marmot::VectorXt< T > F( X );

      // initialize with sigma_n+1 , m_n+1, alphaP_n_+1, ...
      F( idxInternal ) = alphaL;

      int idxCurrentSurface = idxSurface;
      if ( activeSurfaces( 0 ) ) {
        const T       dLambda = X( idxCurrentSurface );
        Vector6t< T > dfL_dStress;
        F.head( nStress ) += Cel * dLambda * dfL_dStress;
        F( idxInternal ) -= dLambda;
        F( idxCurrentSurface ) = fL( stress, alphaL );
        idxCurrentSurface++;
      }
      if ( activeSurfaces( 1 ) ) {
        const T dLambda = X( idxCurrentSurface );
        F.head( nStress ) += Cel * dLambda * dgRankine_dStress( stress );
        F( idxInternal ) -= dLambda;
        F( idxCurrentSurface ) = fR( stress, alphaR );
        idxCurrentSurface++;
      }
      if ( activeSurfaces( 2 ) ) {
        const T dLambda = X( idxCurrentSurface );
        F.head( nStress ) += Cel * dLambda * dgCutOff_dStress( stress );
        F( idxInternal ) -= dLambda;
        F( idxCurrentSurface ) = fT( stress, alphaT );
        idxCurrentSurface++;
      }

      return F;
    }

    Eigen::MatrixXd dFdX( const Eigen::VectorXd&  X,
                          const Marmot::Matrix6d& elasticStiffness,
                          const YieldSurfFlagArr& activeSurfaces );

    template < typename T >
    T fL( const Marmot::Vector6t< T >& stress, T alphaL )
    {
      T res = aL.dot( stress ) + stress.dot( bL.dot( stress ) ) + qL( alphaL ) - 1.0;
      return res;
    }

    template < typename T >
    T fT( const Marmot::Vector6t< T >& stress, T alphaT )
    {
      T res = aT.dot( stress ) + stress.dot( bT.dot( stress ) ) + qT( alphaT ) - 1.0;
      return res;
    }

    template < typename T >
    T fR( const Marmot::Vector6t< T >& stress, T alphaR )
    {
      T res = aR.dot( stress ) + stress.dot( bR.dot( stress ) ) + qR( alphaR ) - 1.0;
      return res;
    }

    template < typename T >
    T qR( const T alphaR )
    {
    }

    Marmot::Vector6d computePlasticStrainIncrement( const Marmot::Vector6d& stressNew,
                                                    const Marmot::Vector6d& trialStress,
                                                    const Marmot::Matrix6d& Cel );

    Eigen::VectorXd createResidualScaleVector( const Eigen::VectorXd& leftHandSideTarget );
    double          ResidualNorm( const Eigen::VectorXd& ResidualX, const Eigen::MatrixXd& scaleMatrix );

    ReturnMapResult returnMappingAttempt( const MaterialState&    trialState,
                                          const Marmot::Matrix6d& elasticStiffness,
                                          const Marmot::Matrix6d& elasticCompliance,
                                          YieldSurfFlagArr&       activeSurfaces );

    ReturnMapResult retryReturnMapping( const MaterialState&    trialState,
                                        const Marmot::Matrix6d& elasticStiffness,
                                        const Marmot::Matrix6d& elasticCompliance );

    MaterialState extractMaterialState( const Eigen::VectorXd& X );
  };

} // namespace Marmot::Materials
  //
