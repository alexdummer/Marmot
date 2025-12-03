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
#include <map>

namespace Marmot::Materials {

  class UnistrandPlasticity {

    using Diag6d   = Eigen::DiagonalMatrix< double, 6 >;
    using Vector9d = Eigen::Vector< double, 9 >;

  public:
    struct ReturnMappingFailedException : std::exception {};

    struct strength {
      double r11t; // x-1 tension
      double r11c; // x-1 compression
      double r22t; // x-2 tension
      double r22c; // x-2 compression
      double r33t; // x-3 tension
      double r33c; // x-3 compression
      double s12;  // shear x1-x2
      double s13;  // shear x1-x3
      double s23;  // shear x2-x3
    };

    enum class failureMode {
      tension1,
      compression1,
      tension2,
      compression2,
      tension3,
      compression3,
      shear12,
      shear13,
      shear23
    };

    struct MaterialState {
      Marmot::Vector6d stress;
      Vector9d         alpha;
    };

    struct AlgorithmicModuli {
      Marmot::Matrix6d dStress_dStrain;
      Marmot::Matrix6d dStrainPlastic_dStrain;
    };

    struct ReturnMapResult {
      MaterialState     materialState;
      Marmot::Vector6d  dStrainPlastic;
      AlgorithmicModuli algorithmicModuli;
    };

    UnistrandPlasticity( const strength& );

    ReturnMapResult performSmartReturnMapping( const MaterialState&    trialState,
                                               const MaterialState&    stateOld,
                                               const Marmot::Matrix6d& elasticStiffness,
                                               const Marmot::Matrix6d& elasticCompliance );

    AlgorithmicModuli computeAlgorithmicModuli( const Eigen::MatrixXd&  dF_dX,
                                                const Marmot::Matrix6d& elasticStiffness,
                                                const Marmot::Matrix6d& elasticCompliance );

    bool checkIfYielding( const MaterialState& trialState )
    {
      const auto res = evaluateYieldFunctions( trialState );
      if ( res.any() )
        return true;
      else
        return false;
    }

  private:
    constexpr static int nStress            = 6;
    constexpr static int nInternalVariables = 9;
    constexpr static int nYieldSurfaces     = 9;

    constexpr static int idxS        = 0;
    constexpr static int idxInternal = idxS + nStress;
    constexpr static int idxSurface  = idxInternal + 1;

    constexpr static int    nMaxInnerNewtonCycles    = 5;
    constexpr static int    nMaxInnerNewtonCyclesAlt = 10;
    constexpr static double innerNewtonTol           = 1e-14;
    constexpr static double innerNewtonRTol          = 1e-10;
    constexpr static double innerNewtonTolAlt        = 1e-8;
    constexpr static double innerNewtonRTolAlt       = 1e-6;

    constexpr static double yieldSurfaceTol = 1e-8;

    const strength& st;

    Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces > yieldSurfCombiManager;

    typedef Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces >::YieldSurfFlagArr
      YieldSurfFlagArr;
    typedef Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces >::YieldSurfResArr
      YieldSurfResArr;

    YieldSurfFlagArr evaluateYieldFunctions( const MaterialState& trialState )
    {
      YieldSurfFlagArr activeSurfaces;
      activeSurfaces.setConstant( false );
      // check for active yield surfaces
      for ( int i = 0; i < nYieldSurfaces; i++ ) {
        if ( yieldFunction( trialState.stress, failureMode( i ), trialState.alpha( i ) ) > yieldSurfaceTol )
          activeSurfaces( i ) = true;
      }
      return activeSurfaces;
    }

    template < typename T >
    Marmot::VectorXt< T > R( const Marmot::VectorXt< T >& X,
                             const MaterialState&         trialState,
                             const Marmot::Matrix6d&      elasticStiffness,
                             const YieldSurfFlagArr&      activeSurfaces )
    {

      // extract sigma_n+1
      const Marmot::Vector6t< T >& stress = X.head( nStress );
      const Matrix6d&              Cel    = elasticStiffness;

      // initialize with sigma_n+1 - sigma_n+1^trial
      Marmot::VectorXt< T > R = VectorXt< T >::Zero( X.size() );
      R.head( nStress )       = stress - trialState.stress;

      int idxCurrentSurface  = idxSurface;
      int idxCurrentInternal = idxInternal;

      for ( int i = 0; i < nYieldSurfaces; i++ ) {
        if ( activeSurfaces( i ) ) {
          // yield surface is active
          // get corresponding internal variable and consistency parameter
          const T alpha   = X( idxCurrentInternal );
          const T dLambda = X( idxCurrentSurface );
          // get cooresponding failure mode
          const auto          fm                   = failureMode( i );
          const Vector6t< T > plasticFlowDirection = a[fm];
          R.head( nStress ) += Cel * dLambda * plasticFlowDirection;
          R( idxCurrentInternal ) = alpha - dLambda - trialState.alpha( i );
          R( idxCurrentSurface )  = yieldFunction( stress, fm, alpha );
          idxCurrentInternal += 2;
          idxCurrentSurface += 2;
        }
      }

      return R;
    }

    // map for strength tensor a
    std::map< failureMode, Vector6d > a{ { failureMode::tension1, { 1 / st.r11t, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                                         { failureMode::compression1, { -1 / st.r11c, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                                         { failureMode::tension2, { 0.0, 1 / st.r22t, 0.0, 0.0, 0.0, 0.0 } },
                                         { failureMode::compression2, { 0.0, -1 / st.r22c, 0.0, 0.0, 0.0, 0.0 } },
                                         { failureMode::tension3, { 0.0, 0.0, 1 / st.r33t, 0.0, 0.0, 0.0 } },
                                         { failureMode::compression3, { 0.0, 0.0, -1 / st.r33c, 0.0, 0.0, 0.0 } },
                                         { failureMode::shear12, { 0.0, 0.0, 0.0, 1 / st.s12, 0.0, 0.0 } },
                                         { failureMode::shear13, { 0.0, 0.0, 0.0, 0.0, 1 / st.s13, 0.0 } },
                                         { failureMode::shear23, { 0.0, 0.0, 0.0, 0.0, 0.0, 1 / st.s23 } } };

    template < typename T >
    T yieldFunction( const Marmot::Vector6t< T >& stress, const failureMode fm, const T alpha )
    {
      // remove sign from shear components

      Marmot::Vector6t< T > stressMod = stress;
      stressMod( 3 )                  = abs( stressMod( 3 ) );
      stressMod( 4 )                  = abs( stressMod( 4 ) );
      stressMod( 5 )                  = abs( stressMod( 5 ) );

      const T res = a[fm].dot( stressMod ) - 1.0;
      return res;
    }

    Marmot::Vector6d computePlasticStrainIncrement( const Marmot::Vector6d& stressNew,
                                                    const Marmot::Vector6d& trialStress,
                                                    const Marmot::Matrix6d& Cel );

    Eigen::VectorXd createResidualScaleVector( const Eigen::VectorXd& leftHandSideTarget );
    double          ResidualNorm( const Eigen::VectorXd& ResidualX, const Eigen::MatrixXd& scaleMatrix );

    ReturnMapResult performReturnMapping( const MaterialState&    trialState,
                                          const MaterialState&    initialGuessState,
                                          const Marmot::Matrix6d& elasticStiffness,
                                          const Marmot::Matrix6d& elasticCompliance );

    ReturnMapResult returnMappingAttempt( const MaterialState&    trialState,
                                          const MaterialState&    initialGuessState,
                                          const Marmot::Matrix6d& elasticStiffness,
                                          const Marmot::Matrix6d& elasticCompliance,
                                          YieldSurfFlagArr&       activeSurfaces );

    MaterialState extractMaterialState( const Eigen::VectorXd&  X,
                                        const MaterialState&    trialState,
                                        const YieldSurfFlagArr& activeSurfaces );

    std::tuple< Eigen::VectorXd, Eigen::MatrixXd > dR_dX( const Eigen::VectorXd&  X,
                                                          const MaterialState&    trialState,
                                                          const Marmot::Matrix6d& elasticStiffness,
                                                          YieldSurfFlagArr&       activeSurfaces );
  };
} // namespace Marmot::Materials
