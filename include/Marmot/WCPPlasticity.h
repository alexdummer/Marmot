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

    using Diag6d = Eigen::DiagonalMatrix< double, 6 >;

  public:
    struct ReturnMappingFailedException : std::exception {};

    enum materialDirection { radial, tangential, longitudinal };

    struct hardeningParameters {
      double beta0;
      double beta1;
      double beta2;
    };

    struct strengthParameters {
      double strengtAtZeroMoisture;
      double moistureSensitivity;
    };

    struct MaterialParameters {
      strengthParameters  strengthCR;
      hardeningParameters betaR;
      strengthParameters  strengthCT;
      hardeningParameters betaT;
      strengthParameters  strengthCL;
      hardeningParameters betaL;
      strengthParameters  strengthTL;
      strengthParameters  strengthSRT;
      strengthParameters  strengthSRL;
      strengthParameters  strengthSTL;
      double              moisture;
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

    ReturnMapResult performReturnMapping( const MaterialState&    trialState,
                                          const Marmot::Matrix6d& elasticStiffness,
                                          const Marmot::Matrix6d& elasticCompliance );

    AlgorithmicModuli computeAlgorithmicModuli( const ReturnMapResult&,
                                                const Marmot::Matrix6d&,
                                                const Marmot::Matrix6d& );
    
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
    constexpr static int nInternalVariables = 3;
    constexpr static int nYieldSurfaces     = 3;

    constexpr static int idxS        = 0;
    constexpr static int idxInternal = idxS + nStress;
    constexpr static int idxSurface  = idxInternal + nInternalVariables;

    const MaterialParameters mp;
    
    Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces > yieldSurfCombiManager;

    typedef Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces >::YieldSurfFlagArr YieldSurfFlagArr;
    typedef Marmot::NumericalAlgorithms::YieldSurfaceCombinationManager< nYieldSurfaces >::YieldSurfResArr  YieldSurfResArr;
    
    YieldSurfFlagArr evaluateYieldFunctions( const MaterialState& trialState )
    {
      YieldSurfFlagArr activeSurfaces;
      activeSurfaces.setConstant( false );
      // check if radial direction is yielding
      double f = yieldFunction( trialState.stress, materialDirection::radial, trialState.alphaR );
      if ( f > 0.0 )
        activeSurfaces( 0 ) = true;
      // check if tangential direction is yielding
      f = yieldFunction( trialState.stress, materialDirection::tangential, trialState.alphaT );
      if ( f > 0.0 )
        activeSurfaces( 1 ) = true;
      // chekc if longitudinal direction is yielding
      f = yieldFunction( trialState.stress, materialDirection::longitudinal, trialState.alphaT );
      if ( f > 0.0 )
        activeSurfaces( 2 ) = true;
      std::cout << "activeSurfaces: " << activeSurfaces << std::endl; 
      return activeSurfaces;

    }

    hardeningParameters getHardeningParameters( const materialDirection dir )
    {
      switch ( dir ) {
      case materialDirection::radial: return mp.betaR;
      case materialDirection::tangential: return mp.betaT;
      case materialDirection::longitudinal: return mp.betaL;
      }
    };


    template < typename T >
    Marmot::VectorXt< T > F( const Marmot::VectorXt< T >& X,
                             const Marmot::Matrix6d&      elasticStiffness,
                             const YieldSurfFlagArr&      activeSurfaces )
    {

      // extract sigma_n+1
      const Marmot::VectorXt< T >& stress = X.head( nStress );
      const Matrix6d&              Cel    = elasticStiffness;

      // initialize with sigma_n+1 , m_n+1, alphaP_n_+1, ...
      Marmot::VectorXt< T > F( X );

      int idxCurrentSurface  = idxSurface;
      int idxCurrentInternal = idxInternal;
      if ( activeSurfaces( 0 ) ) {
        // yield surface  in radial direction is active
        // get corresponding internal variable and consistency
        const T alpha   = X( idxCurrentInternal );
        const T dLambda = X( idxCurrentSurface );
        // set direction
        const auto          dir                  = materialDirection::radial;
        const Vector6t< T > plasticFlowDirection = computeA( dir ) + 2. * computeB( dir ) * stress;
        F.head( nStress ) += Cel * dLambda * plasticFlowDirection;
        F( idxInternal ) -= dLambda;
        F( idxCurrentSurface ) = yieldFunction( stress, dir, alpha );
        idxCurrentInternal += 2;
        idxCurrentSurface += 2;
      }
      if ( activeSurfaces( 1 ) ) {
        // yield surface  in tenagential direction is active
        // get corresponding internal variable and consistency
        const T alpha   = X( idxCurrentInternal );
        const T dLambda = X( idxCurrentSurface );
        // set direction
        const auto dir = materialDirection::tangential;
        // compute plastic flow direction
        const Vector6t< T > plasticFlowDirection = computeA( dir ) + 2. * computeB( dir ) * stress;
        F.head( nStress ) += Cel * dLambda * plasticFlowDirection;
        F( idxInternal ) -= dLambda;
        F( idxCurrentSurface ) = yieldFunction( stress, dir, alpha );
        idxCurrentInternal += 2;
        idxCurrentSurface += 2;
      }
      if ( activeSurfaces( 2 ) ) {
        // yield surface  in longitudinal direction is active
        // get corresponding internal variable and consistency
        const T alpha   = X( idxCurrentInternal );
        const T dLambda = X( idxCurrentSurface );
        // set direction
        const auto dir = materialDirection::longitudinal;
        // compute plastic flow direction
        const Vector6t< T > plasticFlowDirection = computeA( dir ) + 2. * computeB( dir ) * stress;
        F.head( nStress ) += Cel * dLambda * plasticFlowDirection;
        F( idxInternal ) -= dLambda;
        F( idxCurrentSurface ) = yieldFunction( stress, dir, alpha );
        idxCurrentInternal += 2;
        idxCurrentSurface += 2;
      }

      return F;
    }

    Eigen::MatrixXd dFdX( const Eigen::VectorXd&  X,
                          const Marmot::Matrix6d& elasticStiffness,
                          const YieldSurfFlagArr& activeSurfaces );

    Vector6d computeA( const materialDirection dir )
    {

      Vector6d a;
      switch ( dir ) {
      case materialDirection::radial: {
        const double fcR = computeMoistureDependentStrengthParameter( mp.strengthCR );
        a                = { -1. / fcR, 0., 0., 0., 0., 0. };
        break;
      }
      case materialDirection::tangential: {
        const double fcT = computeMoistureDependentStrengthParameter( mp.strengthCT );
        a                = { 0., -1. / fcT, 0., 0., 0., 0. };
        break; }
      case materialDirection::longitudinal: {
        const double fcL = computeMoistureDependentStrengthParameter( mp.strengthCL );
        a                = { 0., 0., -1. / fcL, 0., 0., 0. };
        break; }
      };
      return a;
    }

    Diag6d computeB( const materialDirection dir )
    {

      const double fcR  = computeMoistureDependentStrengthParameter( mp.strengthCR );
      const double fcT  = computeMoistureDependentStrengthParameter( mp.strengthCT );
      const double fcL  = computeMoistureDependentStrengthParameter( mp.strengthCL );
      const double ftL  = computeMoistureDependentStrengthParameter( mp.strengthTL );
      const double fsRL = computeMoistureDependentStrengthParameter( mp.strengthSRL );
      const double fsTL = computeMoistureDependentStrengthParameter( mp.strengthSTL );
      const double fsRT = computeMoistureDependentStrengthParameter( mp.strengthSRT );

      Diag6d b;
      // compute a and b depending on the direction
      switch ( dir ) {
      case materialDirection::radial:
        // clang-format off
        b =  { 0.0000,
               0.4000 / (  fcT * fcT  ),
               0.2500 / (  fcL * ftL  ),
               0.4000 / ( fsRT * fsRT ),
               0.3300 / ( fsRL * fsRL ),
               0.3300 / ( fsTL * fsTL ) };
        // clang-format on
        break;
      case materialDirection::tangential:
        // clang-format off
        b = {  0.4000 / (  fcR * fcR  ),
               0.0000,
               0.2500 / (  fcL * ftL  ),
               0.4000 / ( fsRT * fsRT ),
               0.3300 / ( fsRL * fsRL ),
               0.3300 / ( fsTL * fsTL ) };
        // clang-format on
        break;
      case materialDirection::longitudinal:
        // clang-format off
        b = {  0.3300 / (  fcR * fcR  ),
               0.3300 / (  fcT * fcT  ),
               0.0000,
               0.2500 / ( fsRT * fsRT ),
               0.2500 / ( fsRL * fsRL ),
               0.2500 / ( fsTL * fsTL ) };
        // clang-format on
        break;
      }
      return b;
    }

    template < typename T >
    T yieldFunction( const Marmot::Vector6t< T >& stress, const materialDirection dir, const T alpha )
    {

      const Vector6d a = computeA( dir );
      const Diag6d   b = computeB( dir );

      const T res = a.dot( stress ) - stress.dot( b * stress ) + q( alpha, dir ) - 1.0;
      return res;
    }

    double computeMoistureDependentStrengthParameter( const strengthParameters& s )
    {
      const double f = s.strengtAtZeroMoisture + s.moistureSensitivity * mp.moisture;
      return f;
    }


    template < typename T >
    T q( const T alpha, const materialDirection dir )
    {

      const auto h = getHardeningParameters( dir );
      double     f;
      switch ( dir ) {
      case materialDirection::radial: f = computeMoistureDependentStrengthParameter( mp.strengthCR ); break;
      case materialDirection::tangential: f = computeMoistureDependentStrengthParameter( mp.strengthCT ); break;
      case materialDirection::longitudinal: f = computeMoistureDependentStrengthParameter( mp.strengthCL ); break;
      }

      const T res = ( h.beta0 * mp.moisture + h.beta1 ) * ( 1. - exp( -h.beta2 * alpha ) ) / f;
      return res;
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
