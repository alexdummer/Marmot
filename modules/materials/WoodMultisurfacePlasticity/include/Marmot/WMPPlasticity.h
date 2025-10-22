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

  class WMPPlasticity {

    using Diag6d   = Eigen::DiagonalMatrix< double, 6 >;
    using Vector7d = Eigen::Vector< double, 7 >;

  public:
    struct ReturnMappingFailedException : std::exception {};

    enum failureMode {
      radialTension,
      radialCompression,
      tangentialTension,
      tangentialCompression,
      longitudinalTension,
      longitudinalCompression,
      shear
    };

    struct hp {
      double kappa;
      double eta;
      double zeta;
      double alphaD;
      double alphaMax;
    };

    struct hardening {
      hp rt; // radial tension
      hp rc; // radial compression
      hp tt; // tangential tension
      hp tc; // tangential compression
      hp lt; // longitudinal tension
      hp lc; // longitudinal compression
      hp s;  // longitudinal compression
    };

    struct strength {
      double rt;  // radial tension
      double rc;  // radial compression
      double tt;  // tangential tension
      double tc;  // tangential compression
      double lt;  // longitudinal tension
      double lc;  // longitudinal compression
      double rts; // radial-tangential shear
      double tls; // tangential-longitudinal shear
      double rls; // radial-longitudinal shear
    };

    struct parameters {
      strength  f;
      hardening h;
    };

    struct MaterialState {
      Marmot::Vector6d stress;
      Vector7d         alpha;
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

    WMPPlasticity( const parameters& );

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
    constexpr static int nInternalVariables = 7;
    constexpr static int nYieldSurfaces     = 7;

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

    const parameters mp;

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

    std::map< failureMode, hp > hardeningParameters{ { failureMode::radialTension, mp.h.rt },
                                                     { failureMode::radialCompression, mp.h.rc },
                                                     { failureMode::tangentialTension, mp.h.tt },
                                                     { failureMode::tangentialCompression, mp.h.tc },
                                                     { failureMode::longitudinalTension, mp.h.lt },
                                                     { failureMode::longitudinalCompression, mp.h.lc },
                                                     { failureMode::shear, mp.h.s } };

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
          const Vector6t< T > plasticFlowDirection = a[fm] + 2. * b[fm] * stress;
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
    std::map< failureMode, Vector6d > a{ { failureMode::radialTension, { 1. / mp.f.rt, 0., 0., 0., 0., 0. } },
                                         { failureMode::radialCompression, { -1. / mp.f.rc, 0., 0., 0., 0., 0. } },
                                         { failureMode::tangentialTension, { 0., 1. / mp.f.tt, 0., 0., 0., 0. } },
                                         { failureMode::tangentialCompression, { 0., -1. / mp.f.tc, 0., 0., 0., 0. } },
                                         { failureMode::longitudinalTension, { 0., 0., 1. / mp.f.lt, 0., 0., 0. } },
                                         { failureMode::longitudinalCompression,
                                           { 0., 0., -1. / mp.f.lc, 0., 0., 0. } },
                                         { failureMode::shear,
                                           { 0.20 / mp.f.rc, 0.2 / mp.f.tc, 0.05 / mp.f.lc, 0., 0., 0. } } };
    // map for strength tensor b
    std::map< failureMode, Diag6d > b{ { failureMode::radialTension,
                                         { 0.0000,
                                           0.5000 / ( mp.f.tc * mp.f.tc ),
                                           0.2500 / ( mp.f.lc * mp.f.lt ),
                                           0.5000 / ( mp.f.rts * mp.f.rts ),
                                           0.3300 / ( mp.f.rls * mp.f.rls ),
                                           0.3300 / ( mp.f.tls * mp.f.tls ) } },
                                       { failureMode::radialCompression,
                                         { 0.0000,
                                           0.4000 / ( mp.f.tc * mp.f.tc ),
                                           0.2500 / ( mp.f.lc * mp.f.lt ),
                                           0.4000 / ( mp.f.rts * mp.f.rts ),
                                           0.3300 / ( mp.f.rls * mp.f.rls ),
                                           0.3300 / ( mp.f.tls * mp.f.tls ) } },
                                       { failureMode::tangentialTension,
                                         { 0.5000 / ( mp.f.rc * mp.f.rc ),
                                           0.0000,
                                           0.2500 / ( mp.f.lc * mp.f.lt ),
                                           0.4000 / ( mp.f.rts * mp.f.rts ),
                                           0.3300 / ( mp.f.rls * mp.f.rls ),
                                           0.3300 / ( mp.f.tls * mp.f.tls ) } },
                                       { failureMode::tangentialCompression,
                                         { 0.4000 / ( mp.f.rc * mp.f.rc ),
                                           0.0000,
                                           0.2500 / ( mp.f.lc * mp.f.lt ),
                                           0.4000 / ( mp.f.rts * mp.f.rts ),
                                           0.3300 / ( mp.f.rls * mp.f.rls ),
                                           0.3300 / ( mp.f.tls * mp.f.tls ) } },
                                       { failureMode::longitudinalTension,
                                         { 0.2000 / ( mp.f.rc * mp.f.rc ),
                                           0.2000 / ( mp.f.tc * mp.f.tc ),
                                           0.0000,
                                           0.1000 / ( mp.f.rts * mp.f.rts ),
                                           0.1000 / ( mp.f.rls * mp.f.rls ),
                                           0.1000 / ( mp.f.tls * mp.f.tls ) } },
                                       { failureMode::longitudinalCompression,
                                         { 0.3300 / ( mp.f.rc * mp.f.rc ),
                                           0.3300 / ( mp.f.tc * mp.f.tc ),
                                           0.0000,
                                           0.2500 / ( mp.f.rts * mp.f.rts ),
                                           0.2500 / ( mp.f.rls * mp.f.rls ),
                                           0.2500 / ( mp.f.tls * mp.f.tls ) } },
                                       { failureMode::shear,
                                         { 0.2500 / ( mp.f.rc * mp.f.rc ),
                                           0.2500 / ( mp.f.tc * mp.f.tc ),
                                           0.1000 / ( mp.f.lc * mp.f.lt ),
                                           1.0000 / ( mp.f.rts * mp.f.rts ),
                                           1.0000 / ( mp.f.rls * mp.f.rls ),
                                           1.0000 / ( mp.f.tls * mp.f.tls ) } } };

    template < typename T >
    T yieldFunction( const Marmot::Vector6t< T >& stress, const failureMode fm, const T alpha )
    {
      const T res = a[fm].dot( stress ) + stress.dot( b[fm] * stress ) + q( alpha, fm ) - 1.0;
      return res;
    }

    template < typename T >
    T q( const T alpha, const failureMode fm )
    {

      const auto h = hardeningParameters[fm];
      // hardening/softening
      T res = ( 1. - h.kappa ) * ( 1. - exp( -h.eta * alpha ) );
      if ( double( alpha ) > h.alphaD ) {
        res -= h.zeta * pow( alpha - h.alphaD, 2 ) / ( h.alphaMax - alpha );
      }
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
