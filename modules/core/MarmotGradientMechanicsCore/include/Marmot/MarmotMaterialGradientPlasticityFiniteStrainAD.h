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
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGradientPlasticityFiniteStrain.h"
#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>
#include <functional>

template < int nYieldSurfaces >
class MarmotMaterialGradientPlasticityFiniteStrainAD
  : public MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces > {

public:
  MarmotMaterialGradientPlasticityFiniteStrainAD( const double* matProperties_,
                                                  int           nMaterialProperties_,
                                                  int           materialNumber_ )
    : MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces >( matProperties_,
                                                                      nMaterialProperties_,
                                                                      materialNumber_ )
  {
  }

  struct DeformationAD {
    Fastor::Tensor< autodiff::dual, 3, 3 > F;
  };

  struct incrementAD {
    DeformationAD                                   deformation;
    Eigen::Vector< autodiff::dual, nYieldSurfaces > dLambda;
    Eigen::Vector< autodiff::dual, nYieldSurfaces > laplaceDLambda;
    double                                          time;
    double                                          dT;
  };

  struct responseAD {
    Fastor::Tensor< autodiff::dual, 3, 3 >          tau;
    Eigen::Vector< autodiff::dual, nYieldSurfaces > f;
    double*                                         stateVars;
    autodiff::dual                                  elasticEnergyDensity;
    autodiff::dual                                  dissipation;

    responseAD()
      : tau( 0.0 ),
        f( Eigen::Vector< autodiff::dual, nYieldSurfaces >::Zero() ),
        stateVars( nullptr ),
        elasticEnergyDensity( 0.0 ),
        dissipation( 0.0 )
    {
    }

    responseAD( const typename MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces >::response& res )
      : tau( Marmot::makeDual( res.tau ) ),
        stateVars( res.stateVars ),
        elasticEnergyDensity( res.elasticEnergyDensity ),
        dissipation( res.dissipation )
    {
      for ( int i = 0; i < nYieldSurfaces; ++i )
        f( i ) = res.f( i );
    }
  };

  void computeStress(
    typename MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces >::response&        res,
    typename MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces >::tangents&        tan,
    const typename MarmotMaterialGradientPlasticityFiniteStrain< nYieldSurfaces >::increment& inc ) const override
  {
    using scalar = autodiff::dual;

    Eigen::Map< Eigen::VectorXd > stateVars( res.stateVars, this->getNumberOfRequiredStateVars() );
    const Eigen::VectorXd         stateVarsOld = stateVars;

    auto evaluateAD = [&]( int var_type, int idx ) {
      stateVars = stateVarsOld; // reset state

      DeformationAD defAD;
      defAD.F = Marmot::makeDual( inc.deformation.F );

      Eigen::Vector< scalar, nYieldSurfaces > dLambdaAD;
      Eigen::Vector< scalar, nYieldSurfaces > laplaceDLambdaAD;
      for ( int i = 0; i < nYieldSurfaces; ++i ) {
        dLambdaAD( i )        = inc.dLambda( i );
        laplaceDLambdaAD( i ) = inc.laplaceDLambda( i );
      }

      // Seed
      if ( var_type == 0 ) {
        defAD.F.data()[idx].grad = 1.0;
      }
      else if ( var_type == 1 ) {
        dLambdaAD( idx ).grad = 1.0;
      }
      else if ( var_type == 2 ) {
        laplaceDLambdaAD( idx ).grad = 1.0;
      }

      responseAD resAD;
      resAD.tau                  = Marmot::makeDual( res.tau );
      resAD.elasticEnergyDensity = res.elasticEnergyDensity;
      resAD.dissipation          = res.dissipation;
      resAD.stateVars            = res.stateVars;

      incrementAD incAD;
      incAD.deformation    = defAD;
      incAD.dLambda        = dLambdaAD;
      incAD.laplaceDLambda = laplaceDLambdaAD;
      incAD.time           = inc.time;
      incAD.dT             = inc.dT;

      this->computeStressAD( resAD, incAD );

      return resAD;
    };

    // 0. Base evaluation to get primal values
    {
      auto resAD = evaluateAD( -1, -1 );
      res.tau    = Marmot::makeReal( resAD.tau );
      for ( int i = 0; i < nYieldSurfaces; ++i ) {
        res.f( i ) = double( resAD.f( i ) );
      }
      res.elasticEnergyDensity = double( resAD.elasticEnergyDensity );
      res.dissipation          = double( resAD.dissipation );
    }

    // 1. Differentiate w.r.t F
    for ( int j = 0; j < 9; ++j ) {
      auto resAD = evaluateAD( 0, j );
      for ( int i = 0; i < 9; ++i ) {
        tan.dTau_dF.data()[i * 9 + j] = resAD.tau.data()[i].grad;
      }
      for ( int i = 0; i < nYieldSurfaces; ++i ) {
        tan.df_dF( i, j ) = resAD.f( i ).grad;
      }
    }

    // 2. Differentiate w.r.t dLambda
    for ( int j = 0; j < nYieldSurfaces; ++j ) {
      auto resAD = evaluateAD( 1, j );
      for ( int i = 0; i < 9; ++i ) {
        tan.dTau_ddLambda( i, j ) = resAD.tau.data()[i].grad;
      }
      for ( int i = 0; i < nYieldSurfaces; ++i ) {
        tan.dFddLambda( i, j ) = resAD.f( i ).grad;
      }
    }

    // 3. Differentiate w.r.t laplaceDLambda
    for ( int j = 0; j < nYieldSurfaces; ++j ) {
      auto resAD = evaluateAD( 2, j );
      for ( int i = 0; i < 9; ++i ) {
        tan.dTau_ddLaplacian( i, j ) = resAD.tau.data()[i].grad;
      }
      for ( int i = 0; i < nYieldSurfaces; ++i ) {
        tan.dFddLaplacian( i, j ) = resAD.f( i ).grad;
      }
    }
  }

  virtual void computeStressAD( responseAD& res, const incrementAD& inc ) const = 0;
};
