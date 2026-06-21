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
#include "Marmot/MarmotStateHelpers.h"
#include <Eigen/Core>
#include <Fastor/tensor/Tensor.h>
#include <algorithm>
#include <cmath>
#include <vector>

template < int nYieldSurfaces >
class MarmotMaterialGradientPlasticityFiniteStrain {

protected:
  const double* materialProperties;
  const int     nMaterialProperties;

public:
  const int materialNumber;

  MarmotMaterialGradientPlasticityFiniteStrain( const double* matProperties_,
                                                int           nMaterialProperties_,
                                                int           materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  virtual ~MarmotMaterialGradientPlasticityFiniteStrain() = default;

  struct Deformation {
    Fastor::Tensor< double, 3, 3 > F;
  };

  struct increment {
    Deformation                             deformation;
    Eigen::Vector< double, nYieldSurfaces > dLambda;
    Eigen::Vector< double, nYieldSurfaces > laplaceDLambda;
    double                                  time;
    double                                  dT;
  };

  struct response {
    Fastor::Tensor< double, 3, 3 >          tau;
    Eigen::Vector< double, nYieldSurfaces > f;
    double*                                 stateVars;
    double                                  elasticEnergyDensity;
    double                                  dissipation;
  };

  struct tangents {
    Fastor::Tensor< double, 3, 3, 3, 3 >       dTau_dF          = Fastor::Tensor< double, 3, 3, 3, 3 >( 0.0 );
    Eigen::Matrix< double, 9, nYieldSurfaces > dTau_ddLambda    = Eigen::Matrix< double, 9, nYieldSurfaces >::Zero();
    Eigen::Matrix< double, 9, nYieldSurfaces > dTau_ddLaplacian = Eigen::Matrix< double, 9, nYieldSurfaces >::Zero();
    Eigen::Matrix< double, nYieldSurfaces, 9 > df_dF            = Eigen::Matrix< double, nYieldSurfaces, 9 >::Zero();
    Eigen::Matrix< double, nYieldSurfaces, nYieldSurfaces >
      dFddLambda = Eigen::Matrix< double, nYieldSurfaces, nYieldSurfaces >::Zero();
    Eigen::Matrix< double, nYieldSurfaces, nYieldSurfaces >
      dFddLaplacian = Eigen::Matrix< double, nYieldSurfaces, nYieldSurfaces >::Zero();
  };

  MarmotStateLayoutDynamic stateLayout;

  virtual void computeStress( response& res, tangents& tan, const increment& inc ) const = 0;

  virtual void computeStressExplicit( response& res, const increment& inc ) const
  {
    tangents tan;
    computeStress( res, tan, inc );
  }

  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }
  }

  virtual double getDensity( const double* stateVars ) const = 0;

  virtual std::vector< double > getNonlocalViscosity( const double* stateVars ) const = 0;
};
