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
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

using namespace Eigen;
using namespace Marmot;

template < int nNonlocalVariables >
class MarmotMaterialGeneralGradientEnhancedHypoElastic {

protected:
  const double* materialProperties;
  const int     nMaterialProperties;

public:
  const int materialNumber;
  MarmotMaterialGeneralGradientEnhancedHypoElastic( const double* matProperties_,
                                                    int           nMaterialProperties_,
                                                    int           materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }
  struct increment {
    Marmot::Vector6d                            dStrain;
    Eigen::Vector< double, nNonlocalVariables > K;
    Eigen::Vector< double, nNonlocalVariables > dK;
    double                                      time;
    double                                      dT;
  };

  struct response {
    Marmot::Vector6d                            stress;
    Eigen::Vector< double, nNonlocalVariables > KLocal;
    Eigen::Vector< double, nNonlocalVariables > c;
    double*                                     stateVars;
  };

  struct tangents {
    Marmot::Matrix6d                               dStressddStrain = Marmot::Matrix6d::Zero();
    Eigen::Matrix< double, 6, nNonlocalVariables > dStressddK = Eigen::Matrix< double, 6, nNonlocalVariables >::Zero();
    Eigen::Matrix< double, nNonlocalVariables, 6 >
      dKLocalddStrain = Eigen::Matrix< double, nNonlocalVariables, 6 >::Zero();
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      dKLocalddK = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      dcddK = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      d2cddK2 = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
  };

  MarmotStateLayoutDynamic stateLayout;

  virtual void computeStress( response& res, tangents& tan, const increment& inc ) const = 0;

  virtual void computePlaneStress( response& res, tangents& tan, const increment& inc ) const
  {
    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    Map< VectorXd > stateVars( res.stateVars, stateLayout.totalSize() );

    VectorXd  stateVarsOld = stateVars;
    response  resTemp      = { res.stress, res.KLocal, res.c, res.stateVars };
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
      resTemp = { res.stress, res.KLocal, res.c, res.stateVars };
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

    res = { resTemp.stress, resTemp.KLocal, resTemp.c, resTemp.stateVars };
  }

  /**
   * @brief Initialize the layout of the state variables.
   *
   * This method has to be implemented in derived classes.
   * @warning This method has to be called in the constructor of the derived class.
   */
  virtual void initializeStateLayout() = 0;

  /**
   * @brief Get a view to the state variables.
   * @param stateName Name of the state variable
   * @param stateVars Pointer to the state variable array
   * @return StatView to access the state variable
   */
  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  /**
   * @brief Get the total number of required state variables.
   * @return Total number of required state variables
   */
  int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }
  /**
   * @brief Initialize the state variables at a material point.
   * @param stateVars Pointer to the state variable array
   * @param nStateVars Number of state variables
   *
   * @note The default implementation initializes all state variables to zero.
   */
  virtual void initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }
  }

  virtual double getDensity() { return -1; }
};
