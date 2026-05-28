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
 * Modified for C0-Continuous Penalty-Enhanced Gradient Plasticity
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
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"

/**
 * @brief Base class for C0-continuous penalty-enhanced gradient plasticity material models.
 *
 * @details This abstract class defines the interface for gradient-enhanced plasticity material
 * models using a penalty-based formulation. The formulation introduces a nonlocal cumulative
 * plastic strain field \f$ \bar{\kappa} \f$ that is C0-continuous and interpolated with the
 * same shape functions as the displacement field.
 *
 * The penalty method enforces the constraint \f$ \bar{\kappa} \approx \kappa \f$ (local
 * cumulative plastic strain) through a penalty parameter \f$ \beta \f$. The additional
 * balance equation reads:
 * \f[
 *   \beta \left( \bar{\kappa} - \kappa \right) - l^2 \nabla^2 \bar{\kappa} = 0
 * \f]
 * where \f$ l \f$ is the internal length scale parameter.
 *
 * The material model must provide the stress response as a function of the strain increment
 * and the nonlocal cumulative plastic strain, as well as the local cumulative plastic strain
 * and the consistent algorithmic tangent operators.
 */
class MarmotMaterialC0PenaltyGradientPlasticity {

protected:
  /// @brief Pointer to the array of material properties.
  const double* materialProperties;
  /// @brief Number of material properties.
  const int nMaterialProperties;

public:
  /// @brief Material number (identifier for the material).
  const int materialNumber;

  /**
   * @brief Constructor.
   * @param matProperties_ Pointer to the array of material properties.
   * @param nMaterialProperties_ Number of material properties.
   * @param materialNumber_ Material number (identifier).
   */
  MarmotMaterialC0PenaltyGradientPlasticity( const double* matProperties_,
                                             int           nMaterialProperties_,
                                             int           materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  /// @brief Struct to hold the strain increment and nonlocal variable information.
  struct increment {
    Marmot::Vector6d dStrain;  ///< Increment of the strain tensor in Voigt notation
    double           kappaNL;  ///< Nonlocal cumulative plastic strain at the current increment (total)
    double           dKappaNL; ///< Increment of the nonlocal cumulative plastic strain
    double           time;     ///< Current time
    double           dT;       ///< Time increment
  };

  /// @brief Struct to hold the material response information.
  struct response {
    Marmot::Vector6d stress;     ///< Stress tensor in Voigt notation
    double           kappaLocal; ///< Local cumulative plastic strain at the current increment
    double*          stateVars;  ///< Pointer to the array of state variables
  };

  /// @brief Struct to hold the algorithmic tangent matrices for the material model.
  struct tangents {
    /// @brief Algorithmic tangent: dStress / dStrain
    Marmot::Matrix6d dStressDDStrain = Marmot::Matrix6d::Zero();
    /// @brief Algorithmic tangent: dStress / dKappaNL (6x1)
    Marmot::Vector6d dStressDDKappaNL = Marmot::Vector6d::Zero();
    /// @brief Algorithmic tangent: dKappaLocal / dStrain (1x6 stored as row vector)
    Eigen::Matrix< double, 1, 6 > dKappaLocalDDStrain = Eigen::Matrix< double, 1, 6 >::Zero();
    /// @brief Algorithmic tangent: dKappaLocal / dKappaNL (scalar)
    double dKappaLocalDDKappaNL = 0.0;
  };

  /// @brief Layout of the state variables for the material model.
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @brief Compute the stress response of the material model.
   * @param[in,out] res Response struct to be filled with computed stress and local plastic strain.
   * @param[in,out] tan Tangents struct to be filled with algorithmic tangent matrices.
   * @param[in] inc Increment struct with strain increment and nonlocal plastic strain.
   */
  virtual void computeStress( response& res, tangents& tan, const increment& inc ) const = 0;

  /**
   * @brief Compute the plane stress response of the material model.
   * @param[in,out] res Response struct.
   * @param[in,out] tan Tangents struct for plane stress.
   * @param[in] inc Increment struct.
   */
  virtual void computePlaneStress( response& res, tangents& tan, const increment& inc ) const
  {
    using namespace Marmot;
    using namespace Eigen;

    Map< VectorXd > stateVars( res.stateVars, stateLayout.totalSize() );

    VectorXd  stateVarsOld = stateVars;
    response  resTemp      = { res.stress, res.kappaLocal, res.stateVars };
    increment incTemp      = { inc.dStrain, inc.kappaNL, inc.dKappaNL, inc.time, inc.dT };

    double residual          = 1;
    double tangentCompliance = 1.;
    // assumption of isochoric deformation for initial guess
    double strainIncrement = ( -incTemp.dStrain( 0 ) - incTemp.dStrain( 1 ) );
    incTemp.dStrain( 2 )   = strainIncrement;

    int planeStressCount = 1;
    while ( true ) {

      resTemp   = { res.stress, res.kappaLocal, res.stateVars };
      stateVars = stateVarsOld;
      computeStress( resTemp, tan, incTemp );

      residual = std::abs( resTemp.stress.array().abs()[2] / std::max( resTemp.stress.array().abs().maxCoeff(), 1. ) );
      if ( ( residual < 1e-10 && std::abs( strainIncrement ) < 1e-8 ) || ( planeStressCount > 7 && residual < 1e-5 ) ) {
        break;
      }

      tangentCompliance = 1. / tan.dStressDDStrain( 2, 2 );
      if ( Math::isNaN( tangentCompliance ) || std::abs( tangentCompliance ) > 1e10 ) {
        tangentCompliance = 1e10;
      }

      strainIncrement = -resTemp.stress( 2 ) * tangentCompliance;
      incTemp.dStrain( 2 ) += strainIncrement;

      planeStressCount += 1;
      if ( planeStressCount > 13 ) {
        throw Marmot::StressUpdateFailed( "PlaneStressWrapper requires cutback" );
      }
    }

    res = { resTemp.stress, resTemp.kappaLocal, resTemp.stateVars };
  }

  /**
   * @brief Get a view to the state variables.
   * @param stateName Name of the state variable
   * @param stateVars Pointer to the state variable array
   * @return StateView to access the state variable
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
   */
  virtual void initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }
  }

  /**
   * @brief Get the density of the material.
   * @return Density of the material
   */
  virtual double getDensity( const double* stateVars ) const = 0;

  virtual ~MarmotMaterialC0PenaltyGradientPlasticity() = default;
};
