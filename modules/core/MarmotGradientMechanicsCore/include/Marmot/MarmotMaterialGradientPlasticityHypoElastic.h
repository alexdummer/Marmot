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
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"
#include <algorithm>
#include <cmath>
#include <vector>

template < int nYieldSurfaces >
class MarmotMaterialGradientPlasticityHypoElastic {

protected:
  /// @brief Pointer to the array of material properties.
  const double* materialProperties;
  /// @brief Number of material properties.
  const int nMaterialProperties;

public:
  /// @brief Material number (identifier for the material).
  const int materialNumber;

  /**
   * @brief Constructor for the gradient-plastic hypoelastic material model.
   * @param matProperties_ Pointer to the array of material properties.
   * @param nMaterialProperties_ Number of material properties.
   * @param materialNumber_ Material number (identifier for the material).
   */
  MarmotMaterialGradientPlasticityHypoElastic( const double* matProperties_,
                                               int           nMaterialProperties_,
                                               int           materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  /// @brief Struct to hold the increment information.
  struct increment {
    Marmot::Vector6d                        dStrain;        ///< Increment of the strain tensor in Voigt notation
    Eigen::Vector< double, nYieldSurfaces > dLambda;        ///< Increment of the plastic multipliers
    Eigen::Vector< double, nYieldSurfaces > laplaceDLambda; ///< Laplacian of the increment of the plastic multipliers
    double                                  time;           ///< Current time
    double                                  dT;             ///< Time increment
  };

  /// @brief Struct to hold the material response information.
  struct response {
    Marmot::Vector6d                        stress;    ///< Stress tensor in Voigt notation
    Eigen::Vector< double, nYieldSurfaces > f;         ///< Local yield function values at the current increment
    double*                                 stateVars; ///< Pointer to the array of state variables
    double                                  elasticEnergyDensity; ///< Elastic strain energy density
    double                                  dissipation;          ///< Dissipation (if applicable)
  };

  /// @brief Struct to hold the algorithmic tangent matrices for the material model.
  struct tangents {
    /// @brief Algorithmic tangent matrix relating the increment of stress to the increment of strain.
    Marmot::Matrix6d dStressddStrain = Marmot::Matrix6d::Zero();
    /// @brief Algorithmic tangent matrix relating the increment of stress to the increment of the plastic multipliers.
    Eigen::Matrix< double, 6, nYieldSurfaces > dStressddLambda = Eigen::Matrix< double, 6, nYieldSurfaces >::Zero();
    /// @brief Algorithmic tangent relating the increment of stress to the increment of the laplacian of the plastic
    /// multipliers.
    Eigen::Matrix< double, 6, nYieldSurfaces > dStressddLaplacian = Eigen::Matrix< double, 6, nYieldSurfaces >::Zero();
    /// @brief Algorithmic tangent matrix relating the increment of the local yield function values to the increment of
    /// strain.
    Eigen::Matrix< double, nYieldSurfaces, 6 > dFddStrain = Eigen::Matrix< double, nYieldSurfaces, 6 >::Zero();
    /// @brief Algorithmic tangent matrix relating the increment of the local yield function values to the increment of
    /// the plastic multipliers.
    Eigen::Matrix< double, nYieldSurfaces, nYieldSurfaces >
      dFddLambda = Eigen::Matrix< double, nYieldSurfaces, nYieldSurfaces >::Zero();
    /// @brief Algorithmic tangent matrix relating the increment of the local yield function values to the increment of
    /// the laplacian of the plastic multipliers.
    Eigen::Matrix< double, nYieldSurfaces, nYieldSurfaces >
      dFddLaplacian = Eigen::Matrix< double, nYieldSurfaces, nYieldSurfaces >::Zero();
  };

  /**
   * @brief Layout of the state variables for the material model.
   * @note Must be defined in derived classes to specify the structure and organization of the state variables used in
   * the material model.
   */
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @brief Compute the stress response of the material model including the algorithmic tangents.
   * @param[in,out] res Reference to the response struct to be filled with the computed stress and other response
   * variables.
   * @param[in,out] tan Reference to the tangents struct to be filled with the computed algorithmic tangent matrices.
   * @param[in] inc Reference to the increment struct containing the strain increment, nonlocal variable increments, and
   * time information.
   *
   * This method must be implemented in derived classes to compute the stress response based on the specific material
   * model.
   */
  virtual void computeStress( response& res, tangents& tan, const increment& inc ) const = 0;

  /**
   * @brief Compute the stress response of the material model.
   * @param[in,out] res Reference to the response struct to be filled with the computed stress and other response
   * variables.
   * @param[in] inc Reference to the increment struct containing the strain increment, nonlocal variable increments, and
   * time information.
   *
   * This method can be used when only the stress response is needed without the tangents, e.g., for explicit dynamics.
   * The default implementation calls the `computeStress` method and ignores the tangents, which can be computationally
   * expensive to compute.
   */
  virtual void computeStressExplicit( response& res, const increment& inc ) const
  {
    tangents tan;
    computeStress( res, tan, inc );
  }

  /**
   * @brief Compute the plane stress response of the material model.
   * @param[in,out] res Reference to the response struct to be filled with the computed plane stress and other response
   * variables.
   * @param[in,out] tan Reference to the tangents struct to be filled with the computed algorithmic tangent matrices for
   * plane stress.
   * @param[in] inc Reference to the increment struct containing the strain increment, nonlocal variable increments, and
   * time information.
   *
   * This method provides a default implementation for computing the plane stress response using an iterative approach.
   * It assumes an initial guess of isochoric deformation and iteratively corrects the out-of-plane strain until
   * convergence is achieved. Derived classes can override this method if a different approach for plane stress
   * computation is desired.
   */
  virtual void computePlaneStress( response& res, tangents& tan, const increment& inc ) const
  {
    using namespace Marmot;
    using namespace Eigen;

    Map< VectorXd > stateVars( res.stateVars, stateLayout.totalSize() );

    VectorXd  stateVarsOld = stateVars;
    response  resTemp      = { res.stress, res.f, res.stateVars, res.elasticEnergyDensity, res.dissipation };
    increment incTemp      = { inc.dStrain, inc.dLambda, inc.laplaceDLambda, inc.time, inc.dT };

    double residual          = 1;
    double tangentCompliance = 1.;
    // assumption of isochoric deformation for initial guess
    double strainIncrement = ( -incTemp.dStrain( 0 ) - incTemp.dStrain( 1 ) );
    incTemp.dStrain( 2 )   = strainIncrement;

    int planeStressCount = 1;
    while ( true ) {

      // set old response
      resTemp = { .stress               = res.stress,
                  .f                    = res.f,
                  .stateVars            = res.stateVars,
                  .elasticEnergyDensity = res.elasticEnergyDensity,
                  .dissipation          = res.dissipation };
      // set old state variables
      stateVars = stateVarsOld;
      // compute stress
      computeStress( resTemp, tan, incTemp );

      // evauate residual
      residual = std::abs( resTemp.stress.array().abs()[2] / std::max( resTemp.stress.array().abs().maxCoeff(), 1. ) );
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
        throw Marmot::StressUpdateFailed( "PlaneStressWrapper requires cutback" );
      }
    }

    res = resTemp;
  }

  /**
   * @brief Compute the plane stress response of the material model without computing tangents.
   * @param[in,out] res Reference to the response struct to be filled with the computed plane stress and other response
   * variables.
   * @param[in] inc Reference to the increment struct containing the strain increment, nonlocal variable increments, and
   * time information.
   *
   * This method can be used when only the plane stress response is needed without the tangents, e.g., for explicit
   * dynamics. The default implementation calls the `computePlaneStress` method ignoring the tangents.
   * This method maybe overridden in derived classes if a different approach for plane stress computation is desired.
   */
  virtual void computePlaneStressExplicit( response& res, const increment& inc ) const
  {
    tangents tan;
    computePlaneStress( res, tan, inc );
  }

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
  /**
   * @brief Get the density of the material.
   * @param stateVars Pointer to the array of state variables
   * @return Density of the material
   *
   * This method must be implemented in derived classes to return the density of the material, which is required for
   * dynamic analyses.
   */
  virtual double getDensity( const double* stateVars ) const = 0;

  /**
   * @brief Get the maximum wave speed for the current response state.
   * @param currentResponse Current response state
   * @return Maximum wave speed
   * @details The default implementation computes the 3D algorithmic tangent and returns
   *          `sqrt(max(C_ii) / rho)` with `C_ii` from the Voigt tangent diagonal entries.
   */
  virtual double getMaximumWaveSpeed( const response& currentResponse ) const
  {
    const int nStateVars = getNumberOfRequiredStateVars();

    std::vector< double > stateVarsCopy( nStateVars, 0.0 );
    if ( currentResponse.stateVars != nullptr && nStateVars > 0 ) {
      std::copy_n( currentResponse.stateVars, nStateVars, stateVarsCopy.begin() );
    }

    response responseCopy  = currentResponse;
    responseCopy.stateVars = stateVarsCopy.data();

    tangents  tan;
    increment inc;
    inc.dStrain        = Marmot::Vector6d::Zero();
    inc.dLambda        = Eigen::Vector< double, nYieldSurfaces >::Zero();
    inc.laplaceDLambda = Eigen::Vector< double, nYieldSurfaces >::Zero();
    inc.time           = 0.0;
    inc.dT             = 1.0;

    computeStress( responseCopy, tan, inc );

    const double maxStiffnessDiagonal = std::max( { tan.dStressddStrain( 0, 0 ),
                                                    tan.dStressddStrain( 1, 1 ),
                                                    tan.dStressddStrain( 2, 2 ),
                                                    tan.dStressddStrain( 3, 3 ),
                                                    tan.dStressddStrain( 4, 4 ),
                                                    tan.dStressddStrain( 5, 5 ) } );
    const double density              = getDensity( responseCopy.stateVars );

    return density > 0.0 ? std::sqrt( std::max( 0.0, maxStiffnessDiagonal ) / density ) : 0.0;
  }

  /**
   * @brief Get the nonlocal viscosity of the material.
   * @param stateVars Pointer to the array of state variables
   * @return Vector containing the nonlocal viscosity values for each nonlocal variable
   *
   * This method must be implemented in derived classes to return the nonlocal viscosity values, which are required for
   * dynamic analyses involving nonlocal variables. The returned vector should have a length equal to
   * `nYieldSurfaces`.
   */
  virtual std::vector< double > getNonlocalViscosity( const double* stateVars ) const = 0;
};
