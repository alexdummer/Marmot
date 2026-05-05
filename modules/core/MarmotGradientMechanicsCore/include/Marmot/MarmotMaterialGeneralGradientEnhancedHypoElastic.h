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

/**
 * @brief Base class for general gradient-enhanced hypoelastic material models.
 * @details This class defines the interface for general gradient-enhanced hypoelastic material models,
 * including methods for computing stress, plane stress response, and managing state variables.
 * The template parameter `nNonlocalVariables` specifies the number of nonlocal variables used in the material model.
 *
 * In addition to the standard balance of linear momentum each nonlocal variable introduces an additional balance
 * equation, which is solved simultaneously with the balance of linear momentum. This balance equation is defined as:
 * \f[ \kappa_i - \nabla ( c( \kappa_i ) \nabla \kappa_i ) = f_i (\boldsymbol \varepsilon, \kappa_i ) \f]
 * where \f$ \kappa_i \f$ is the nonlocal variable, \f$ c( \kappa_i ) \f$ is the nonlocal interaction function, and \f$
 * f_i \f$ is the local driving variable for the nonlocal variable \f$ \kappa_i \f$. The nonlocal interaction function
 * \f$ g( \kappa_i ) \f$ defines the influence of the nonlocal variable on its own gradient and can be used to model
 * phenomena such as damage-dependent interactions. Also phase-field models can be implemented in this framework by
 * defining the nonlocal variable as the phase-field variable and the local driving variable as a function of the strain
 * tensor.
 */
template < int nNonlocalVariables >
class MarmotMaterialGeneralGradientEnhancedHypoElastic {

protected:
  /// @brief Pointer to the array of material properties.
  const double* materialProperties;
  /// @brief Number of material properties.
  const int nMaterialProperties;

public:
  /// @brief Material number (identifier for the material).
  const int materialNumber;

  /**
   * @brief Constructor for the general gradient-enhanced hypoelastic material model.
   * @param matProperties_ Pointer to the array of material properties.
   * @param nMaterialProperties_ Number of material properties.
   * @param materialNumber_ Material number (identifier for the material).
   */
  MarmotMaterialGeneralGradientEnhancedHypoElastic( const double* matProperties_,
                                                    int           nMaterialProperties_,
                                                    int           materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  /// @brief Struct to hold the increment information.
  struct increment {
    Marmot::Vector6d                            dStrain; ///< Increment of the strain tensor in Voigt notation
    Eigen::Vector< double, nNonlocalVariables > K;       ///< Nonlocal variables at the current increment
    Eigen::Vector< double, nNonlocalVariables > dK;      ///< Increment of the nonlocal variables
    double                                      time;    ///< Current time
    double                                      dT;      ///< Time increment
  };

  /// @brief Struct to hold the material response information.
  struct response {
    Marmot::Vector6d                            stress;    ///< Stress tensor in Voigt notation
    Eigen::Vector< double, nNonlocalVariables > KLocal;    ///< Local driving variables at the current increment
    Eigen::Vector< double, nNonlocalVariables > c;         ///< Nonlocal interaction parameters at the current increment
    double*                                     stateVars; ///< Pointer to the array of state variables
  };

  /// @brief Struct to hold the algorithmic tangent matrices for the material model.
  struct tangents {
    ///> Algorithmic tangent matrix relating the increment of stress to the increment of strain.
    Marmot::Matrix6d dStressddStrain = Marmot::Matrix6d::Zero();
    ///> Algorithmic tangent matrix relating the increment of stress to the increment of nonlocal variables.
    Eigen::Matrix< double, 6, nNonlocalVariables > dStressddK = Eigen::Matrix< double, 6, nNonlocalVariables >::Zero();
    Eigen::Matrix< double, nNonlocalVariables, 6 >
      ///> Algorithmic tangent matrix relating the increment of local driving variables to the increment of strain.
      dKLocalddStrain = Eigen::Matrix< double, nNonlocalVariables, 6 >::Zero();
    ///> Algorithmic tangent matrix relating the increment of local driving variables to the increment of nonlocal
    /// variables.
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      dKLocalddK = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
    ///> Algorithmic tangent matrix relating the increment of nonlocal interaction parameters to the increment of
    /// strain.
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      dcddK = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
    ///> Algorithmic tangent matrix relating the increment of nonlocal interaction parameters to the increment of
    /// nonlocal variables.
    Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >
      d2cddK2 = Eigen::Matrix< double, nNonlocalVariables, nNonlocalVariables >::Zero();
  };

  /// @brief Layout of the state variables for the material model.
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @brief Compute the stress response of the material model.
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
    using namespace Marmot::ContinuumMechanics::VoigtNotation;

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
