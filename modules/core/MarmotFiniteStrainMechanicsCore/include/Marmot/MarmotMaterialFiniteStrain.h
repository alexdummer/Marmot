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
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotStateHelpers.h"
#include <Fastor/tensor/Tensor.h>

class MarmotMaterialFiniteStrain {

  /**
   * @class MarmotMaterialFiniteStrain
   * @brief Abstract basic class for mechanical materials in the finite strain regime
   */
protected:
  const double* materialProperties;
  const int     nMaterialProperties;

public:
  const int materialNumber;
  MarmotMaterialFiniteStrain( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
    : materialProperties( matProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  /// Default destructor
  virtual ~MarmotMaterialFiniteStrain() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @struct ConstitutiveResponse
   * @brief Constitutive response of a material at given state.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   *  Contains stress, density and elastic energy density.
   */
  template < int nDim >
  struct ConstitutiveResponse {
    Fastor::Tensor< double, nDim, nDim > tau;                  ///< Kirchhoff stress
    double                               rho;                  ///< mass density
    double                               elasticEnergyDensity; ///< elastic energy per unit volume
    double*                              stateVars;            ///< pointer to state variables
  };

  /**
   * @struct AlgorithmicModuli
   * @brief Algorithmic tangent moduli of a material.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   * Contains the algorithmic tangent moduli \f$\frac{\partial \boldsymbol{\tau}}{\partial \boldsymbol{F}}\f$
   * with respect to the deformation gradient \f$\boldsymbol{F}\f$.
   * */
  template < int nDim >
  struct AlgorithmicModuli {
    Fastor::Tensor< double, nDim, nDim, nDim, nDim > dTau_dF; ///< tangent operator w.r.t. deformation gradient
  };

  /**
   * @struct Deformation
   * @brief Represents the deformation state of a material.
   *
   * This struct holds the deformation gradient \f$\boldsymbol{F}\f$ which describes the local deformation of a
   * material.
   */
  template < int nDim >
  struct Deformation {
    Fastor::Tensor< double, nDim, nDim > F; ///< deformation gradient
  };

  /**
   * @struct TimeIncrement
   * @brief Represents a time increment in a simulation.
   *
   * This struct holds information about the current time
   * and the time step size (dT).
   */
  struct TimeIncrement {
    const double time; ///< time at the beginning of the increment
    const double dT;   ///< size of the time increment
  };

  /**
   * @struct StateSensitivities
   * @brief Sensitivity matrices required for analytical substepping.
   * * All matrices represent derivatives of flattened arrays (Row-Major Fastor layout).
   */
  struct StateSensitivities {
    /** * @brief Jacobian of State_new w.r.t Deformation Gradient F_new.
     * Size: nStateVars x 9
     */
    Eigen::MatrixXd dState_dF;

    /** * @brief Jacobian of State_new w.r.t State_old.
     * Size: nStateVars x nStateVars
     */
    Eigen::MatrixXd dState_dStateOld;

    /** * @brief Jacobian of Stress w.r.t State_old.
     * Size: 9 x nStateVars
     */
    Eigen::MatrixXd dStress_dStateOld;
  };

  /**
   * @brief Updates the material state.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] tangents AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * Computes the Kirchhoff \f$\boldsymbol{\tau}\f$ stress
   * from the deformation gradient \f$\boldsymbol{F}\f$ and the
   * time increment \f$\Delta t\f$.
   * It further updates the mass density \f$\rho\f$ and the elastic energy density.
   * Additionally, computes the algorithmic tangent moduli
   * \f$\frac{\partial \boldsymbol{\tau}}{\partial \boldsymbol{F}}\f$.
   *
   * */
  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const Deformation< 3 >&,
                              const TimeIncrement& ) const = 0;

  /**
   * @brief Extended stress update computing sensitivities for substepping.
   * * The default implementation uses Central Finite Differences.
   */
  virtual void computeStressWithSensitivities( ConstitutiveResponse< 3 >& response,
                                               AlgorithmicModuli< 3 >&    tangents,
                                               StateSensitivities&        sensitivities,
                                               const Deformation< 3 >&    deformation,
                                               const TimeIncrement&       timeIncrement ) const
  {
    using namespace Eigen;
    using namespace Marmot::NumericalAlgorithms::Differentiation;

    // 1. Compute the baseline stress and algorithmic tangent (dTau/dF)
    // This updates response.stateVars to the *new* state.
    // We must be careful: to compute derivatives w.r.t Old State, we need the Old State.
    // Since computeStress updates in-place, we must assume response.stateVars contains
    // the OLD state when this function is called.

    int nState = getNumberOfRequiredStateVars();

    // Backup the OLD state (at t_n)
    std::vector< double > stateOld( nState );
    std::memcpy( stateOld.data(), response.stateVars, nState * sizeof( double ) );

    /* std::cout << "deformation.F: \n" << deformation.F << std::endl; */
    /* std::cout << "stateOld: \n" << Map<const VectorXd>(stateOld.data(), nState) << std::endl; */
    // Run the actual update once to get Stress and New State (at t_n+1)
    computeStress( response, tangents, deformation, timeIncrement );

    /* std::cout << "Base Stress after update: \n" << response.tau << std::endl; */

    // If there are no state variables (elastic), sensitivities are zero/empty
    if ( nState == 0 )
      return;

    // --- Numerical Differentiation Setup ---

    // Prepare sensitivity matrices
    sensitivities.dState_dF.resize( nState, 9 );
    sensitivities.dState_dStateOld.resize( nState, nState );
    sensitivities.dStress_dStateOld.resize( 9, nState );

    // We need temporary storage for the FD routines to avoid corrupting the actual results
    // or the input pointers.

    // 2. Compute d(State_new) / d(F_new)
    // We perturb F, run computeStress (starting from fixed StateOld), and measure State_new.
    auto func_dState_dF = [&]( const VectorXd& F_vec ) -> VectorXd {
      Deformation< 3 > defPerturbed;
      std::memcpy( defPerturbed.F.data(), F_vec.data(), 9 * sizeof( double ) );

      ConstitutiveResponse< 3 > respPert;
      std::vector< double >     stateTemp = stateOld; // Copy from OLD
      respPert.stateVars                  = stateTemp.data();

      AlgorithmicModuli< 3 > tanPert; // dummy

      // Run update
      computeStress( respPert, tanPert, defPerturbed, timeIncrement );

      // Return NEW state
      Map< VectorXd > res( respPert.stateVars, nState );
      return res;
    };

    Map< const VectorXd > F_map( deformation.F.data(), 9 );
    sensitivities.dState_dF = forwardDifference( func_dState_dF, F_map );

    // print dState_dF if det is nan
    if ( sensitivities.dState_dF.hasNaN() ) {
      std::cout << "Computed dState/dF has NaN!" << std::endl;
      std::cout << "F: \n" << deformation.F << std::endl;
      std::cout << "StateOld: \n" << Map< const VectorXd >( stateOld.data(), nState ) << std::endl;
    }
    /* std::cout << "Computed dState/dF via FD: \n" << sensitivities.dState_dF << std::endl; */

    // 3. Compute d(State_new) / d(State_old) AND d(Stress) / d(State_old)
    // We perturb State_old, run computeStress (fixed F), and measure State_new and Stress.

    // Since centralDifference returns one matrix, we do this in two passes or combine outputs.
    // Let's do two passes for clarity, though it doubles the cost of this specific block.
    // Given the prompt constraints, we use the provided interface.

    // A. d(State_new) / d(State_old)
    auto func_dState_dStateOld = [&]( const VectorXd& StateOld_vec ) -> VectorXd {
      ConstitutiveResponse< 3 > respPert;
      std::vector< double >     stateTemp( nState );
      std::memcpy( stateTemp.data(), StateOld_vec.data(), nState * sizeof( double ) );
      respPert.stateVars = stateTemp.data();

      AlgorithmicModuli< 3 > tanPert; // dummy

      // Run update with perturbed old state
      computeStress( respPert, tanPert, deformation, timeIncrement );

      Map< VectorXd > res( respPert.stateVars, nState );
      return res;
    };

    Map< const VectorXd > StateOld_map( stateOld.data(), nState );
    sensitivities.dState_dStateOld = forwardDifference( func_dState_dStateOld, StateOld_map );

    // print dState_dF if det is nan
    if ( sensitivities.dState_dStateOld.hasNaN() ) {
      std::cout << "Computed dState/dStateOld has NaN!" << std::endl;
      std::cout << "F: \n" << deformation.F << std::endl;
      std::cout << "StateOld: \n" << Map< const VectorXd >( stateOld.data(), nState ) << std::endl;
    }
    // B. d(Stress) / d(State_old)
    auto func_dStress_dStateOld = [&]( const VectorXd& StateOld_vec ) -> VectorXd {
      ConstitutiveResponse< 3 > respPert;
      std::vector< double >     stateTemp( nState );
      std::memcpy( stateTemp.data(), StateOld_vec.data(), nState * sizeof( double ) );
      respPert.stateVars = stateTemp.data();

      AlgorithmicModuli< 3 > tanPert; // dummy

      computeStress( respPert, tanPert, deformation, timeIncrement );

      // Return Stress (flattened 9x1)
      Map< VectorXd > res( respPert.tau.data(), 9 );
      return res;
    };

    sensitivities.dStress_dStateOld = forwardDifference( func_dStress_dStateOld, StateOld_map );
    if ( sensitivities.dStress_dStateOld.hasNaN() ) {
      std::cout << "Computed dStress/dStateOld has NaN!" << std::endl;
      std::cout << "F: \n" << deformation.F << std::endl;
      std::cout << "StateOld: \n" << Map< const VectorXd >( stateOld.data(), nState ) << std::endl;
    }
  }

  /**
   * @brief Computes the Kirchhoff stress given the deformation, time increment, and eigen deformation.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] tangents AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   * @param[in] eigenDeformation Tuple representing eigen deformation in each spatial direction.
   *
   */
  virtual void computeStress( ConstitutiveResponse< 3 >&                  response,
                              AlgorithmicModuli< 3 >&                     tangents,
                              const Deformation< 3 >&                     deformation,
                              const TimeIncrement&                        timeIncrement,
                              const std::tuple< double, double, double >& eigenDeformation ) const;

  /**
   * @brief Compute stress under plane strain conditions.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * It uses the general 3D computeStress function for a plane strain Deformation.
   * The algorithmic tangent is modified according to plane strain conditions.
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    algorithmicModuli,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement ) const;
  /**
   * @brief Compute stress under plane strain conditions with eigen deformation.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   * @param[in] eigenDeformation Tuple representing eigen deformation in each spatial direction.
   *
   * It uses the general 3D computeStress function for a plane strain Deformation.
   * The algorithmic tangent is modified according to plane strain conditions.
   */
  virtual void computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                   AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                   const Deformation< 3 >&                     deformation,
                                   const TimeIncrement&                        timeIncrement,
                                   const std::tuple< double, double, double >& eigenDeformation ) const;
  /**
   * @brief Compute stress under plane stress conditions.
   * @param[inout] response ConstitutiveResponse instance
   * @param[out] algorithmicModuli AlgorithmicModuli instance
   * @param[in] deformation Deformation instance
   * @param[in] timeIncrement TimeIncrement instance
   *
   * It uses the general 3D computeStress function and iteratively finds the out-of-plane deformation.
   * The algorithmic tangent is modified according to plane stress conditions.
   */
  virtual void computePlaneStress( ConstitutiveResponse< 2 >& response,
                                   AlgorithmicModuli< 2 >&    algorithmicModuli,
                                   const Deformation< 2 >&    deformation,
                                   const TimeIncrement&       timeIncrement ) const;

  /**
   * @brief Find the eigen deformation that corresponds to a given eigen stress.
   * @param initialGuess Initial guess for the eigen deformation.
   * @param eigenStress Target eigen stress.
   * @return Eigen deformation that corresponds to the given eigen stress.
   *
   * This function iteratively finds the eigen deformation that corresponds to a given eigen stress.
   * This is used e.g. for geostatic stress initialization.
   */
  std::tuple< double, double, double > findEigenDeformationForEigenStress(
    const std::tuple< double, double, double >& initialGuess,
    const std::tuple< double, double, double >& eigenStress,
    double*                                     stateVars ) const;

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
};
