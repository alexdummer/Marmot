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
#include <Eigen/Core>
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"

/**
 * @brief Base class for finite strain materials using automatic differentiation.
 *
 * This class extends MarmotMaterialFiniteStrain by providing an interface for
 * computing stresses and algorithmic tangent moduli using automatic
 * differentiation. Derived material models implement computeStressAD() to
 * supply the AD-based constitutive response.
 */
class MarmotMaterialFiniteStrainAD : public MarmotMaterialFiniteStrain {

public:
  MarmotMaterialFiniteStrainAD( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
    : MarmotMaterialFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ )
  {
  }
  /**
   * @struct ConstitutiveResponseAD
   * @brief Constitutive response of a material at given state.
   * @tparam nDim Number of spatial dimensions (2 or 3).
   *
   *  Contains stress, density and elastic energy density.
   */
  template < int nDim >
  struct ConstitutiveResponseAD {
    Fastor::Tensor< autodiff::dual, nDim, nDim > tau;                  ///< Kirchhoff stress
    double                                       rho;                  ///< mass density
    double                                       elasticEnergyDensity; ///< elastic energy per unit volume
    double*                                      stateVars;            ///< pointer to state variables
  };

  /**
   * @struct DeformationAD
   * @brief Represents the deformation state of a material.
   *
   * This struct holds the deformation gradient \f$\boldsymbol{F}\f$ which describes the local deformation of a
   * material.
   */
  template < int nDim >
  struct DeformationAD {
    Fastor::Tensor< autodiff::dual, nDim, nDim > F; ///< deformation gradient
  };
  void computeStress( ConstitutiveResponse< 3 >& response,
                      AlgorithmicModuli< 3 >&    tangents,
                      const Deformation< 3 >&    deformation,
                      const TimeIncrement&       timeIncrement ) const override
  {
    using namespace Marmot; // Assuming Marmot::AutomaticDifferentiation is in this namespace
    using namespace FastorStandardTensors;
    using scalar = autodiff::dual;

    Eigen::Map< Eigen::VectorXd > stateVars( response.stateVars, this->getNumberOfRequiredStateVars() );

    // Remember old state to reset during potential multiple AD evaluations
    const Eigen::VectorXd stateVarsOld = stateVars;
    const Tensor33d       tauOld       = response.tau;

    std::function< Tensor33t< scalar >( const Tensor33t< scalar >& ) > computeTauAD =
      [&]( const Tensor33t< scalar >& F_ ) {
        // Reset stateVars to old state
        stateVars = stateVarsOld;

        // Explicitly convert old (double-valued) stress to a dual-valued tensor for AD
        Tensor33t< scalar > tauAD = makeDual( tauOld );

        // Construct AD state
        DeformationAD< 3 > deformationAD;
        deformationAD.F = F_;

        ConstitutiveResponseAD< 3 > responseAD;
        responseAD.tau                  = tauAD;
        responseAD.rho                  = response.rho;
        responseAD.elasticEnergyDensity = response.elasticEnergyDensity;
        responseAD.stateVars            = stateVars.data();

        // Compute stress utilizing the child class's AD implementation
        this->computeStressAD( responseAD, deformationAD, timeIncrement );

        // Extract and return the dual stress tensor for differentiation
        Tensor33t< scalar > res = responseAD.tau;

        return res;
      };

    // Compute Kirchhoff stress (tau) and algorithmic tangent (dTau_dF) with autodiff
    std::tie( response.tau, tangents.dTau_dF ) = Marmot::AutomaticDifferentiation::dF_dT( computeTauAD, deformation.F );
  }

  /**
   * @brief AD-based constitutive update for finite strain materials.
   *
   * This method is called repeatedly from computeStress() during seeded
   * automatic-differentiation evaluations performed by
   * Marmot::AutomaticDifferentiation::dF_dT. For each seeded evaluation,
   * the state variables in response.stateVars are reset to the same
   * "old" state before computeStressAD() is invoked.
   *
   * Implementations MUST update response.stateVars based only on the
   * primal values of the dual numbers passed in (e.g. deformation.F,
   * response.tau, etc.), and MUST NOT depend on or modify state based
   * on AD gradient components. After dF_dT completes, response.stateVars
   * will contain the state corresponding to the last seeded evaluation,
   * which is correct only under this assumption.
   */
  virtual void computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                const DeformationAD< 3 >&    deformation,
                                const TimeIncrement&         timeIncrement ) const = 0;
};
