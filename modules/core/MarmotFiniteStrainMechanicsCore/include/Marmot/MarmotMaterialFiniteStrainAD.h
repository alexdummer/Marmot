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
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"

class MarmotMaterialFiniteStrainAD : public MarmotMaterialFiniteStrain {

public:
  MarmotMaterialFiniteStrainAD( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
    : MarmotMaterialFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ )
  {
  }
  /**
   * @struct ConstitutiveResponse
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
   * @struct Deformation
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

        // Promote old stress to dual (Fastor usually handles this implicit cast natively)
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

  virtual void computeStressAD( ConstitutiveResponseAD< 3 >& response,
                                const DeformationAD< 3 >&    deformation,
                                const TimeIncrement&         timeIncrement ) const = 0;
};
