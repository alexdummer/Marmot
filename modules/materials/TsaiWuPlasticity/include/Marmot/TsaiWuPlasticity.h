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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"

namespace Marmot::Materials {

  /// @brief An implementation of Tsai-Wu Plasticity model with isotropic hardening and orthotropic elasticity.
  class TsaiWuPlasticityModel : public MarmotMaterialHypoElastic {

    /// @brief Elastic properties
    const double& E1;
    const double& E2;
    const double& E3;
    const double& nu12;
    const double& nu23;
    const double& nu13;
    const double& G12;
    const double& G23;
    const double& G13;

    /// @brief Strength parameters
    const double& T1;
    const double& C1;
    const double& T2;
    const double& C2;
    const double& T3;
    const double& C3;
    const double& S12;
    const double& S23;
    const double& S13;

    /// @brief Hardening parameters
    const double& HLin;
    const double& deltaYieldStress;
    const double& delta;

    /// @brief Material coordinate system directions
    const Eigen::Vector3d direction1;
    const Eigen::Vector3d direction2;

  public:
    TsaiWuPlasticityModel( const double* materialProperties, const int nMaterialProperties, const int materialLabel );

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStressDDStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

    /**
     * @brief Get material density.
     * @return Density value.
     * @throw std::runtime_error if density is not defined.
     */
    double getDensity( const double* stateVars ) const override;

    void initializeStateLayout()
    {
      stateLayout.add( "kappa", 1 );
      stateLayout.finalize();
    }

  private:
    /// @brief Local coordinate system of the material
    Eigen::Matrix3d localCoordinateSystem;

    /// @brief Initial elastic stiffness in local coordinate system
    Eigen::Matrix< double, 6, 6 > Cel;

    /// @brief Initial elastic stiffness in global coordinate system
    Eigen::Matrix< double, 6, 6 > CelGlobal;

    /// @brief Tsai-Wu coefficient matrix P
    Eigen::Matrix< double, 6, 6 > P;

    /// @brief Tsai-Wu coefficient vector q
    Eigen::Matrix< double, 6, 1 > q;
  };

} // namespace Marmot::Materials
