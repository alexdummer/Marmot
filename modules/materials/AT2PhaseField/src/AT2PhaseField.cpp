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

#include "Marmot/AT2PhaseField.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotPhaseFieldEnergyDegradation.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  AT2PhaseField::AT2PhaseField( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 >( materialProperties, nMaterialProperties, materialNumber )
  {
    initializeStateLayout();
  }

  void AT2PhaseField::computeStress( response& res, tangents& tan, const increment& inc ) const
  {
    // material properties
    const double& E  = materialProperties[0];
    const double& nu = materialProperties[1];
    const double& Gc = materialProperties[2];
    const double& l  = materialProperties[3];

    // precompute elastic stiffness
    const Matrix6d C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    // access state variables
    double& H      = stateLayout.getAs< double& >( res.stateVars, "maxCrackDrivingForce" );
    auto    strain = stateLayout.getAs< Eigen::Map< Eigen::Matrix< double, 6, 1 > > >( res.stateVars, "strain" );

    // accumulate total strain
    const Vector6d eps = strain + inc.dStrain;

    // phase-field value at the current Gauss point (nonlocal variable K)
    const double phi = inc.K( 0 );

    // quadratic degradation function: g(phi) = (1-phi)^2
    const auto [g, dg_dphi, d2g_dphi2] = PhaseField::EnergyDegradationFunctions::SecondOrderDerived::quadratic( phi );

    // positive elastic strain energy density (no tension-compression split)
    const double psiPlus = 0.5 * eps.dot( C * eps );

    // irreversibility: update crack driving force history variable
    const bool   loading = psiPlus > H;
    const double H_new   = loading ? psiPlus : H;

    // compute degraded stress
    res.stress = g * ( C * eps );

    // local crack driving force (right-hand side of the phase-field equation)
    //   phi - l^2 * Delta(phi) = KLocal = 2*l/Gc * (1-phi) * H
    res.KLocal( 0 ) = ( 2. * l / Gc ) * ( 1. - phi ) * H_new;

    // gradient enhancement coefficient c = l^2
    res.c( 0 ) = l * l;

    // --- tangent moduli ---

    // d(stress) / d(dStrain) = g(phi) * C
    tan.dStressddStrain = g * C;

    // d(stress) / d(phi) = g'(phi) * C * eps  (stored as 6x1 column)
    tan.dStressddK.col( 0 ) = dg_dphi * ( C * eps );

    // d(KLocal) / d(dStrain):  2l/Gc * (1-phi) * dH/d(eps)
    //   Active loading: dH/d(eps) = d(psi+)/d(eps) = C * eps
    //   Unloading / elastic: dH/d(eps) = 0
    if ( loading ) {
      tan.dKLocalddStrain.row( 0 ) = ( 2. * l / Gc ) * ( 1. - phi ) * ( C * eps ).transpose();
    }
    // else: already zero (default-initialized in the tangents struct)

    // d(KLocal) / d(phi) = -2*l/Gc * H_new
    tan.dKLocalddK( 0, 0 ) = -( 2. * l / Gc ) * H_new;

    // dcddK and d2cddK2 remain zero (c = l^2 is constant)

    // update state variables
    strain = eps;
    H      = H_new;
  }

} // namespace Marmot::Materials
