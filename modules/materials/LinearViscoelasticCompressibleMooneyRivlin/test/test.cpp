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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <string>

using namespace Marmot::Testing;
using namespace Marmot::Solvers;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

// -----------------------------------------------------------------------
// Material property helpers
// -----------------------------------------------------------------------

// Layout: [C1, C2, K, nMaxwell, (gamma1, tau1, ...)]

static std::vector< double > getElasticMooneyRivlinProps()
{
  // C1=500, C2=200, K=3000, nMaxwell=0
  return { 500.0, 200.0, 3000.0, 0.0 };
}

static std::vector< double > getViscoelasticMooneyRivlinProps()
{
  // C1=500, C2=200, K=3000, nMaxwell=1, gamma=0.3, tau=10
  return { 500.0, 200.0, 3000.0, 1.0, 0.3, 10.0 };
}

static MarmotMaterialPointSolverFiniteStrain makeSolver( const std::string&     matName,
                                                         std::vector< double >& matProps )
{
  auto solveropts = MarmotMaterialPointSolverFiniteStrain::SolverOptions();
  return MarmotMaterialPointSolverFiniteStrain( matName,
                                               matProps.data(),
                                               static_cast< int >( matProps.size() ),
                                               solveropts );
}

static MarmotMaterialPointSolverFiniteStrain::Step makeStep( const Tensor33d& gradUIncrement,
                                                             double           timeStart,
                                                             double           timeEnd,
                                                             double           dT )
{
  MarmotMaterialPointSolverFiniteStrain::Step step;
  step.gradUIncrementTarget        = gradUIncrement;
  step.stressIncrementTarget       = Tensor33d( 0.0 );
  step.isGradUComponentControlled  = Tensor33t< bool >( true );
  step.isStressComponentControlled = Tensor33t< bool >( false );
  step.timeStart                   = timeStart;
  step.timeEnd                     = timeEnd;
  step.dTStart                     = dT;
  step.dTMax                       = dT;
  return step;
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

// Test I-1: F=I gives zero Kirchhoff stress (elastic Mooney–Rivlin)
void testUndeformedResponse()
{
  const std::string matName  = "LINEARVISCOELASTICCOMPRESSIBLEMOONEYRIVLIN";
  auto              matProps = getElasticMooneyRivlinProps();
  auto              solver   = makeSolver( matName, matProps );

  solver.addStep( makeStep( Tensor33d( 0.0 ), 0.0, 1.0, 1.0 ) );
  solver.solve();

  throwExceptionOnFailure(
    checkIfEqual( solver.getHistory().back().stress, Tensor33d( 0.0 ), 1e-10 ),
    "I-1: Undeformed configuration - stress should be zero in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-2: Uniaxial stretch F_11=1.1 matches Mooney–Rivlin reference (elastic, nMaxwell=0)
void testUniaxialElasticResponse()
{
  const std::string matName  = "LINEARVISCOELASTICCOMPRESSIBLEMOONEYRIVLIN";
  auto              matProps = getElasticMooneyRivlinProps();
  auto              solver   = makeSolver( matName, matProps );

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.1; // F_11 = 1.1
  solver.addStep( makeStep( gradU, 0.0, 1.0, 1.0 ) );
  solver.solve();

  auto finalStress = solver.getHistory().back().stress;

  // Reference: Mooney-Rivlin C1=500, C2=200, K=3000 (D1=2/K), F=diag(1.1,1,1)
  // Computed via numerical differentiation of psi(C) = C1*(I1bar-3) + C2*(I2bar-3) + 1/D1*(0.5*(J^2-1)-lnJ)
  // tau_11 = 495.698, tau_22 = tau_33 = 224.651
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 495.698;
  stressTarget( 1, 1 ) = 224.651;
  stressTarget( 2, 2 ) = 224.651;

  throwExceptionOnFailure(
    checkIfEqual( finalStress, stressTarget, 1e-2 ),
    "I-2: Uniaxial elastic response failed in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-3: Kirchhoff stress tensor is symmetric for arbitrary deformation
void testStressTensorSymmetry()
{
  const std::string matName  = "LINEARVISCOELASTICCOMPRESSIBLEMOONEYRIVLIN";
  auto              matProps = getElasticMooneyRivlinProps();
  auto              solver   = makeSolver( matName, matProps );

  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.01;
  gradU( 0, 1 ) = 0.06;
  gradU( 0, 2 ) = -0.03;
  gradU( 1, 0 ) = 0.06;
  gradU( 1, 1 ) = 0.02;
  gradU( 1, 2 ) = 0.04;
  gradU( 2, 0 ) = -0.03;
  gradU( 2, 1 ) = 0.04;
  gradU( 2, 2 ) = -0.05;
  solver.addStep( makeStep( gradU, 0.0, 1.0, 1.0 ) );
  solver.solve();

  auto tau = solver.getHistory().back().stress;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      throwExceptionOnFailure(
        checkIfEqual( tau( i, j ), tau( j, i ), 1e-10 ),
        "I-3: Stress symmetry failed for tau(" + std::to_string( i ) + "," + std::to_string( j ) + ") in " +
          std::string( __PRETTY_FUNCTION__ ) );
}

// Test I-4: Pure rotation gives zero Kirchhoff stress
void testPureRotationZeroStress()
{
  const std::string matName  = "LINEARVISCOELASTICCOMPRESSIBLEMOONEYRIVLIN";
  auto              matProps = getElasticMooneyRivlinProps();

  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d F( 0.0 );
    F( 0, 0 ) = cos( phi );
    F( 0, 1 ) = -sin( phi );
    F( 1, 0 ) = sin( phi );
    F( 1, 1 ) = cos( phi );
    F( 2, 2 ) = 1.0;

    auto solver = makeSolver( matName, matProps );
    solver.addStep( makeStep( F - Spatial3D::I, 0.0, 1.0, 1.0 ) );
    solver.solve();

    throwExceptionOnFailure(
      checkIfEqual( solver.getHistory().back().stress, Tensor33d( 0.0 ), 1e-10 ),
      "I-4: Pure rotation (phi=" + std::to_string( phi_deg ) + ") should give zero stress in " +
        std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test I-5: Objectivity: tau(Q*F) = Q * tau(F) * Q^T
void testObjectivity()
{
  const std::string matName  = "LINEARVISCOELASTICCOMPRESSIBLEMOONEYRIVLIN";
  auto              matProps = getElasticMooneyRivlinProps();

  Tensor33d F( 0.0 );
  F( 0, 0 ) = 1.01;
  F( 0, 1 ) = 0.06;
  F( 0, 2 ) = -0.03;
  F( 1, 0 ) = 0.06;
  F( 1, 1 ) = 1.02;
  F( 1, 2 ) = 0.04;
  F( 2, 0 ) = -0.03;
  F( 2, 1 ) = 0.04;
  F( 2, 2 ) = 0.95;

  auto solverRef = makeSolver( matName, matProps );
  solverRef.addStep( makeStep( F - Spatial3D::I, 0.0, 1.0, 1.0 ) );
  solverRef.solve();
  Tensor33d tauRef = solverRef.getHistory().back().stress;

  for ( int phi_deg = 30; phi_deg <= 180; phi_deg += 30 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );

    Tensor33d Q( 0.0 );
    Q( 0, 0 ) = cos( phi );
    Q( 0, 1 ) = -sin( phi );
    Q( 1, 0 ) = sin( phi );
    Q( 1, 1 ) = cos( phi );
    Q( 2, 2 ) = 1.0;

    Tensor33d F_rot  = einsum< ik, kj, to_ij >( Q, F );
    auto      solver = makeSolver( matName, matProps );
    solver.addStep( makeStep( F_rot - Spatial3D::I, 0.0, 1.0, 1.0 ) );
    solver.solve();
    Tensor33d tauRot = solver.getHistory().back().stress;

    // Expected: tau_rot = Q * tau_ref * Q^T
    Tensor33d tauExpected = einsum< iI, IJ, jJ, to_ij >( Q, tauRef, Q );

    throwExceptionOnFailure(
      checkIfEqual( tauRot, tauExpected, 1e-8 ),
      "I-5: Objectivity failed for phi=" + std::to_string( phi_deg ) + " in " +
        std::string( __PRETTY_FUNCTION__ ) );
  }
}

// Test I-6: Viscoelastic relaxation — only the deviatoric part relaxes.
// Long-term stress = tau_vol + (1-gamma)*tau_dev
void testViscoelasticRelaxation()
{
  const std::string matName  = "LINEARVISCOELASTICCOMPRESSIBLEMOONEYRIVLIN";
  auto              matProps = getViscoelasticMooneyRivlinProps();
  auto              solver   = makeSolver( matName, matProps );

  // Step 1: Apply instantaneous deformation (tiny time)
  Tensor33d gradU( 0.0 );
  gradU( 0, 0 ) = 0.1; // F_11 = 1.1
  solver.addStep( makeStep( gradU, 0.0, 1e-4, 1e-4 ) );

  // Step 2: Hold deformation for 1000 >> tau=10 (fully relaxed)
  solver.addStep( makeStep( Tensor33d( 0.0 ), 1e-4, 1000.0, 100.0 ) );

  solver.solve();

  auto stressRelaxed = solver.getHistory().back().stress;

  // Reference: Mooney-Rivlin C1=500, C2=200, K=3000, F=diag(1.1,1,1)
  // Elastic: tau_11=495.698, tau_22=224.651
  // Deviatoric split: PK2 = diag(409.668, 224.651, 224.651)
  //   PK2vol = trace/3 = 858.97/3 = 286.323
  //   PK2dev_11 = 123.345, PK2dev_22 = -61.672
  // After full relaxation (gamma=0.3):
  //   PK2_relax = PK2vol + (1-gamma)*PK2dev
  //   tau_relax_11 = 1.21 * PK2_relax_11 = 450.924
  //   tau_relax_22 = PK2_relax_22 = 243.153
  Tensor33d stressTarget( 0.0 );
  stressTarget( 0, 0 ) = 450.924;
  stressTarget( 1, 1 ) = 243.153;
  stressTarget( 2, 2 ) = 243.153;

  throwExceptionOnFailure(
    checkIfEqual( stressRelaxed, stressTarget, 1.0 ),
    "I-6: Viscoelastic relaxation - long-term stress wrong in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testUndeformedResponse,       // I-1: F=I gives zero stress
    testUniaxialElasticResponse,  // I-2: Uniaxial elastic stretch reference values
    testStressTensorSymmetry,     // I-3: Kirchhoff stress symmetry
    testPureRotationZeroStress,   // I-4: Pure rotation gives zero stress
    testObjectivity,              // I-5: Objectivity tau(Q*F) = Q*tau(F)*Q^T
    testViscoelasticRelaxation,   // I-6: Viscoelastic relaxation (deviatoric only)
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
