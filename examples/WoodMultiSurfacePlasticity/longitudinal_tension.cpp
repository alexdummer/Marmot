#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"

int main()
{

  std::string materialName = "WOODMULTISURFACEPLASTICITY";

  std::vector< double > materialProperties;
  materialProperties.push_back( 9000 );   // E1 - Youngs modulus in direction 1
  materialProperties.push_back( 4500.0 ); // E2 - Youngs modulus in direction 2
  materialProperties.push_back( 130 );    // E3 - Youngs modulus in direction 3

  materialProperties.push_back( 0.21 );   // NU21 - Poisson's ratio
  materialProperties.push_back( 0.0036 ); // NU13 - Poisson's ratio
  materialProperties.push_back( 0.009 );  // NU23 - Poisson's ratio

  materialProperties.push_back( 1112 );   // G12 - Shear modulus
  materialProperties.push_back( 160 );    // G13 - Shear modulus
  materialProperties.push_back( 100 );    // G23 - Shear modulus

  // direction R
  materialProperties.push_back( 1.0 ); // aR1 - direction R vector component 1
  materialProperties.push_back( 0.0 ); // aR2 - direction R vector component 2
  materialProperties.push_back( 0.0 ); // aR3 - direction R vector component 3

  // direction T
  materialProperties.push_back( 0.0 ); // aT1 - direction T vector component 1
  materialProperties.push_back( 1.0 ); // aT2 - direction T vector component 2
  materialProperties.push_back( 0.0 ); // aT3 - direction T vector component 3

  // strength parameters
  materialProperties.push_back( 21.8 ); // ftR - tensile strength in direction R
  materialProperties.push_back( 26.3 ); // fcR - compressive strength in direction R
  materialProperties.push_back( 7.0 );  // ftT - tensile strength in direction T
  materialProperties.push_back( 9.8 );  // fcT - compressive strength in direction T
  materialProperties.push_back( 0.42 ); // ftL - tensile strength in direction L
  materialProperties.push_back( 0.42 ); // fcL - compressive strength in direction L

  materialProperties.push_back( 7.5 );  // fsRT - shear strength in plane RT
  materialProperties.push_back( 1.8 );  // fsRL - shear strength in plane RL
  materialProperties.push_back( 1.4 );  // fsTL - shear strength in plane TL

  // plastic hardening parameters
  materialProperties.push_back( 1. );  // kappa - radial tension hardening parameter
  materialProperties.push_back( 0.5 ); // eta - radial tension hardening parameter
  materialProperties.push_back( 0.0 ); // zeta - radial tension hardening parameter
  materialProperties.push_back( 1.0 ); // alphaD - radial tension hardening parameter
  materialProperties.push_back( 1.0 ); // alphaMax - radial tension hardening parameter

  materialProperties.push_back( 1. );  // kappa - radial compression hardening parameter
  materialProperties.push_back( 0.5 ); // eta - radial compression hardening parameter
  materialProperties.push_back( 0.0 ); // zeta - radial compression hardening parameter
  materialProperties.push_back( 1.0 ); // alphaD - radial compression hardening parameter
  materialProperties.push_back( 1.0 ); // alphaMax - radial compression hardening parameter

  materialProperties.push_back( 1. );  // kappa - tangential tension hardening parameter
  materialProperties.push_back( 0.5 ); // eta - tangential tension hardening parameter
  materialProperties.push_back( 0.0 ); // zeta - tangential tension hardening parameter
  materialProperties.push_back( 1.0 ); // alphaD - tangential tension hardening parameter
  materialProperties.push_back( 1.0 ); // alphaMax - tangential tension hardening parameter

  materialProperties.push_back( 1. );  // kappa - tangential compression hardening parameter
  materialProperties.push_back( 0.5 ); // eta - tangential compression hardening parameter
  materialProperties.push_back( 0.0 ); // zeta - tangential compression hardening parameter
  materialProperties.push_back( 1.0 ); // alphaD - tangential compression hardening parameter
  materialProperties.push_back( 1.0 ); // alphaMax - tangential compression hardening parameter

  materialProperties.push_back( 1 );   // kappa - longitudinal tension hardening parameter
  materialProperties.push_back( 50 );  // eta - longitudinal tension hardening parameter
  materialProperties.push_back( 0.0 ); // zeta - longitudinal tension hardening parameter
  materialProperties.push_back( 1.0 ); // alphaD - longitudinal tension hardening parameter
  materialProperties.push_back( 1.0 ); // alphaMax - longitudinal tension hardening parameter

  materialProperties.push_back( 1 );   // kappa - longitudinal compression hardening parameter
  materialProperties.push_back( 0.5 ); // eta - longitudinal compression hardening parameter
  materialProperties.push_back( 0.0 ); // zeta - longitudinal compression hardening parameter
  materialProperties.push_back( 1.0 ); // alphaD - longitudinal compression hardening parameter
  materialProperties.push_back( 1.0 ); // alphaMax - longitudinal compression hardening parameter

  materialProperties.push_back( 1 );   // kappa - RT shear hardening parameter
  materialProperties.push_back( 0.5 ); // eta - RT shear hardening parameter
  materialProperties.push_back( 0.0 ); // zeta - RT shear hardening parameter
  materialProperties.push_back( 1.0 ); // alphaD - RT shear hardening parameter
  materialProperties.push_back( 1.0 ); // alphaMax - RT shear hardening parameter

  int nMaterialProperties = materialProperties.size();

  std::cout << "Number of material properties: " << nMaterialProperties << std::endl;

  auto options = MarmotMaterialPointSolverHypoElastic::SolverOptions();

  MarmotMaterialPointSolverHypoElastic solver( materialName, &materialProperties[0], nMaterialProperties, options );

  // define a step
  MarmotMaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, false, false, false, false, false };
  step.isStressComponentControlled = { false, true, true, true, true, true };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.01, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.dTStart                     = .05;

  // add step to solver
  solver.addStep( step );

  solver.solve();

  solver.exportHistoryToCSV( "longitudinal_tension_history.csv" );

  return 0;
}
