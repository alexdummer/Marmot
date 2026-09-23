#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"

int main()
{

  std::string materialName = "TSAIWUPLASTICITY";

  std::vector< double > materialProperties;

  // Elements 0-8: Orthotropic elasticity
  materialProperties.push_back( 9000 );   // E1
  materialProperties.push_back( 4500.0 ); // E2
  materialProperties.push_back( 130 );    // E3
  materialProperties.push_back( 0.21 );   // NU12
  materialProperties.push_back( 0.009 );  // NU23
  materialProperties.push_back( 0.0036 ); // NU13
  materialProperties.push_back( 1112 );   // G12
  materialProperties.push_back( 100 );    // G23
  materialProperties.push_back( 160 );    // G13

  // Elements 9-17: Tsai-Wu strengths
  materialProperties.push_back( 21.8 ); // T1
  materialProperties.push_back( 26.3 ); // C1
  materialProperties.push_back( 7.0 );  // T2
  materialProperties.push_back( 9.8 );  // C2
  materialProperties.push_back( 0.42 ); // T3
  materialProperties.push_back( 0.42 ); // C3
  materialProperties.push_back( 7.5 );  // S12
  materialProperties.push_back( 1.4 );  // S23
  materialProperties.push_back( 1.8 );  // S13

  // Elements 18-20: Hardening parameters
  materialProperties.push_back( 0000.0 ); // HLin
  materialProperties.push_back( 0.0 );    // deltaYieldStress
  materialProperties.push_back( 00.0 );   // delta

  // Elements 21-26: coordinate systems
  materialProperties.push_back( 1.0 ); // direction1 x
  materialProperties.push_back( 0.0 ); // direction1 y
  materialProperties.push_back( 0.0 ); // direction1 z
  materialProperties.push_back( 0.0 ); // direction2 x
  materialProperties.push_back( 0.0 ); // direction2 y
  materialProperties.push_back( 1.0 ); // direction2 z

  // Elements 27-28: Density and Damping Coefficient
  materialProperties.push_back( 1.5e-9 ); // density
  materialProperties.push_back( 0.0 );    // damping coefficient

  int nMaterialProperties = materialProperties.size();

  std::cout << "Number of material properties: " << nMaterialProperties << std::endl;

  using namespace Marmot::Solvers;
  auto options = MarmotMaterialPointSolverHypoElastic::SolverOptions();

  MarmotMaterialPointSolverHypoElastic solver( materialName, &materialProperties[0], nMaterialProperties, options );

  // define a step
  MarmotMaterialPointSolverHypoElastic::Step step;

  // Tension in direction 1
  step.isStrainComponentControlled = { false, true, false, false, false, false };
  step.isStressComponentControlled = { true, false, true, true, true, true };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.01, 0.02, -0.002, 0.0, 0.0, 0.0 };
  step.dTStart                     = .05;
  solver.addStep( step );

  step.isStrainComponentControlled = { false, true, false, false, false, false };
  step.isStressComponentControlled = { true, false, true, true, true, true };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { -0.02, -0.04, -0.002, 0.0, 0.0, 0.0 };
  step.dTStart                     = .05;
  solver.addStep( step );
  /* // Shear in plane 1-2 */
  /* step.isStrainComponentControlled = { true, true, true, true, true, true }; */
  /* step.isStressComponentControlled = { false, false, false, false, false, false }; */
  /* step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; */
  /* step.strainIncrementTarget       = { -0.01, 0.002, 0.002, 0.01, 0.0, 0.0 }; */
  /* step.dTStart                     = .05; */
  /* solver.addStep( step ); */

  solver.solve();

  solver.exportHistoryToCSV( "tsai_wu_plasticity_history.csv" );

  return 0;
}
