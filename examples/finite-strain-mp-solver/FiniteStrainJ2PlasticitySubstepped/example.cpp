#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialPointSolverFiniteStrain.h"
// Ensure the registration unit is linked or included if using static registration
// #include "Marmot/MaterialRegistration.cpp"

int main()
{
  using namespace Marmot::Solvers;
  using namespace Marmot::FastorStandardTensors;

  // 1) Define the material model
  // Must match the key used in MarmotMaterialFiniteStrainFactory::registerMaterial
  std::string materialName = "FINITESTRAINJ2PLASTICITY_SUBSTEPPED";

  // Define Properties
  // [0]: nSubsteps (Specific to SubsteppingMaterial)
  // [1...8]: Standard FiniteStrainJ2Plasticity properties
  double properties[] = {
    20,   // nSubsteps (20 substeps per global increment)
    35,   // K (Bulk Modulus)
    15,   // G (Shear Modulus)
    0.45, // fy (Initial Yield Stress)
    0.75, // fyInf (Saturated Yield Stress)
    0.5,  // eta (Saturation rate)
    0.1,  // H (Linear Hardening)
    1.0,  // implementationType (1 = Full Analytical Return Mapping)
    0.0   // density
  };
  int nProps = 9;

  // 2) Configure solver options
  // With a consistent analytical tangent from the substepper,
  // we expect quadratic convergence (typically < 5 iterations).
  MarmotMaterialPointSolverFiniteStrain::SolverOptions options;
  options.maxIterations       = 15;
  options.residualTolerance   = 1e-8;
  options.correctionTolerance = 1e-8;

  // 3) Create solver instance
  // The factory will instantiate SubsteppingMaterial, which internally instantiates FiniteStrainJ2Plasticity
  MarmotMaterialPointSolverFiniteStrain solver( materialName, properties, nProps, options );

  // 5) Define Loading Steps

  // --- Step 1: Loading (Uniaxial extension to 10%) ---
  MarmotMaterialPointSolverFiniteStrain::Step step1;
  step1.timeStart = 0.0;
  step1.timeEnd   = 1.0;
  step1.dTStart   = 1e-0; // 10 global increments
  step1.dTMin     = 1e-1;
  step1.dTMax     = 0.5;

  step1.gradUIncrementTarget        = Tensor33d( 0.0 );
  step1.stressIncrementTarget       = Tensor33d( 0.0 );
  step1.isGradUComponentControlled  = Tensor33t< bool >( false );
  step1.isStressComponentControlled = Tensor33t< bool >( true );

  // Control gradU_11 (stretch to 10%)
  step1.gradUIncrementTarget( 0, 0 )        = 1;
  step1.isGradUComponentControlled( 0, 0 )  = true;
  step1.isStressComponentControlled( 0, 0 ) = false;

  // Mixed control for lateral contraction (Plane Stress approximation via Solver)
  // We keep stress 22 and 33 controlled to 0.0 (default true)

  solver.addStep( step1 );

  // --- Step 2: Unloading (Return to 0% global strain) ---
  // This validates that F_n is correctly updated and stored between steps
  MarmotMaterialPointSolverFiniteStrain::Step step2 = step1;
  step2.timeStart                                   = 1.0;
  step2.timeEnd                                     = 2.0;

  // Reverse direction: target increment is -0.10 to go back to 0
  step2.gradUIncrementTarget( 0, 0 ) = -1;

  /* solver.addStep( step2 ); */

  // 6) Run the solver
  try {
    solver.solve();
  }
  catch ( const std::exception& e ) {
    // Catch exceptions thrown by the substepper (e.g., if inner Newton fails)
    std::cerr << "\n!!! SOLVER FAILED !!!\nError: " << e.what() << std::endl;
    return 1;
  }

  // 7) Output the results
  // The CSV will contain extra columns for the substepper state (Substepping_F_n components)
  solver.exportHistoryToCSV( "j2_substepped_results.csv" );

  return 0;
}
