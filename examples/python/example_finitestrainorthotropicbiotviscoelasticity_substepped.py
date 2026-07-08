import marmot
import numpy as np

print("Running example for FiniteStrain material: FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY_SUBSTEPPED")

# Material properties extracted from C++ tests
properties = np.array([1.0, 0.3], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY_SUBSTEPPED", properties, options)

# Setup a loading step
step = marmot.solvers.FiniteStrainSolver.Step()
step.gradUIncrementTarget = np.array([[0.01, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64)

step.isGradUComponentControlled = np.array([[True, True, True], [True, True, True], [True, True, True]])
step.isStressComponentControlled = np.logical_not(step.isGradUComponentControlled)
step.stressIncrementTarget = np.zeros((3, 3), dtype=np.float64)

solver.addStep(step)
solver.solve()

history = solver.getHistory()
if len(history) > 0:
    print(f"Final Stress XX: {history[-1].stress[0,0]}")
