import marmot
import numpy as np

print("Running example for FiniteStrain material: FINITESTRAINJ2PLASTICITY_SUBSTEPPED")

# Material properties extracted from C++ tests
properties = np.array([10.0, 175000.0, 80800.0, 260.0, 580.0, 9.0, 70.0, 1.0], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("FINITESTRAINJ2PLASTICITY_SUBSTEPPED", properties, options)

# Setup a loading step
step = marmot.solvers.FiniteStrainSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1
step.gradUIncrementTarget = np.zeros((3, 3), dtype=np.float64)
step.isGradUComponentControlled = np.array([[False, True, True], [True, False, True], [True, True, False]])
step.isStressComponentControlled = np.logical_not(step.isGradUComponentControlled)
s = np.zeros((3, 3), dtype=np.float64)
s[0, 0] = 10.0
step.stressIncrementTarget = s

solver.addStep(step)
solver.solve()

history = solver.getHistory()

import os
import marmot

reference_file = os.path.splitext(__file__)[0] + "_reference.npz"
marmot.compareHistoryToReferenceSolution(history, reference_file)
