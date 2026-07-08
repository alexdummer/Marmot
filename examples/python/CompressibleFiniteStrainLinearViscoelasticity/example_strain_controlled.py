import marmot
import numpy as np

print("Running example for FiniteStrain material: COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY")

# Material properties extracted from C++ tests
properties = np.array([0.0, 0.0, 3500.0, 1500.0, 1.0, 0.3, 10.0], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY", properties, options)

# Setup a loading step
step = marmot.solvers.FiniteStrainSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1
step.gradUIncrementTarget = np.array([[0.01, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64)

step.isGradUComponentControlled = np.array([[True, True, True], [True, True, True], [True, True, True]])
step.isStressComponentControlled = np.logical_not(step.isGradUComponentControlled)
step.stressIncrementTarget = np.zeros((3, 3), dtype=np.float64)

solver.addStep(step)
solver.solve()

history = solver.getHistory()
if len(history) == 0:
    raise RuntimeError("Solver history is empty!")

import os

reference_file = os.path.splitext(__file__)[0] + "_reference.npz"

stress = np.array([h.stress for h in history])
if hasattr(history[0], "F"):
    strain = np.array([h.F for h in history])
else:
    strain = np.array([h.strain for h in history])

if not os.path.exists(reference_file):
    np.savez(reference_file, stress=stress, strain=strain)
    print(f"Reference solution created: {reference_file}")
else:
    ref = np.load(reference_file)
    np.testing.assert_allclose(stress, ref["stress"], atol=1e-6, rtol=1e-5)
    np.testing.assert_allclose(strain, ref["strain"], atol=1e-6, rtol=1e-5)
    print(f"Reference solution passed: {reference_file}")
